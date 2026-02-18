# Automated Solvent Box Generation and Equilibration

This toolkit builds and equilibrates GROMACS solvent boxes for arbitrary solutes and solvents, including ion placement and local ion–solute rearrangement.

Scripts:

- `solvent.sh` – main driver: builds solvated, ionized box and prepares MD inputs.
- `swap_ion_solvent.py` – Python helper: moves ions from bulk into the solute environment.
- `solvent_run.sh` – runs NVT/NPT equilibration and basic analysis (RMSD, density).


## 1. Requirements

- GROMACS 2019.4 (or compatible) available as `gmx` and `mdrun_gpu_mpi`.[^1][^2]
- Bash, `awk`, `sed`, `grep`.
- Python 3 + MDAnalysis + NumPy for `swap_ion_solvent.py`.[^3]
- Cluster launcher `mBASH` (or adjust `solvent.sh` to your scheduler).[^2]

Environment modules and MPI settings in the scripts are cluster-specific; adapt the `module add ...` and `OMPI_MCA_*` lines as needed.[^1][^2]

## 2. Input files

In each project directory (e.g. `0a/`, `1a/`), you must provide:

- `Infos.dat` – steering file with key–value pairs, e.g.:[^2]

```text
Project 0a
Solvent THF
Resname UNK
Boxsize 6.5
Ions -neutral
Fixed UNK
Posre 50000
Anchor UNK:CA
Temperature 293
NVT 5000000
NPT 5000000
```

    - `Project` – basename of the solute PDB (`Project.pdb`).
    - `Solvent` – solvent residue name and `.gro` prefix (expects `Solvent.gro` in the directory, e.g. `THF.gro`).[^2]
    - `Resname` – solute residue name used in `conf_UNK.itp` etc. (not directly used in the scripts yet, but kept for clarity).[^2]
    - `Boxsize` – cubic box edge in nm.
    - `Ions` – arguments passed to `gmx genion` (e.g. `-neutral`, `-np 10 -nn 10`, etc.).[^2]
    - `Fixed` – residue name whose atoms are position-restrained and frozen during minimization (e.g. `UNK`).[^2]
    - `Posre` – force constant for position restraints (kJ mol$^{-1}$ nm$^{-2}$).[^2]
    - `Anchor` – solute selector for ion swapping, either `RESNAME` or `RESNAME:ATOM` (e.g. `UNK` or `UNK:CA`).[^3][^2]
    - `Temperature` – target temperature (K).
    - `NVT`, `NPT` – number of MD steps for NVT / NPT.[^2]

Additional required files in the same directory:

- `<Project>.pdb` – solute structure.[^2]
- `<Solvent>.gro` – pre-equilibrated pure solvent box used as `-cs` for `gmx solvate`.[^2]
- `topol.top`, `conf_UNK.itp`, solvent `.itp` files, and force field directory as usual for a GROMACS system.[^4][^5][^2]

The script also expects templates in a shared script directory (configured as `sd` in `solvent.sh`):[^2]

- `minim.mdp`, `nvt.mdp`, `npt.mdp`
- `swap_ion_solvent.py`
- `solvent_run.sh`

You can copy these into your repo and adjust `sd=...` accordingly.

## 3. Workflow overview

`solvent.sh` performs the full workflow:[^2]

1. Read `Infos.dat`.
2. Build a cubic box around the solute.
3. Fill with chosen solvent (`gmx solvate`).
4. Add ions (`gmx genion`) with user-specified options.
5. Use `swap_ion_solvent.py` to move ions from bulk into the solute environment (within 3–5 Å).
6. Build index groups and position restraints.
7. Patch `.mdp` files for POSRES and correct temperature coupling groups.
8. Run energy minimization.
9. Submit NVT/NPT equilibration via `solvent_run.sh`.

`solvent_run.sh` then:[^1]

- Runs NVT and NPT MD.
- Computes RMSD vs time (`rmsd_nvt.xvg`, `rmsd_npt.xvg`).
- Extracts density vs time from NPT (`density_npt.xvg`).


## 4. Script details

### 4.1 `swap_ion_solvent.py`

Purpose: translate a set of ions from bulk to the solute vicinity by swapping their center-of-mass positions with solvent molecules near the solute.[^3]

Usage:

```bash
python3 swap_ion_solvent.py input.gro output.gro SOLUTE_SPEC SOLVENT_RES ION1 [ION2 ...]
```

- `input.gro` – solvated, ionized configuration (`solv_ions.gro`).[^2]
- `output.gro` – output after swapping.
- `SOLUTE_SPEC` – solute selection:
    - `UNK` → `resname UNK`
    - `UNK:CA` → `resname UNK and name CA` (anchors on a specific atom).[^3]
- `SOLVENT_RES` – solvent residue name (e.g. `THF` or `SOL`).[^3][^2]
- `ION1 [ION2 ...]` – ion residue names to move (e.g. `NA CL`).[^3][^2]

Algorithm:[^3]

- Load `input.gro` via MDAnalysis.
- Build selections:
    - Solute: `solute_sel` according to `SOLUTE_SPEC`.
    - Solvent: `resname SOLVENT_RES`.
    - Ions: `resname ION1 or resname ION2 ...`.
- Find solvent residues within `cutoff` Å (default 3.0; 5.0 in `solvent.sh` call) of the solute.
- For up to `nswap = min(#nearby solvent residues, #ion residues)`:
    - Compute COMs of a solvent residue and an ion residue.
    - Translate all atoms of the solvent residue to the ion COM, and vice versa.
- Write the modified structure to `output.gro`.[^3]

This effectively relocates ions close to the solute, while preserving solvent and ion COM distributions globally.

### 4.2 `solvent.sh`

```
Main pipeline script; run **inside a project directory** that contains `Infos.dat`, `<Project>.pdb`, `<Solvent>.gro`, and topology files.[^2]
```

Key steps:[^2]

1. **Read `Infos.dat`:**
    - Shell loop parses keys and values into variables (`Solvent`, `Project`, `Resname`, `Boxsize`, `Ions`, `Fixed`, `Posre`, `Anchor`, `Temperature`, `NVT`, `NPT`).[^2]
    - Missing mandatory keys cause an immediate error (`: "${Project:?...}"`).[^2]
    - If `Solvent` is empty, defaults to `Water`.[^2]
2. **Prepare starting structures:**
    - Convert `<Project>.pdb` to `solute.gro`.[^2]
    - Center in a cubic box of size `Boxsize` as `boxed.gro`.[^2]
3. **Solvate:**
    - `gmx solvate -cp boxed.gro -cs <Solvent>.gro -o solv.gro -p topol.top`.[^2]
4. **Add ions:**
    - Pre-GROMPP: `gmx grompp -f minim.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1`.[^2]
    - Run `gmx genion` with the user-specified `Ions` options and fixed ion names `-pname NA -nname CL`, replacing solvent molecules (you select the solvent group interactively via heredoc `${Solvent}`).[^2]
    - Result: `solv_ions.gro`, `topol.top` updated with ions.[^2]
5. **Ion relocation near solute:**
    - Copy `swap_ion_solvent.py` from `sd` and run:[^2]

```bash
python3 swap_ion_solvent.py solv_ions.gro solv_ions_swapped.gro "$Anchor" "$Solvent" "NA CL"
```

    - Replace `solv_ions.gro` with swapped version.[^2]
6. **Index and position restraints:**
    - Run `gmx make_ndx` once to count default groups and determine group indices.[^2]
    - Create:
        - `${Fixed}_posres`: all atoms of residue `Fixed`.[^2]
        - `non-${Solvent}`: everything that is not solvent.[^2]
    - Verify both groups exist in `index.ndx`.[^2]
    - Generate position restraints for `${Fixed}_posres`:

```bash
gmx genrestr -f solv_ions.gro -n index.ndx -o posre_${Fixed}.itp -fc $Posre $Posre $Posre
```

    - Ensure `topol.top` includes `posre_UNK.itp` after `conf_UNK.itp` via an `awk` patch if not already present.[^2]
7. **Prepare `.mdp` files:**
    - Copy `minim.mdp`, `nvt.mdp`, `npt.mdp` from `sd`.[^2]
    - Replace placeholders:
        - `TTT` → `Temperature` (K).
        - `XXX` → `NVT` / `NPT` step counts.[^2]
    - Add `-DPOSRES` to `define` lines (or create `define = -DPOSRES` if absent) in all three `.mdp`s.[^2]
    - Append `freezegrps = Fixed` and `freezedim = Y Y Y` to `minim.mdp` to freeze that residue during minimization.[^2]
    - In `nvt.mdp` and `npt.mdp`, adjust temperature-coupling groups:

```text
tc-grps = non-Solvent Solvent
```

where `Solvent` is taken from `Infos.dat` (`non-${Solvent} ${Solvent}`).[^2]
8. **Energy minimization:**
    - `gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -n index.ndx -r solv_ions.gro -o em.tpr -maxwarn 1`
    - `gmx mdrun -deffnm em`[^2]
9. **Equilibration launch:**
    - Copy `solvent_run.sh` from `sd` and submit via your cluster wrapper:[^2]

```bash
cp $sd/solvent_run.sh .
bash ~/bin/mBASH solvent_run -Nnod tornado -Npr 24 -Time 30
```


Adapt the last line to your own scheduler (e.g. `sbatch solvent_run.sh`).

### 4.3 `solvent_run.sh`

Runs NVT and NPT equilibration and basic analysis.[^1]

Steps:

1. Load GROMACS module and set MPI environment variables (cluster-specific).[^1]
2. NVT equilibration:[^1]

```bash
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -n index.ndx -r em.gro -maxwarn 1
mpiexec --bind-to none mdrun_gpu_mpi -deffnm nvt
```

3. NVT RMSD vs time:[^1]

```bash
gmx rms -s nvt.tpr -f nvt.xtc -o rmsd_nvt.xvg -tu ps << EOF
0
0
EOF
```

4. NPT equilibration:[^1]

```bash
gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -n index.ndx -r nvt.gro -maxwarn 1
mpiexec --bind-to none mdrun_gpu_mpi -deffnm npt
```

5. NPT RMSD vs time:[^1]

```bash
gmx rms -s npt.tpr -f npt.xtc -o rmsd_npt.xvg -tu ps << EOF
0
0
EOF
```

6. Density vs time from NPT:[^1]

```bash
gmx energy -f npt.edr -o density_npt.xvg << EOF
Density
EOF
```


Outputs:

- `nvt.tpr`, `nvt.xtc`, `nvt.edr`, `nvt.gro`, `rmsd_nvt.xvg`, `energy_nvt.xvg`
- `npt.tpr`, `npt.xtc`, `npt.edr`, `npt.gro`, `rmsd_npt.xvg`, `density_npt.xvg`[^1]


## 5. Typical usage

1. Prepare a project directory, e.g. `0a/`, with:
    - `Infos.dat`
    - `0a.pdb`
    - `THF.gro` (or other solvent)
    - `topol.top`, `conf_UNK.itp`, solvent `.itp`s, force field files
2. From inside that directory:

```bash
/path/to/solvent.sh
```

3. The script will:
    - Build and ionize the solvent box.
    - Move ions into the solute’s local shell.
    - Generate index and position restraints.
    - Run EM and then submit NVT/NPT via `solvent_run.sh`.
4. When equilibration finishes, inspect:
    - `rmsd_npt.xvg` – solute RMSD during NPT.
    - `density_npt.xvg` – box density vs time.

You can then use the final `npt.gro` / `npt.cpt` as starting points for production MD.

<div align="center">⁂</div>

[^1]: solvent_run.sh

[^2]: solvent.sh

[^3]: swap_ion_solvent.py

[^4]: conf_UNK.itp.txt

[^5]: THF.itp.txt

