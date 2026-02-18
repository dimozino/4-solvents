# Generation of snapshots for ASEC

This toolkit prepares and runs MD simulation starting from a representative snapshot to generate ASEC environment for solute. ASEC is the averaged solvent configuration generated from 100 snapshots along MD. 
Since the solute (e.g. the molecular system of interest that will be considered inside ASEC environment) can consists of several molecules, the script also merges all solute molecules into a single one (named UNK).  
The toolkit consists of:

- `asec_gensnapshots.sh` – main driver that prepares MD simulation. [file:311]
- `merge_itp.py` – merges multiple `.itp` files into a single UNK moleculetype. [file:313]
- `merge_3c.py` – merges selected residues into a single UNK residue in a GRO file and adjusts the topology. [file:312]
- `asec_mdrun.sh` – runs the NPT equilibration/production for the ASEC system. [file:310]

These scripts assume an existing equilibrated solvent box and an averaged-distance snapshot `<Project>_avg_dist.gro` produced earlier (e.g. via `extract_average.sh`).

---

## 1. `asec_gensnapshots.sh`

### Purpose

Automates preparation of MD run for ASEC snapshots generation:

1. Reads configuration from `Infos.dat`.
2. Copies the averaged snapshot and base topology into a new `asec/` directory.
3. Merges the specified `.itp` components into a single `conf_<Resname>.itp` using `merge_itp.py`. [file:311][file:313]
4. Merges the solute residues into a single UNK residue in the GRO/topology using `merge_3c.py`. [file:311][file:312]
5. Builds index groups and position restraints for the solute (very strong constraints, since we are sampling environment here). [file:311]
6. Prepares for and submits an ASEC NPT run via `asec_mdrun.sh`. [file:311][file:310]

### Expected inputs

Run this script in a project directory that contains: [file:311]

- `Infos.dat` – must include at least:

  ```text
  Project 0a
  Solvent THF
  Resname UNK
  Boxsize 6.5
  Temperature 293
  Ions -neutral
  Fixed UNK
  Posre 50000
  NVT 5000000
  NPT 5000000
  ASEC resname UNK or resid 1163 or resname CL
  ITP conf_UNK.itp THF.itp CL.itp
```

Fields used by `asec_gensnapshots.sh`:

- `Project` – basename (`Project`), used for `Project_avg_dist.gro`. [file:311]
- `Solvent` – solvent residue name; also expects `Solvent.itp` and `Solvent.gro` in the parent workflow. [file:311]
- `Resname` – solute residue name; affects `conf_${Resname}.itp`. [file:311]
- `Boxsize`, `Ions`, `Fixed`, `Posre`, `Temperature` – carried over for consistency/logging and restraint generation. [file:311]
- `ASEC` - VMD-style selection of solute residues that need to be merged. [file:311]
- `ITP` – space-separated list of `.itp` files to merge; first one must be `conf_${Resname}.itp`. [file:311][file:313]
- `<Project>_avg_dist.gro` – representative snapshot from the NPT trajectory (e.g. `0a_avg_dist.gro`). [file:311]
- `topol.top` – topology corresponding to this snapshot. [file:311]
- All `.itp` files listed in `ITP` in `Infos.dat`. [file:311]

The script also uses a script/template directory defined at the top: [file:311]

```bash
sd=/home/ibsbnicigt/nikolaev_d/script_library/4-solvents/
```

This directory must contain `merge_itp.py`, `merge_3c.py`, `npt_asec.mdp`, and `asec_mdrun.sh`. [file:311]

### Workflow

1. **Environment setup**: load GROMACS module, MPI env vars, and Conda environment. [file:311]
2. **Read `Infos.dat`** into shell variables: [file:311]
    - `Solvent`, `Project`, `Resname`, `Boxsize`, `Ions`, `Fixed`, `Posre`, `Temperature`, `ITP_FILES`.
    - Enforce that required variables are set; if `Solvent` is empty, default to `Water`. [file:311]
3. **Create and enter `asec/` directory**: [file:311]

```bash
mkdir -p asec
cd asec
```

4. **Copy base files**: [file:311]
    - `<Project>_avg_dist.gro` → working directory.
    - `topol.top`
    - `Infos.dat`
5. **Copy `.itp` files** listed under `ITP` in `Infos.dat`: [file:311]

```bash
for itp in $ITP_FILES; do
    cp ../"$itp" .
done
cp ../${Solvent}.itp .
```

6. **Merge `.itp` files** into a single `conf_${Resname}.itp` using `merge_itp.py`: [file:311][file:313]

```bash
python3 merge_itp.py $ITP_FILES
mv output.itp conf_${Resname}.itp
```

    - The first file in `ITP_FILES` must be `conf_${Resname}.itp`.
    - The merged file is always written as `output.itp` then renamed. [file:313]
7. **Merge ASEC residues and topology** using `merge_3c.py`: [file:311][file:312]

```bash
cp ${Project}_avg_dist.gro NA_avg_dist.gro
python3 merge_3c.py
mv topol_merged.top topol.top
mv NA_avg_dist_merged.gro ${Project}_avg_dist.gro
```

After this step, solute residues (defined in `Infos.dat` via the `ASEC` line) are merged into a single UNK residue in the GRO and topology. [file:312]
8. **Create index groups** for ASEC restraints and temperature coupling: [file:311]
    - First count default groups:

```bash
echo "q" | gmx make_ndx -f ${Project}_avg_dist.gro -o index_temp.ndx 2>&1 | tee make_ndx_output.txt
NGROUPS=$(grep -c "^ *[0-9]" make_ndx_output.txt || echo 0)
rm -f index_temp.ndx make_ndx_output.txt
```

    - Then create:

```bash
gmx make_ndx -f ${Project}_avg_dist.gro -o index.ndx << EOF
r ${Fixed}
!r ${Solvent}
name ${NGROUPS} ${Fixed}_posres
name $((NGROUPS + 1)) non-${Solvent}
q
EOF
```

This yields groups:
        - `${Fixed}_posres` – atoms of the ASEC complex to restrain.
        - `non-${Solvent}` – everything except solvent, used in `tc-grps` in `.mdp`. [file:311]
    - The script verifies that both groups exist in `index.ndx`. [file:311]
9. **Generate position restraints**: [file:311]

```bash
gmx genrestr -f ${Project}_avg_dist.gro -n index.ndx \
    -o posre_${Fixed}.itp -fc "$Posre" "$Posre" "$Posre" << EOF
${Fixed}_posres
EOF
```

10. **Prepare ASEC MD**: [file:311]

```bash
cp $sd/npt_asec.mdp .
cp $sd/asec_mdrun.sh .
bash ~/bin/mBASH asec_mdrun -Nnod tornado -Npr 24 -Time 30
```

The last line is cluster-specific; adapt to your scheduler (e.g. `sbatch asec_mdrun.sh`).

---

## 2. `merge_itp.py`

### Purpose

Merge several `.itp` files defining UNK and additional components (e.g. THF, CL, etc.) into one combined UNK moleculetype with correctly offset atom indices. [file:313]

The resulting combined file is always `output.itp` and is usually renamed to `conf_UNK.itp` by the calling script. [file:311][file:313]

### Usage

```bash
python3 merge_itp.py conf_UNK.itp [extra1.itp extra2.itp ...]
```

- `conf_UNK.itp` – base UNK moleculetype file; must contain `[ moleculetype ]` with the name `UNK`. [file:313]
- `extra*.itp` – any number of additional `.itp` files; for each, the first `[ moleculetype ]` block is merged. [file:313]

The script writes `output.itp` in-place. [file:313]

### Algorithm

1. **Read base UNK file** and extract its moleculetype block (`[ moleculetype ]` + sections until the next `[ moleculetype ]` or EOF). [file:313]
2. **Count atoms** in the UNK `[ atoms ]` section to determine initial offset. [file:313]
3. **Normalize UNK moleculetype header** so that the moleculetype name is `UNK` (even if originally different). [file:313]
4. **For each extra `.itp` file**: [file:313]
    - Find the first `[ moleculetype ]` block and infer its name. [file:313]
    - Extract its full block.
    - Count atoms in its `[ atoms ]` section.
    - Offset all atom indices (and references) by the cumulative atom count from previous blocks for the sections:
        - `[ atoms ]`, `[ bonds ]`, `[ angles ]`, `[ dihedrals ]`, `[ pairs ]`. [file:313]
    - Strip its `[ moleculetype ]` header and name line, leaving only the sections. [file:313]
    - Append to an internal list and update the cumulative atom count. [file:313]
5. **Assemble final `.itp`**: [file:313]
    - Keep everything in `conf_UNK.itp` before the UNK moleculetype block.
    - Insert the normalized UNK block.
    - Append the bodies (sections) from all extras, now with shifted indices.
    - Append everything after the UNK block from `conf_UNK.itp`.
    - Write to `output.itp`. [file:313]

---

## 3. `merge_3c.py`

### Purpose

Given:

- An averaged snapshot `NA_avg_dist.gro` and
- Its topology `topol.top` and
- A definition of the ASEC complex in `Infos.dat`,

this script:

1. Merges all solute-selected residues (e.g. UNK plus specific residues by `resid`) into a single UNK residue in the GRO file. [file:312]
2. Renumbers atom indices sequentially.
3. Adjusts the `[ molecules ]` section in `topol.top` to subtract molecules that have been merged into UNK. [file:312]

Outputs:

- `NA_avg_dist_merged.gro` – updated coordinates.
- `topol_merged.top` – updated topology. [file:312]

These are then renamed by `asec_gensnapshots.sh`. [file:311][file:312]

### Inputs

- `Infos.dat` – must contain an `ASEC` line, for example: [file:312]

```text
ASEC resname UNK or resid 1163 or resname CL
```

- `NA_avg_dist.gro` – initial GRO (a copy of `<Project>_avg_dist.gro`). [file:312]
- `topol.top` – topology corresponding to `NA_avg_dist.gro`. [file:312]


### ASEC selection

The script parses the `ASEC` line in `Infos.dat` and extracts: [file:312]

- All `resname X` tokens → set `resnames_from_infos`.
- All `resid N` tokens → set `resids_from_infos` (as integers).

Example: `ASEC resname UNK or resid 1163` selects:

- Any residue with resname `UNK`.
- The residue with resid `1163`. [file:312]


### GRO merging steps

1. **Parse GRO file**: [file:312]
    - Read title, `natoms`, atom lines, and box line.
    - For each atom line, parse:

```python
resnr   = int(line[0:5])
resname = line[5:10].strip()
atom    = line[10:15].strip()
atomnr  = int(line[15:20])
coords  = line[20:].rstrip("\n")
```

    - Group atoms by `(resnr, resname)` into `residue_atoms`. [file:312]
2. **Identify base UNK residue**: [file:312]
    - Find all residue IDs with `resname == "UNK"`, take the smallest `resnr` as `unk_resid`. [file:312]
    - This residue is the base into which others will be merged (`base_key = (unk_resid, "UNK")`). [file:312]
3. **Build ASEC selection**: [file:312]
    - For each `(resnr, resname)`:
        - Select if `resname` is in `resnames_from_infos` or `resnr` is in `resids_from_infos`.
    - Ensure `base_key` is always included. [file:312]
4. **Construct merged residue and the rest**: [file:312]
    - **Merged UNK residue**:

```python
merged_unk_atoms = []
for key in sorted(selected_keys):
    for (resnr, resname, atom, atomnr, coords) in residue_atoms[key]:
        merged_unk_atoms.append((unk_resid, "UNK", atom, atomnr, coords))
```

    - **Other residues**: all atoms in residues not in `selected_keys` are left unchanged. [file:312]
    - Concatenate:

```python
new_atoms = merged_unk_atoms + other_atoms
```

5. **Renumber atom indices** to 1..N, preserving residue ids and names: [file:312]

```python
renum_atoms = []
for i, (resnr, resname, atom, atomnr, coords) in enumerate(new_atoms, start=1):
    renum_atoms.append((resnr, resname, atom, i, coords))
```

6. **Write `NA_avg_dist_merged.gro`** with the same title, new atom list, and original box line. [file:312]

### Topology update

1. Count how many residues of each resname (except UNK) are in the merged set: [file:312]

```python
merged_resname_counts = defaultdict(int)
for (resnr, resname) in selected_keys:
    if resname != "UNK":
        merged_resname_counts[resname] += 1
```

2. Parse `topol.top` line by line, focusing on the `[ molecules ]` section: [file:312]
    - For each molecule entry `name   count`:
        - If `name == "UNK"` → leave count unchanged.
        - If `name` is in `merged_resname_counts`:
            - Set `new_count = count - merged_resname_counts[name]`.
            - If `new_count < 0`, abort with an error. [file:312]
        - Otherwise, keep line unchanged.
3. Write the updated topology to `topol_merged.top`. [file:312]

This ensures that the molecules whose residues were merged into UNK are removed from the topology counts, while the overall atom content matches the new GRO. [file:312]

---

## 4. `asec_mdrun.sh`

### Purpose

Run the NPT ASEC MD on the prepared snapshot/system in the `asec/` directory. [file:310]

It:

1. Reads `Infos.dat` (for logging consistency).
2. Builds `npt_asec.tpr` from `npt_asec.mdp`, `${Project}_avg_dist.gro`, `topol.top`, and `index.ndx`. [file:310]
3. Runs `mdrun_gpu_mpi` to perform the ASEC MD. [file:310]

### Usage

From inside `asec/` (after `asec_gensnapshots.sh` has run):

```bash
bash asec_mdrun.sh
```

or via your scheduler (as set in `asec_gensnapshots.sh`, e.g. through `mBASH`). [file:311][file:310]

### Steps

1. **Environment setup**: loads GROMACS module, MPI env vars, and Conda environment. [file:310]
2. **Read `Infos.dat`** for `Solvent`, `Project`, `Resname`, `Boxsize`, `Ions`, `Fixed`, `Posre`, `Temperature` (printed for logging). [file:310]
3. **Prepare NPT ASEC tpr**: [file:310]

```bash
gmx grompp \
    -f npt_asec.mdp \
    -c ${Project}_avg_dist.gro \
    -p topol.top \
    -o npt_asec.tpr \
    -n index.ndx \
    -r ${Project}_avg_dist.gro \
    -maxwarn 1
```

4. **Run MD**: [file:310]

```bash
mpiexec --bind-to none mdrun_gpu_mpi -deffnm npt_asec
```


Outputs:

- `npt_asec.tpr`, `npt_asec.xtc`, `npt_asec.edr`, `npt_asec.gro` and associated log/trajectory files. [file:310]

You can then use this ASEC trajectory for downstream QM/MM or field calculations centered on the merged UNK complex.

```
<span style="display:none">[^1][^2][^3][^4]</span>

<div align="center">⁂</div>

[^1]: asec_mdrun.sh
[^2]: asec_gensnapshots.sh
[^3]: merge_3c.py
[^4]: merge_itp.py```

