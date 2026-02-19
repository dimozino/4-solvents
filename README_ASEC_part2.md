# ASEC Post‑processing and QM Input Preparation

This toolkit processes ASEC MD trajectories to:

- Extract evenly spaced snapshots.
- Build equal-sized ASEC environments around the UNK complex.
- Assign point charges and merge frames into a single `.pc` file.
- Extract the QM subsystem geometry for subsequent QM/MM or pure QM runs.

Scripts:

- `asec_extract.sh` – main driver: frame extraction, ASEC neighbor selection, charge mapping, QM part extraction. [file:316]
- `asec_extract.py` – selects full solvent/residue shells around ASEC atoms and equalizes atom counts across frames. [file:315]
- `asec_merge.py` – converts neighbor PDBs plus residue charge tables into a merged point-charge file. [file:318]
- `qmpart_extract.py` – extracts the UNK residue into a QM geometry file (`qmpart.xyz`). [file:317]
- `equalize_asecs.sh` – batch processing over multiple systems listed in `folders.txt`, including equalization and merging. [file:319]

---

## 1. `asec_extract.sh`

### Purpose

Runs the full ASEC extraction pipeline for a **single** system with an existing NPT ASEC trajectory:

1. Analyzes `npt_asec.xtc` to determine a time step for ~100 evenly spaced frames. [file:316]
2. Extracts snapshots as separate PDBs (`snapshots/frame_*.pdb`). [file:316]
3. Renames the first NUMA residue to UNK (via an external helper script). [file:316]
4. Uses `asec_extract.py` to build equal‑size neighbor shells around the ASEC group. [file:316][file:315]
5. Uses `asec_merge.py` + charge tables to generate a merged point‑charge file. [file:316][file:318]
6. Uses `qmpart_extract.py` to extract the QM subsystem geometry. [file:316][file:317]

### Expected inputs

In the working directory: [file:316]

- `npt_asec.xtc`, `npt_asec.tpr` – ASEC trajectory and corresponding TPR.
- `index.ndx` – must contain ASEC group named `UNK` (used downstream by `asec_extract.py`). [file:315]
- `Infos.dat` – must contain a `Cutoff` line, e.g. `Cutoff 5.0` (Å), used by `asec_extract.py`. [file:315]
- Access to scripts and data in `sd=/home/ibsbnicigt/nikolaev_d/script_library/4-solvents/`:
  - `rename_first_NUMA_to_UNK.py`
  - `asec_extract.py`
  - `THF_charges.txt`, `NA_charges.txt`
  - `asec_merge.py`
  - `qmpart_extract.py` [file:316]

### Workflow

1. **Environment**: loads GROMACS module, MPI variables, and activates Conda. [file:316]

2. **Check trajectory and estimate dt**: [file:316]

   ```bash
   gmx check -f npt_asec.xtc 2> check.log
```

Then a temporary Python script (`get_dt.py`) parses `check.log` for a line like
`Last frame  N  time  T` and sets:

$$
\text{dt} \approx T / 100
$$

rounded to an integer (ps). [file:316]

3. **Extract ~100 frames as PDBs**: [file:316]

```bash
mkdir -p snapshots
echo System | gmx trjconv -s npt_asec.tpr -f npt_asec.xtc \
    -o snapshots/frame_.pdb -dt "$dt" -sep
```

4. **Rename NUMA to UNK**: copies and runs `rename_first_NUMA_to_UNK.py`, which creates `snapshots_renamed/`; the script replaces the old directory: [file:316]

```bash
cp $sd/rename_first_NUMA_to_UNK.py .
python3 rename_first_NUMA_to_UNK.py
rm -r snapshots
mv snapshots_renamed/ snapshots
```

5. **ASEC neighbor extraction**: [file:316][file:315]

```bash
cp $sd/asec_extract.py .
python3 asec_extract.py
rm neighbors_eq_fullres/frame_0.pdb
```

This produces `neighbors_eq_fullres/frame_*.pdb`, each containing the same number of atoms (neighbors within `Cutoff` Å).
6. **Charge mapping and point‑charge file**: [file:316][file:318]

```bash
cp $sd/THF_charges.txt .
cp $sd/NA_charges.txt .
cp $sd/asec_merge.py .
python3 asec_merge.py
mv merged_frames.dat field_asec.pc
```

7. **QM subsystem extraction**: [file:316][file:317]

```bash
cp $sd/qmpart_extract.py .
python3 qmpart_extract.py
```

Produces `qmpart.xyz` from `snapshots/frame_0.pdb`.

---

## 2. `asec_extract.py`

### Purpose

From a set of ASEC trajectory snapshots (`snapshots/frame_*.pdb`):

- Identify all residues within `Cutoff` Å of the ASEC group (UNK) in each frame. [file:315]
- Enforce that **all frames have the same number of neighbor atoms**, by selecting the nearest full residues up to a global target count. [file:315]
- Write trimmed neighbor PDBs to `neighbors_eq_fullres/frame_*.pdb`. [file:315]


### Inputs

In the working directory: [file:315]

- `index.ndx` – must contain a group named `UNK` (ASEC atoms). [file:315]
- `Infos.dat` – must include a `Cutoff` line, e.g. `Cutoff 5.0`. [file:315]
- Snapshot files: `snapshots/frame_*.pdb` (from `asec_extract.sh`). [file:315]


### Algorithm overview

1. **Read cutoff** from `Infos.dat`: [file:315]
    - Scan for line starting with `Cutoff`.
    - Use the second token as a float; abort if not found.
2. **Get ASEC atom indices** from `index.ndx` group `UNK`. [file:315]
3. **Parse each PDB** (`ATOM`/`HETATM` records) into per‑atom dicts (serial, name, resname, chain, resid, xyz). [file:315]
4. **Map ASEC coordinates** by matching `serial` with ASEC indices. [file:315]
5. **Group by residue** (`(chain, resid, resname)`). [file:315]
6. For each residue that is not purely ASEC:
    - Compute minimal squared distance from any atom in the residue (excluding ASEC atoms) to any ASEC atom. [file:315]
    - If $\text{min} d^2 \le \text{Cutoff}^2$, mark residue as neighbor and store `mindist2`. [file:315]
7. For each frame:
    - Sort neighbor residues by `mindist2`.
    - Compute cumulative atom counts by including full residues in order; store this **per‑frame budget**. [file:315]
8. **Global target** atom count = minimum over per‑frame budgets; abort if zero. [file:315]
9. For each frame, **select residues** in order of proximity until adding another residue would exceed `target_atoms` (with a small guard condition). [file:315]
10. Build `out_atoms` list (all atoms from chosen residues), renumber serials from 1, and write PDB: [file:315]
    - First line: `REMARK neighbors of ASEC within X A, full residues, N atoms`.
    - Then all atom records (updated serials).
    - Final line: `END`. [file:315]

Result: All `neighbors_eq_fullres/frame_*.pdb` have **identical atom counts** and consist of complete residues around the ASEC. [file:315]

---

## 3. `asec_merge.py`

### Purpose

Convert all equalized ASEC neighbor frames plus per‑residue charge tables into a merged point‑charge file for QM calculations:

- Reads `neighbors_eq_fullres/frame_*.pdb`. [file:318]
- Loads charges from `*_charges.txt` (e.g. `THF_charges.txt`, `NA_charges.txt`). [file:318]
- For every atom in every frame, associates a partial charge and coordinate. [file:318]
- Writes them in a combined file `merged_frames.dat` (later renamed to `field_asec.pc`). [file:318]


### Inputs

- PDBs: `neighbors_eq_fullres/frame_*.pdb`. [file:318]
- Charge tables: one or more `RESNAME_charges.txt` files in the current directory, with format: [file:318]

```text
; comment lines starting with #
ATOMNAME  charge_in_e
...
```

Example (`THF_charges.txt`): [file:320]

```text
C1   0.123
C2   0.000
...
```


### Algorithm

1. **Load charge tables**: [file:318]
    - For each `*_charges.txt`, derive `resname` from filename (prefix before `_charges.txt`).
    - Parse non‑empty, non‑comment lines; map `atom_name → charge`. [file:318]
    - Build `charge_tables[resname] = {atom_name: q}`.
2. **Read all frames**: [file:318]
    - `pdb_files = sorted(glob.glob("neighbors_eq_fullres/frame_*.pdb"))`.
    - For each PDB, collect all `ATOM`/`HETATM` entries as tuples `(resname, atom_name, x, y, z)`. [file:318]
    - Enforce that **all frames have the same number of atoms** (`expected_natoms`). [file:318]
3. **Map to charges** per atom: [file:318]
    - For each atom:
        - Look up `resname` in `charge_tables` (error if missing).
        - Look up `atom_name` in that residue’s table (error if missing). [file:318]
        - Append `(q, x, y, z)` to `frame_vals`. [file:318]
    - Append `frame_vals` to global `all_coords`. [file:318]
4. **Write merged file**: [file:318]

```text
N
q1 x1 y1 z1
q2 x2 y2 z2
...
```

    - First line: total number of entries across all frames (`len(all_coords)`). [file:318]
    - Charges are scaled by 1/100: `q/100.0` in the output. [file:318]

The calling scripts typically rename `merged_frames.dat` to `field_asec.pc`. [file:316][file:319]

---

## 4. `qmpart_extract.py`

### Purpose

Extract the QM subsystem geometry (UNK residue) from a reference snapshot into XYZ format.

- Input: `snapshots/frame_0.pdb`. [file:317]
- Output: `qmpart.xyz`. [file:317]


### Behavior

1. Reads `snapshots/frame_0.pdb`. [file:317]
2. For each `ATOM`/`HETATM` line:
    - Keep only atoms with `resname == "UNK"`. [file:317]
    - Read PDB columns for atom name and coordinates. [file:317]
3. Maps PDB atom names to element symbols for XYZ: [file:317]
    - `CL*` → `Cl`
    - `SI*` → `Si`
    - Else: first letter of uppercase atom name (e.g. `C`, `N`, `O`). [file:317]
4. Writes `qmpart.xyz`: [file:317]

```text
N
UNK extracted from snapshots/frame_0.pdb
Elem  x  y  z
...
```

where `N` is the number of UNK atoms. [file:317]

Errors out if no UNK atoms are found. [file:317]

---

## 5. `equalize_asecs.sh`

### Purpose

Batch ASEC equalization and field/QM extraction over multiple systems listed in `folders.txt`.

For each project `x`:

1. Enforce that all `x/asec/neighbors_eq_fullres/frame_*.pdb` share a common atom count across projects. [file:319]
2. In each `x/asec/`, run `asec_merge.py` and `qmpart_extract.py`. [file:319]
3. Write per‑system field and QM files: `${x}_field_asec.pc`, `${x}_qmpart.xyz`. [file:319]

### Inputs

- `folders.txt` – list of project directories (`x`, e.g. `0a`, `1a`, …), one per line. [file:319]
- In each `x/asec/neighbors_eq_fullres/`:
    - `frame_1.pdb ... frame_100.pdb`. [file:319]
- Script directory `sd` with `asec_merge.py`, `qmpart_extract.py`, and `*_charges.txt`. [file:319]


### Step 1: find global smallest atom count

For each `x`: [file:319]

- Count atoms in `x/asec/neighbors_eq_fullres/frame_1.pdb` using `grep -Ec '^(ATOM |HETATM)'`. [file:319]
- Keep the minimum positive count as `smallest_number`. [file:319]


### Step 2: trim all frames to `smallest_number`

For each `x` and for `frame_1.pdb` to `frame_100.pdb`: [file:319]

- Use an `awk` filter:

```awk
BEGIN { n=0 }
/^ATOM / || /^HETATM/ {
    n++
    if (n <= maxatoms) print
    next
}
/^END/ { next }  # skip existing END
{ print }
END { print "END" }
```

- This keeps only the first `maxatoms = smallest_number` ATOM/HETATM lines, preserves other lines (e.g. REMARK), and ensures a single `END` per file. [file:319]

After this pass, **all frames in all projects** have the same atom count. [file:319]

### Step 3: per‑project merging and QM extraction

For each `x`: [file:319]

- `cd x/asec/`
- Copy `asec_merge.py` and run it; rename output:

```bash
cp $sd/asec_merge.py .
python3 asec_merge.py
mv merged_frames.dat field_asec.pc
```

- Copy and run `qmpart_extract.py`, then rename outputs: [file:319]

```bash
cp $sd/qmpart_extract.py .
python3 qmpart_extract.py
mv field_asec.pc ${x}_field_asec.pc
mv qmpart.xyz ${x}_qmpart.xyz
```

- Return to the original directory (`currdir`). [file:319]

Outputs are left in each `x/asec/` directory, named per project. [file:319]

```
<span style="display:none">[^1][^2][^3][^4][^5][^6][^7]</span>

<div align="center">⁂</div>

[^1]: asec_extract.py
[^2]: asec_extract.sh
[^3]: qmpart_extract.py
[^4]: asec_merge.py
[^5]: equalize_asecs.sh
[^6]: THF_charges.txt
[^7]: rename_first_NUMA_to_UNK.py```

