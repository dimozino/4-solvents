## Overview

`xyz_psf_to_pdb.py` combines:

- Cartesian coordinates from an `.xyz` file, and
- atom naming/residue metadata from a CHARMM PSF file

to produce a PDB file suitable for visualization or further MD input preparation.[^2][^3]

It assumes a one-to-one, ordered correspondence between atoms in the XYZ and PSF files.

***

## Requirements

- Python ≥ 3.7
- No external dependencies (uses only `sys`).[^1]

***

## Input Files

### XYZ file

The script expects a standard XYZ format:[^2]

1. Line 1: number of atoms, `N`.
2. Line 2: comment (ignored).
3. Lines 3…(N+2): one atom per line:

```text
Elem   x   y   z
```

where:
    - `Elem` is the element symbol (e.g. `C`, `O`, `H`).
    - `x,y,z` are Cartesian coordinates (typically in Å).[^2]

The script reads exactly `N` atomic records starting at line 3 and raises an error if the number of parsed coordinates does not match the header.[^2]

### PSF file

The script reads atom metadata from the `!NATOM` section of a CHARMM/XPLOR-style PSF:[^3]

- It searches for the line containing `!NATOM`.
- The first integer on that line is interpreted as the number of atoms, `N`.
- The next `N` lines must have standard PSF atom records:

```text
idx  segid  resid  resname  atomname  atomtype  charge  mass  ...
```


For each atom, the script uses:

- `idx` – atom index (becomes PDB serial number).
- `segid` – segment ID (first character becomes chain ID).
- `resid` – residue index.
- `resname` – residue name.
- `atomname` – atom name (goes into PDB atom name field).[^3]

The atom count in PSF must equal the atom count in the XYZ.

***

## Assumptions and Matching Logic

- Atoms are **aligned by index**: the first atom in the XYZ corresponds to PSF atom index 1, the second to index 2, etc.[^3][^2]
- There is no attempt to match by element or name; the order must be consistent.
- Typically, this is used for a single small-molecule residue (e.g. `UNL`), but multiple residues/segments are supported as long as ordering matches.[^3]

***

## Usage

Command line:

```bash
python xyz_psf_to_pdb.py 3-S0.xyz mol.psf output.pdb
```

Arguments:[^1]

1. `mol.xyz` – input XYZ file with coordinates.
2. `mol.psf` – CHARMM PSF providing atom names, residue info, segment IDs.
3. `output.pdb` – path to the PDB file to write.

If the number of arguments is not exactly three, the script prints:

```text
Usage: python xyz_psf_to_pdb.py mol.xyz mol.psf output.pdb
```

and exits.[^1]

***

## Output PDB Format

For each atom, the script writes an `ATOM` record with standard PDB formatting:[^1][^3]

- Record type: `ATOM  `
- Serial number: PSF `idx`.
- Atom name: PSF `atomname`.
- Residue name: PSF `resname`.
- Chain ID: first character of PSF `segid` (falls back naturally if `segid` is a single character, e.g. `A`).[^3]
- Residue sequence number: PSF `resid`.
- Coordinates: `x, y, z` from the XYZ file, printed with 3 decimals.[^2]
- Occupancy: fixed to `1.00`.
- B-factor: fixed to `0.00`.
- Element: taken from the first column (`Elem`) of the XYZ line, right-justified in columns 77–78.[^2]

An `END` record is written at the end of the file.[^1]

Example `ATOM` line structure:

```text
ATOM      1  C0  UNL A   1      0.641   1.771  -4.579  1.00  0.00           C
```

- `C0` from PSF atomname, `UNL`/`1` from PSF residue, `C` from XYZ element, coordinates from XYZ.[^3][^2]

***

## Error Handling

The script raises `ValueError` in two main cases:[^1][^2][^3]

- XYZ header atom count does not match the number of coordinate lines actually parsed.
- PSF and XYZ have different atom counts when writing the PDB.

These checks help catch misaligned or corrupted input before generating an inconsistent PDB.

***

## Typical Workflow

1. Optimize or obtain a ligand geometry in quantum chemistry / builder, export as XYZ.
2. Build a CHARMM PSF for the same ligand (same atom order!).[^2][^3]
3. Run:

```bash
python xyz_psf_to_pdb.py mol.xyz mol.psf ligand.pdb
```

4. Inspect `ligand.pdb` in VMD/PyMOL, then use it as input for GROMACS topology building in combination with your `.itp` from `charmm2gmx_itp.py`.[^3][^1][^2]

If you want, I can add optional flags to override chain ID/residue name or to support multiple XYZ frames.

<div align="center">⁂</div>

[^1]: charmm2gmx_itp.py

[^2]: 3-S0.xyz.txt

[^3]: mol.psf.txt

