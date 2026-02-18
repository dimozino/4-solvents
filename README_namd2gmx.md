## Overview

`charmm2gmx_itp.py` converts a CHARMM-style topology and parameter set for a single small molecule into a GROMACS `.itp` file.[^1]
It reads:

- `mol.psf`: CHARMM PSF containing atoms, bonds, angles, dihedrals, impropers.[^2]
- `ff.par`: Custom CHARMM parameter file with BONDS, THETAS, PHI and NONBONDED terms for this molecule.[^3]

It writes a GROMACS topology (`.itp`) with:

- `[ atomtypes ]` from `NONBONDED`
- `[ moleculetype ]`
- `[ atoms ]`
- `[ bonds ]`
- `[ angles ]`
- `[ dihedrals ]` (proper, function type 9)
- Optional impropers (placeholders)

***

## Requirements

- Python ≥ 3.7
- No external libraries beyond the standard library (`sys`, `math`, `collections`).[^1]

***

## Input File Requirements

### PSF (`mol.psf`)

The script assumes a standard CHARMM/XPLOR PSF layout:[^2]

- `!NATOM` section with lines
`idx segid resid resname atomname atomtype charge mass ...`
- `!NBOND`, `!NTHETA`, `!NPHI`, `!NIMPHI` sections with integer index lists for bonds, angles, dihedrals, impropers.[^2]

Atom **names** (e.g. `C0`, `O5`) must correspond to the names used in `ff.par` (BONDS/THETAS/PHI/NONBONDED).[^3][^2]

### Parameter file (`ff.par`)

The script expects the following blocks and formats:[^3][^1]

- `BONDS`
`A  B  Kb  r0`
where `Kb` in kcal/mol/Å², `r0` in Å.
- `THETAS`
`A  B  C  Ktheta  theta0`
where `Ktheta` in kcal/mol/rad², `theta0` in degrees.
- `PHI` (proper torsions)
`A  B  C  D  k  n  delta`
possibly repeated for the same quadruple, one line per multiplicity, with `k` in kcal/mol, `delta` in degrees.[^3]
- `NONBONDED`
One header line starting with `NONBONDED` and parameters line:
`type  0.000  eps  Rmin  0.000  eps14  Rmin14`
where `eps` in kcal/mol, `Rmin` in Å.[^3]
- `END` to terminate the file.

The atom names `A,B,C,D,type` must match PSF atom names.[^2][^3]

***

## Units and Conversions

The script converts CHARMM units to GROMACS as follows:[^1]

- Bond force constant
$K_b$: kcal/mol/Å² → kJ/mol/nm² via factor 836.8.
- Bond length
$r_0$: Å → nm via × 0.1.
- Angle force constant
$K_\theta$: kcal/mol/rad² → kJ/mol/rad² via factor 8.368.
- Angle equilibrium angle
$\theta_0$: kept in degrees.
- Dihedral amplitude
$k$: kcal/mol → kJ/mol via × 4.184.
- Nonbonded epsilon
`eps`: kcal/mol → kJ/mol via × (−4.184), preserving attractive sign.[^1]
- Nonbonded Rmin
`Rmin`: Å → nm via × 0.1, then converted to $\sigma$ assuming
$R_{\min} = 2^{1/6}\sigma$, so $\sigma = R_{\min} \cdot 2 / 2^{1/6}$.[^1]

***

## Usage

From the command line:

```bash
python charmm2gmx_itp.py mol.psf ff.par conf_UNK.itp
```

Arguments:[^1]

1. `mol.psf` – CHARMM PSF file for a single small-molecule residue.
2. `ff.par` – CHARMM-style parameter file described above.
3. `output.itp` – Path to output GROMACS topology.

The script will exit with a short usage message if you do not provide exactly three arguments.[^1]

***

## Output Structure

The generated `.itp` contains:[^1]

### [ atomtypes ]

- One line per unique name from `NONBONDED` (e.g. `C0`, `O5`).[^3][^1]
- Columns: `name at.num mass charge ptype sigma epsilon`.
- Element and mass guessed from the first character of the name (`C`, `H`, `O`, etc.).[^1]
- `sigma` and `epsilon` are derived from `ff.par` as above.


### [ moleculetype ]

- Name: residue name from the first PSF atom (e.g. `UNL`).[^2][^1]
- `nrexcl` is fixed to 3 (standard for all-atom FFs).[^1]


### [ atoms ]

- Each PSF atom becomes one `[ atoms ]` entry.
- GROMACS atom **type** and atom **name** both use PSF `atomname` (e.g. `C0`).[^1]
- Residue index and name from PSF `resid`, `resname`; charge and mass directly from PSF.[^2][^1]


### [ bonds ]

- All PSF `!NBOND` entries are emitted.[^2][^1]
- Parameters look up in `BONDS` by sorted atom **names** `(A,B)`. If missing, a fallback `(kb=0, r0=0.1 nm)` is written and should be corrected manually.[^3][^1]


### [ angles ]

- All PSF `!NTHETA` entries are emitted.[^2][^1]
- Parameters look up in `THETAS` by `(A,B,C)` and also reversed `(C,B,A)` if not found.[^3][^1]
- If still missing, a generic angle is written with `theta0 = 109.5°`, `ktheta = 0.0` (should be fixed manually).[^1]


### [ dihedrals ]

- All PSF `!NPHI` entries are emitted once; multiple PHI lines per quadruple become multiple GROMACS torsions.[^2][^3][^1]
- Dihedral parameters are matched by atom **names** `(A,B,C,D)` or the reverse `(D,C,B,A)`.
- Each PHI term is written as function type 9:
`ai aj ak al  9  phi0  k  mult`
where `phi0=delta` in degrees, `k` in kJ/mol, `mult=n`.[^1]
- If no PHI entry exists, a single zeroed torsion is written for that quadruple.


### Impropers

- PSF `!NIMPHI` entries are written as function type 4 impropers with placeholder parameters (`180.0 0.0`).[^2][^1]
- You should supply proper improper parameters manually if needed.

***

## Algorithm and Matching Logic

- **PSF parsing** uses section tags `!NATOM`, `!NBOND`, `!NTHETA`, `!NPHI`, `!NIMPHI` and reads the specified number of entries; connectivity is stored as integer indices.[^2][^1]
- **Parameter matching** is done by PSF `atomname` strings:
    - Bonds: `BONDS` keyed by sorted `(name_i, name_j)`.
    - Angles: `THETAS` keyed by `(name_i, name_j, name_k)` and reversed.
    - Dihedrals: `PHI` keyed by `(name_i, name_j, name_k, name_l)` or reversed.[^3][^1]

This means that if you change atom names in the PSF, you must update `ff.par` accordingly.

***

## Typical Workflow

1. Generate `mol.psf` for the ligand in CHARMM/XPLOR format.
2. Prepare `ff.par` with all bonded and nonbonded parameters for the same atom names.
3. Run:

```bash
python charmm2gmx_itp.py mol.psf ff.par conf_UNK.itp
```

4. Inspect `conf_UNK.itp`:
    - Check for any bonds/angles/dihedrals with zero force constants.
    - Verify nonbonded `sigma`/`epsilon` against your CHARMM parameters.
5. Include the `.itp` into your GROMACS system topology with `#include "conf_UNK.itp"`.

***

## Limitations and Notes

- Designed for a **single small-molecule residue**, not full proteins or multi-residue ligands.[^1]
- Nonbonded conversion assumes the `NONBONDED` block uses the specific layout in your `ff.par`; if that changes, the parser must be adapted.[^3][^1]
- Improper dihedral parameters are not read from `IMPHI`; only connectivity is transferred with dummy parameters.[^1]
- Dihedral mapping relies on GROMACS function type 9 (CHARMM-style periodic torsions); ensure your force field settings are compatible.

If you want, I can also produce a short in-script `--help` text or a LaTeX-style methods paragraph describing this conversion for a paper.

<div align="center">⁂</div>

[^1]: charmm2gmx_itp.py

[^2]: mol.psf.txt

[^3]: ff.par.txt

