#!/usr/bin/env python3
"""
xyz_psf_to_pdb.py

Use an XYZ structure and a CHARMM PSF file to generate a PDB file.

Usage:
    python xyz_psf_to_pdb.py 3-S0.xyz mol.psf output.pdb
"""

import sys

def read_xyz(xyz_path):
    with open(xyz_path, "r") as f:
        lines = f.readlines()

    # First line: natoms, second line: comment
    natoms = int(lines[0].strip())
    coords = []
    for line in lines[2:2 + natoms]:
        if not line.strip():
            continue
        parts = line.split()
        elem = parts[0]
        x, y, z = map(float, parts[1:4])
        coords.append((elem, x, y, z))
    if len(coords) != natoms:
        raise ValueError(f"XYZ file {xyz_path} has {natoms} atoms declared but {len(coords)} coords parsed.")
    return coords

def read_psf_atoms(psf_path):
    with open(psf_path, "r") as f:
        lines = f.readlines()

    atoms = []
    i = 0
    n = len(lines)

    # Skip to !NATOM
    while i < n and "!NATOM" not in lines[i]:
        i += 1
    if i == n:
        raise ValueError("!NATOM not found in PSF.")
    natom = int(lines[i].split()[0])
    i += 1

    for _ in range(natom):
        fields = lines[i].split()
        idx = int(fields[0])
        segid = fields[1]
        resid = int(fields[2])
        resname = fields[3]
        atomname = fields[4]
        atomtype = fields[5]
        charge = float(fields[6])
        mass = float(fields[7])
        atoms.append(
            dict(
                index=idx,
                segid=segid,
                resid=resid,
                resname=resname,
                atomname=atomname,
                atomtype=atomtype,
                charge=charge,
                mass=mass,
            )
        )
        i += 1

    return atoms

def write_pdb(pdb_path, atoms, coords):
    if len(atoms) != len(coords):
        raise ValueError(f"Atom count mismatch: PSF has {len(atoms)}, XYZ has {len(coords)}")

    with open(pdb_path, "w") as out:
        for i, (atom, (elem, x, y, z)) in enumerate(zip(atoms, coords), start=1):
            # Standard PDB ATOM record formatting
            # Columns:
            # 1-6  Record name "ATOM  "
            # 7-11 Serial
            # 13-16 Atom name
            # 17   altLoc
            # 18-20 Residue name
            # 22   Chain ID
            # 23-26 Residue sequence
            # 31-38 x, 39-46 y, 47-54 z
            # 77-78 element
            serial = atom["index"]
            atomname = atom["atomname"]
            resname = atom["resname"]
            resid = atom["resid"]
            chain = atom["segid"][0] if atom["segid"] else "A"

            out.write(
                f"ATOM  {serial:5d} {atomname:>4s} {resname:>3s} {chain:1s}"
                f"{resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00          {elem:>2s}\n"
            )
        out.write("END\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python xyz_psf_to_pdb.py 3-S0.xyz mol.psf output.pdb")
        sys.exit(1)

    xyz_path = sys.argv[1]
    psf_path = sys.argv[2]
    pdb_path = sys.argv[3]

    coords = read_xyz(xyz_path)
    atoms = read_psf_atoms(psf_path)
    write_pdb(pdb_path, atoms, coords)

if __name__ == "__main__":
    main()

