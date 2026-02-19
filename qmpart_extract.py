#!/usr/bin/env python3
import math

pdb_in  = "snapshots/frame_0.pdb"
xyz_out = "qmpart.xyz"
target_resname = "UNK"

atoms = []

with open(pdb_in) as f:
    for line in f:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        resname = line[17:20].strip()
        if resname != target_resname:
            continue
        name = line[12:16].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        # Map PDB atom name to a simple element symbol for XYZ
        # Take first letter, uppercase; special-case CL, SI, etc. if needed
        up = name.upper()
        if up.startswith("CL"):
            elem = "Cl"
        elif up.startswith("SI"):
            elem = "Si"
        else:
            elem = up[0]

        atoms.append((elem, x, y, z))

if not atoms:
    raise SystemExit(f"No atoms with resname {target_resname} found in {pdb_in}")

with open(xyz_out, "w") as f:
    f.write(f"{len(atoms)}\n")
    f.write(f"{target_resname} extracted from {pdb_in}\n")
    for elem, x, y, z in atoms:
        f.write(f"{elem}  {x: .6f}  {y: .6f}  {z: .6f}\n")

print(f"Wrote {xyz_out} with {len(atoms)} atoms")

