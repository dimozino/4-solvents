#!/usr/bin/env python3
import glob
import os
import re

infos_file = "Infos.dat"
pdb_glob   = "snapshots/frame_*.pdb"
out_dir    = "snapshots_renamed"
target_resname = "UNK"

os.makedirs(out_dir, exist_ok=True)

# 1. Read NUMA from Infos.dat
with open(infos_file) as f:
    text = f.read()

m = re.search(r"^NUMA\s+(\d+)", text, flags=re.MULTILINE)
if not m:
    raise SystemExit("Could not find 'NUMA N' line in Infos.dat")
NUMA = int(m.group(1))
print(f"NUMA = {NUMA}")

# 2. Process each snapshot
pdb_files = sorted(glob.glob(pdb_glob))
if not pdb_files:
    raise SystemExit(f"No PDB files matched {pdb_glob}")

for pdb in pdb_files:
    with open(pdb) as f:
        lines = f.readlines()

    atom_count = 0
    new_lines = []

    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            atom_count += 1
            if atom_count <= NUMA:
                # Change residue name (columns 18-20) to UNK
                # PDB is 1-based: resname at columns 18-20 -> indices 17:20
                new_line = line[:17] + f"{target_resname:>3s}" + line[20:]
                new_lines.append(new_line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    out_name = os.path.join(out_dir, os.path.basename(pdb))
    with open(out_name, "w") as f:
        f.writelines(new_lines)

    print(f"Wrote {out_name}, renamed first {min(atom_count, NUMA)} atoms to resname {target_resname}")

