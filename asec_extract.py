#!/usr/bin/env python3

import os
import glob
import math

index_file = "index.ndx"
infos_file = "Infos.dat"
pdb_glob = "snapshots/frame_*.pdb"
out_dir = "neighbors_eq_fullres"
asec_group_name = "UNK"

# --- 0. Read cutoff from Infos.dat ---
cutoff = None
if os.path.isfile(infos_file):
    with open(infos_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0] == "Cutoff" and len(parts) >= 2:
                cutoff = float(parts[1])
                break

if cutoff is None:
    raise SystemExit(f"Cutoff not set in {infos_file} (expected line like 'Cutoff 5.0')")

print(f"Using cutoff = {cutoff} Å from {infos_file}")

os.makedirs(out_dir, exist_ok=True)

# --- 1. ASEC atoms from index.ndx ---
asec_atoms = set()
with open(index_file) as f:
    current = None
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith('[') and line.endswith(']'):
            current = line[1:-1].strip()
            continue
        if current == asec_group_name:
            for tok in line.split():
                asec_atoms.add(int(tok))

if not asec_atoms:
    raise SystemExit(f"No atoms in group {asec_group_name} in {index_file}")
print(f"ASEC atoms: {len(asec_atoms)}")

# --- 2. Helpers ---
def read_pdb_atoms(fname):
    atoms = []
    with open(fname) as f:
        for line in f:
            if line.startswith(("ATOM ", "HETATM")):
                serial  = int(line[6:11])
                name    = line[12:16].strip()
                resname = line[17:20].strip()
                chain   = line[21].strip()
                resid   = int(line[22:26])
                x       = float(line[30:38])
                y       = float(line[38:46])
                z       = float(line[46:54])
                atoms.append(dict(
                    line=line, serial=serial, name=name,
                    resname=resname, chain=chain, resid=resid,
                    x=x, y=y, z=z
                ))
    return atoms

pdb_files = sorted(glob.glob(pdb_glob))
if not pdb_files:
    raise SystemExit(f"No PDB files matched {pdb_glob}")
print(f"Found {len(pdb_files)} PDB snapshots")

# --- 3. Per-frame neighbor residues, full residues only ---
frames_data = []  # list of dicts per frame

for pdb in pdb_files:
    atoms = read_pdb_atoms(pdb)
    serial_map = {a["serial"]: a for a in atoms}

    # ASEC coordinates
    asec_coords = []
    for aid in asec_atoms:
        if aid in serial_map:
            aa = serial_map[aid]
            asec_coords.append((aa["x"], aa["y"], aa["z"]))
    if not asec_coords:
        raise SystemExit(f"No ASEC atom coordinates in {pdb}")

    # Group atoms by residue
    resid_to_atoms = {}
    for a in atoms:
        key = (a["chain"], a["resid"], a["resname"])
        resid_to_atoms.setdefault(key, []).append(a)

    # Compute min distance of each residue (excluding ASEC atoms) to ASEC
    neighbor_residues = {}   # key -> list of atoms
    resid_mindist2    = {}   # key -> min squared distance to ASEC

    for key, ratoms in resid_to_atoms.items():
        # skip residues that are entirely ASEC
        if all(a["serial"] in asec_atoms for a in ratoms):
            continue

        mind2 = float("inf")
        for a in ratoms:
            if a["serial"] in asec_atoms:
                continue
            for (ax, ay, az) in asec_coords:
                dx = a["x"] - ax
                dy = a["y"] - ay
                dz = a["z"] - az
                d2 = dx*dx + dy*dy + dz*dz
                if d2 < mind2:
                    mind2 = d2

        if mind2 <= cutoff*cutoff:
            neighbor_residues[key] = ratoms
            resid_mindist2[key] = mind2

    total_atoms = sum(len(v) for v in neighbor_residues.values())
    print(f"{os.path.basename(pdb)}: {len(neighbor_residues)} neighbor residues, {total_atoms} atoms")

    frames_data.append(dict(
        pdb=pdb,
        neighbor_residues=neighbor_residues,
        resid_mindist2=resid_mindist2
    ))

# --- 4. Per-frame maximal atom budget using full residues ---
per_frame_budgets = []

for fd in frames_data:
    neighbor_residues = fd["neighbor_residues"]
    resid_mindist2    = fd["resid_mindist2"]

    sorted_res = sorted(neighbor_residues.keys(), key=lambda k: resid_mindist2[k])

    atom_count = 0
    for key in sorted_res:
        n_atoms = len(neighbor_residues[key])
        atom_count += n_atoms
    per_frame_budgets.append(atom_count)

target_atoms = min(per_frame_budgets)
if target_atoms == 0:
    raise SystemExit("At least one frame has zero neighbor atoms; cannot equalize.")

print(f"Global target neighbor atoms (full residues): {target_atoms}")

# --- 5. For each frame, pick nearest full residues up to target_atoms ---
for fd in frames_data:
    pdb = fd["pdb"]
    neighbor_residues = fd["neighbor_residues"]
    resid_mindist2    = fd["resid_mindist2"]

    sorted_res = sorted(neighbor_residues.keys(), key=lambda k: resid_mindist2[k])

    chosen_keys = []
    atom_count = 0

    for key in sorted_res:
        n_atoms = len(neighbor_residues[key])
        if atom_count + n_atoms > target_atoms and chosen_keys:
            break
        chosen_keys.append(key)
        atom_count += n_atoms
        if atom_count >= target_atoms:
            break

    out_atoms = []
    for key in chosen_keys:
        out_atoms.extend(neighbor_residues[key])

    # Renumber atom serials
    for i, a in enumerate(out_atoms, start=1):
        line = a["line"]
        a["out_line"] = f"{line[:6]}{i:5d}{line[11:]}"

    out_name = os.path.join(out_dir, os.path.basename(pdb))
    with open(out_name, "w") as f:
        f.write(f"REMARK neighbors of ASEC within {cutoff} A, full residues, {len(out_atoms)} atoms\n")
        for a in out_atoms:
            f.write(a["out_line"])
        f.write("END\n")

    print(f"Wrote {out_name}: {len(chosen_keys)} residues, {len(out_atoms)} atoms")

print("Done: all PDBs in neighbors_eq_fullres have identical atom counts and complete residues.")

