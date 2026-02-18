#!/usr/bin/env python3
import re
from collections import defaultdict

infos_file   = "Infos.dat"
gro_in       = "NA_avg_dist.gro"
gro_out      = "NA_avg_dist_merged.gro"
top_in       = "topol.top"
top_out      = "topol_merged.top"

# ---------- 1. Parse Infos.dat ASEC line ----------

with open(infos_file) as f:
    infos_txt = f.read()

# Find the ASEC line
asec_line = None
for line in infos_txt.splitlines():
    if line.strip().startswith("ASEC"):
        asec_line = line.strip()
        break

if not asec_line:
    raise SystemExit("No 'ASEC' line found in Infos.dat")

# Extract all resname/resid selectors from the ASEC line.
# Accept patterns like: resname UNK, resid 889, joined by 'or'
resname_patterns = re.findall(r"resname\s+(\S+)", asec_line)
resid_patterns   = re.findall(r"resid\s+(\d+)", asec_line)
resids_from_infos = {int(x) for x in resid_patterns}
resnames_from_infos = set(resname_patterns)

print("ASEC selectors in Infos.dat:")
print(f"  resname: {sorted(resnames_from_infos)}")
print(f"  resid:   {sorted(resids_from_infos)}")

# ---------- 2. Parse GRO and build residue-level info ----------

with open(gro_in) as f:
    lines = f.readlines()

title     = lines[0]
natoms    = int(lines[1].strip())
atom_lines = lines[2:2+natoms]
box_line  = lines[2+natoms]

def parse_gro_atom(line):
    resnr   = int(line[0:5])
    resname = line[5:10].strip()
    atom    = line[10:15].strip()
    atomnr  = int(line[15:20])
    coords  = line[20:].rstrip("\n")
    return resnr, resname, atom, atomnr, coords

def format_gro_atom(resnr, resname, atom, atomnr, coords):
    return f"{resnr:5d}{resname:>5s}{atom:>5s}{atomnr:5d}{coords}\n"

parsed = [parse_gro_atom(l) for l in atom_lines]

# Group by (resnr, resname) to define residues
residue_atoms = defaultdict(list)  # (resnr, resname) -> list of atoms
for rec in parsed:
    resnr, resname, atom, atomnr, coords = rec
    residue_atoms[(resnr, resname)].append(rec)

# Identify the "base" UNK residue we merge into
unk_resids = sorted(r for (r, rn) in residue_atoms.keys() if rn == "UNK")
if not unk_resids:
    raise SystemExit("No UNK residue found in GRO (needed as the base for merging).")
unk_resid = unk_resids[0]
base_key = (unk_resid, "UNK")
print(f"Base UNK residue: resid={unk_resid}")

# Determine which residues are selected by ASEC criteria
selected_keys = set()

for (resnr, resname) in residue_atoms.keys():
    # Selected if matching any resname from Infos.dat
    if resname in resnames_from_infos:
        selected_keys.add((resnr, resname))
    # Selected if matching any explicit resid from Infos.dat
    if resnr in resids_from_infos:
        selected_keys.add((resnr, resname))

if base_key not in selected_keys:
    # Ensure the base UNK residue is always included in the merged set
    selected_keys.add(base_key)

print("Residues selected for merging (resid, resname):")
for k in sorted(selected_keys):
    print("  ", k)

# ---------- 3. Build merged residue + rest, then renumber atoms ----------

# 3.1. Build merged UNK residue: all selected atoms go into UNK with resid=unk_resid, resname=UNK
merged_unk_atoms = []
for key in sorted(selected_keys):
    for (resnr, resname, atom, atomnr, coords) in residue_atoms[key]:
        merged_unk_atoms.append((unk_resid, "UNK", atom, atomnr, coords))

# 3.2. All other residues remain unchanged
other_atoms = []
for key, atom_list in residue_atoms.items():
    if key in selected_keys:
        continue
    other_atoms.extend(atom_list)

# 3.3. Concatenate: merged UNK (first), then the rest
new_atoms = merged_unk_atoms + other_atoms

# 3.4. Renumber atom indices 1..N; residue ids/names remain as set above
renum_atoms = []
for i, (resnr, resname, atom, atomnr, coords) in enumerate(new_atoms, start=1):
    renum_atoms.append((resnr, resname, atom, i, coords))

with open(gro_out, "w") as f:
    f.write(title)
    f.write(f"{natoms:5d}\n")
    for rec in renum_atoms:
        f.write(format_gro_atom(*rec))
    f.write(box_line)

print(f"Wrote merged GRO to {gro_out}")

# ---------- 4. Modify topology: UNK unchanged, others decremented ----------

with open(top_in) as f:
    top_lines = f.readlines()

in_mols = False
new_top = []
# Count how many residues of each resname (other than UNK) were merged
merged_resname_counts = defaultdict(int)
for (resnr, resname) in selected_keys:
    if resname != "UNK":
        merged_resname_counts[resname] += 1

print("Residue name counts to subtract from topology:")
for rn, cnt in merged_resname_counts.items():
    print(f"  {rn}: {cnt}")

for ln in top_lines:
    stripped = ln.strip()
    if stripped.lower().startswith("[ molecules ]"):
        in_mols = True
        new_top.append(ln)
        continue
    if in_mols and stripped.startswith("[") and "]" in stripped:
        in_mols = False
        new_top.append(ln)
        continue

    if in_mols and stripped and not stripped.startswith(";"):
        parts = stripped.split()
        if len(parts) >= 2:
            name, count_s = parts[0], parts[1]
            try:
                count = int(count_s)
            except ValueError:
                new_top.append(ln)
                continue

            if name == "UNK":
                # UNK count unchanged
                new_top.append(ln)
                continue

            if name in merged_resname_counts:
                new_count = count - merged_resname_counts[name]
                if new_count < 0:
                    raise SystemExit(
                        f"Topology would give negative count for {name} "
                        f"({count} - {merged_resname_counts[name]})."
                    )
                new_top.append(f"{name:8s}{new_count:d}\n")
                continue

    new_top.append(ln)

with open(top_out, "w") as f:
    f.writelines(new_top)

print(f"Wrote merged topology to {top_out}")

