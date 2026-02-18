#!/usr/bin/env python3

import sys
import re

if len(sys.argv) < 2:
    print("Usage: merge_itp.py conf_UNK.itp [extra1.itp extra2.itp ...]")
    sys.exit(1)

unk_itp = sys.argv[1]
extra_itps = sys.argv[2:]
out_itp = "output.itp"

def read_lines(fname):
    with open(fname) as f:
        return f.readlines()

def extract_moleculetype_block(lines, name):
    """Return (start, end) indices of the moleculetype block with given name."""
    start = None
    for i, ln in enumerate(lines):
        if ln.strip().lower().startswith("[ moleculetype ]"):
            j = i + 1
            while j < len(lines) and (lines[j].strip() == "" or lines[j].lstrip().startswith(";")):
                j += 1
            if j < len(lines) and lines[j].split():
                if lines[j].split()[0] == name:
                    start = i
                    break
    if start is None:
        raise RuntimeError(f"moleculetype {name} not found")
    end = len(lines)
    for k in range(start + 1, len(lines)):
        if lines[k].strip().lower().startswith("[ moleculetype ]"):
            end = k
            break
    return start, end

def count_atoms_in_block(block_lines):
    """Count number of atoms in [ atoms ] section of given block."""
    n = 0
    in_atoms = False
    for ln in block_lines:
        stripped = ln.strip()
        if stripped.lower().startswith("[ atoms ]"):
            in_atoms = True
            continue
        if in_atoms and stripped.startswith("[") and "]" in stripped and not stripped.lower().startswith("[ atoms ]"):
            break
        if in_atoms:
            if stripped == "" or stripped.startswith(";"):
                continue
            parts = stripped.split()
            try:
                int(parts[0])
                n += 1
            except (ValueError, IndexError):
                pass
    return n

def offset_indices_in_block(block_lines, offset):
    """
    Offset atom indices inside a moleculetype block by 'offset' for sections:
    [ atoms ], [ bonds ], [ angles ], [ dihedrals ], [ pairs ].
    """
    out = []
    current_sec = None
    sec_re = re.compile(r"^\s*\[\s*(\w+)\s*\]")

    for ln in block_lines:
        m = sec_re.match(ln)
        if m:
            current_sec = m.group(1).lower()
            out.append(ln)
            continue

        if current_sec in ("atoms", "bonds", "angles", "dihedrals", "pairs"):
            s = ln.strip()
            if s == "" or s.startswith(";"):
                out.append(ln)
                continue
            parts = ln.split()
            if current_sec == "atoms":
                try:
                    parts[0] = f"{int(parts[0]) + offset:6d}"
                except ValueError:
                    pass
            elif current_sec == "bonds":
                for i in range(2):
                    try:
                        parts[i] = f"{int(parts[i]) + offset:6d}"
                    except ValueError:
                        pass
            elif current_sec == "angles":
                for i in range(3):
                    try:
                        parts[i] = f"{int(parts[i]) + offset:6d}"
                    except ValueError:
                        pass
            elif current_sec == "pairs":
                for i in range(2):
                    try:
                        parts[i] = f"{int(parts[i]) + offset:6d}"
                    except ValueError:
                        pass
            elif current_sec == "dihedrals":
                for i in range(4):
                    try:
                        parts[i] = f"{int(parts[i]) + offset:6d}"
                    except ValueError:
                        pass
            out.append(" ".join(parts) + "\n")
        else:
            out.append(ln)

    return out

def strip_mt_header(block_lines):
    """Strip [ moleculetype ] header and name from a block."""
    out = []
    skipping = False
    mt_seen = False
    for ln in block_lines:
        if ln.strip().lower().startswith("[ moleculetype ]"):
            skipping = True
            mt_seen = True
            continue
        if skipping:
            if ln.strip() and not ln.lstrip().startswith(";"):
                skipping = False
                continue
            continue
        out.append(ln)
    if not mt_seen:
        return block_lines
    return out

# Read UNK file and extract its moleculetype
unk_lines = read_lines(unk_itp)
u_start, u_end = extract_moleculetype_block(unk_lines, "UNK")
unk_block = unk_lines[u_start:u_end]

n_total = count_atoms_in_block(unk_block)

# Prepare UNK block header (we keep name "UNK" as merged type)
merged_mt_name = "UNK"
out_unk_block = []
seen_mt = False
for ln in unk_block:
    if ln.strip().lower().startswith("[ moleculetype ]"):
        out_unk_block.append(ln)
        seen_mt = True
        continue
    if seen_mt and ln.strip() and not ln.lstrip().startswith(";"):
        parts = ln.split()
        parts[0] = merged_mt_name
        out_unk_block.append(" ".join(parts) + "\n")
        seen_mt = False
    else:
        out_unk_block.append(ln)

# Process each extra itp: extract its first moleculetype block, offset, strip header, append
extra_bodies = []
for extra in extra_itps:
    lines = read_lines(extra)
    # Guess moleculetype name: first [ moleculetype ] then following name
    mt_name = None
    for i, ln in enumerate(lines):
        if ln.strip().lower().startswith("[ moleculetype ]"):
            j = i + 1
            while j < len(lines) and (lines[j].strip() == "" or lines[j].lstrip().startswith(";")):
                j += 1
            if j < len(lines) and lines[j].split():
                mt_name = lines[j].split()[0]
            break
    if mt_name is None:
        raise RuntimeError(f"No [ moleculetype ] found in {extra}")
    s, e = extract_moleculetype_block(lines, mt_name)
    block = lines[s:e]
    n_atoms = count_atoms_in_block(block)
    block_off = offset_indices_in_block(block, n_total)
    body = strip_mt_header(block_off)
    extra_bodies.append(body)
    n_total += n_atoms

# Build final output
out_lines = []
out_lines.extend(unk_lines[:u_start])
out_lines.extend(out_unk_block)
for body in extra_bodies:
    out_lines.extend(body)
out_lines.extend(unk_lines[u_end:])

with open(out_itp, "w") as f:
    f.writelines(out_lines)

print(f"Wrote merged itp to {out_itp}")

