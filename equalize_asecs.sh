#!/usr/bin/env bash
set -euo pipefail

currdir="$(pwd)"
sd=/home/ibsbnicigt/nikolaev_d/script_library/4-solvents/
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

FOLDER_LIST="folders.txt"

if [[ ! -f "$FOLDER_LIST" ]]; then
    echo "ERROR: $FOLDER_LIST not found" >&2
    exit 1
fi

smallest_number=999999999

echo "=== Step 1: determine smallest_number from frame_1.pdb ==="

while read -r x; do
    [[ -z "$x" ]] && continue
    dir="$x/asec/neighbors_eq_fullres"
    pdb="$dir/frame_1.pdb"

    if [[ ! -f "$pdb" ]]; then
        echo "  [$x] $pdb not found, skipping" >&2
        continue
    fi

    natoms=$(grep -Ec '^(ATOM  |HETATM)' "$pdb" || echo 0)
    echo "  [$x] frame_1.pdb atoms = $natoms"

    if (( natoms > 0 && natoms < smallest_number )); then
        smallest_number=$natoms
    fi
done < "$FOLDER_LIST"

if (( smallest_number == 999999999 )); then
    echo "ERROR: no valid frame_1.pdb found" >&2
    exit 1
fi

echo "Global smallest_number = $smallest_number"
echo "=== Step 2: trim frame_1.pdb ... frame_100.pdb to $smallest_number atoms ==="

# Step 2: use AWK per file, no nested while
while read -r x; do
    [[ -z "$x" ]] && continue
    dir="$x/asec/neighbors_eq_fullres"

    if [[ ! -d "$dir" ]]; then
        echo "  [$x] directory $dir not found, skipping" >&2
        continue
    fi

    echo "  [$x] processing $dir"

    for i in $(seq 1 100); do
        pdb="$dir/frame_${i}.pdb"
        if [[ ! -f "$pdb" ]]; then
            echo "    $pdb not found, skipping"
            continue
        fi

        tmp="${pdb}.tmp"

        awk -v maxatoms="$smallest_number" '
            BEGIN { n=0 }
            /^ATOM  / || /^HETATM/ {
                n++
                if (n <= maxatoms) print
                next
            }
            /^END/ { next }   # skip existing END
            { print }
            END { print "END" }
        ' "$pdb" > "$tmp"

        mv "$tmp" "$pdb"

        new_atoms=$(grep -Ec '^(ATOM  |HETATM)' "$pdb" || echo 0)
        echo "    trimmed $(basename "$pdb") to $new_atoms atoms"
    done

done < "$FOLDER_LIST"

FOLDER_LIST="folders.txt"

while read -r x; do
    [[ -z "$x" ]] && continue

    dir="$x/asec/"
    cd $dir/
    cp $sd/asec_merge.py .
    python3 asec_merge.py

    mv merged_frames.dat field_asec.pc

    cp $sd/qmpart_extract.py .
    python3 qmpart_extract.py
    
    mv field_asec.pc ${x}_field_asec.pc
    mv qmpart.xyz ${x}_qmpart.xyz 

    cd $currdir

done < "$FOLDER_LIST"




echo "Done."



