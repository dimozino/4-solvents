#!/usr/bin/env bash
module add gromacs/2019.4/gcc
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"

sd=/home/ibsbnicigt/nikolaev_d/script_library/4-solvents/
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base



# 1) Capture full check output (stderr too)
gmx check -f npt_asec.xtc 2> check.log

# 2) Compute dt for ~100 frames
cat > get_dt.py << 'EOF'
import re, math
log = open("check.log").read()
m = re.search(r"Last frame\s+\d+\s+time\s+(\d+\.?\d*)", log)
if not m:
    raise SystemExit("No 'Last frame ... time T' found in check.log")
t_end = float(m.group(1))   # ps
dt = t_end / 100.0
print(int(round(dt)))
EOF

dt=$(python3 get_dt.py)

# 3) Extract frames, each to its own PDB
mkdir -p snapshots
echo System | gmx trjconv -s npt_asec.tpr -f npt_asec.xtc -o snapshots/frame_.pdb -dt "$dt" -sep


cp $sd/rename_first_NUMA_to_UNK.py .
python3 rename_first_NUMA_to_UNK.py 
rm -r snapshots
mv snapshots_renamed/ snapshots

cp $sd/asec_extract.py .
python3 asec_extract.py 
rm neighbors_eq_fullres/frame_0.pdb 
cp $sd/THF_charges.txt .
cp $sd/NA_charges.txt .


cp $sd/asec_merge.py .
python3 asec_merge.py 

mv merged_frames.dat field_asec.pc

cp $sd/qmpart_extract.py .
python3 qmpart_extract.py

