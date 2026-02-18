#!/bin/bash

Project=`grep "Project" Infos.dat | awk '{ print $2 }'`
sd=/home/ibsbnicigt/nikolaev_d/script_library/4-solvents/
module add gromacs/2019.4/gcc
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"

gmx distance -s npt.tpr -f npt.xtc -select 'resid 1 and name C03 plus resname CL and name CL' -oall dist_UNK_CL.xvg

cp $sd/find_time_closest_3000.py .

source ~/miniconda3/etc/profile.d/conda.sh
conda activate base
# Get time (ps) of frame with distance closest to average
tavg=$(python3 find_time_closest_3000.py dist_UNK_CL.xvg 1 2.0)

# Extract that frame from the trajectory
gmx trjconv -s npt.tpr -f npt.xtc -dump "$tavg" -o ${Project}_avg_dist.gro << EOF 
0 
EOF
