#!/usr/bin/env bash
module add gromacs/2019.4/gcc
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate base
sd=/home/ibsbnicigt/nikolaev_d/script_library/4-solvents/
set -euo pipefail

#----------Read Infos.dat-------------------------------------

if [[ ! -f Infos.dat ]]; then
    echo "ERROR: Infos.dat not found in current directory." >&2
    exit 1
fi

declare Solvent Project Resname Boxsize Ions Fixed Posre

while read -r key value rest; do
    case "$key" in
        Solvent)  Solvent="$value" ;;
        Project)  Project="$value" ;;
        Resname)  Resname="$value" ;;
        Boxsize)  Boxsize="$value" ;;
        Ions)
            Ions="$value $rest"
            ;;
        Fixed)    Fixed="$value" ;;
        Posre)    Posre="$value" ;;
        Temperature)    Temperature="$value" ;;
    esac
done < Infos.dat

: "${Project:?Project not set in Infos.dat}"
: "${Resname:?Resname not set in Infos.dat}"
: "${Boxsize:?Boxsize not set in Infos.dat}"
: "${Ions:?Ions not set in Infos.dat}"
: "${Fixed:?Fixed not set in Infos.dat}"
: "${Posre:?Posre not set in Infos.dat}"

# If Solvent is empty, default to Water
if [[ -z "${Solvent:-}" ]]; then
    Solvent="Water"
fi

echo "Solvent = $Solvent"
echo "Project = $Project"
echo "Resname = $Resname"
echo "Boxsize = $Boxsize"
echo "Ions    = $Ions"
echo "Fixed   = $Fixed"
echo "Posre   = $Posre"
echo "Temperature = $Temperature"



gmx grompp -f npt_asec.mdp -c ${Project}_avg_dist.gro  -p topol.top -o npt_asec.tpr -n index.ndx -r ${Project}_avg_dist.gro -maxwarn 1
mpiexec --bind-to none mdrun_gpu_mpi -deffnm npt_asec

