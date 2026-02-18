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

declare Solvent Project Resname Boxsize Ions Fixed Posre ITP_FILES

while read -r key value rest; do
  case "$key" in
    Solvent)    Solvent="$value" ;;
    Project)    Project="$value" ;;
    Resname)    Resname="$value" ;;
    Boxsize)    Boxsize="$value" ;;
    Ions)       Ions="$value $rest" ;;
    Fixed)      Fixed="$value" ;;
    Posre)      Posre="$value" ;;
    Temperature) Temperature="$value" ;;
    ITP)        ITP_FILES="$value $rest" ;;
  esac
done < Infos.dat

: "${Project:?Project not set in Infos.dat}"
: "${Resname:?Resname not set in Infos.dat}"
: "${Boxsize:?Boxsize not set in Infos.dat}"
: "${Ions:?Ions not set in Infos.dat}"
: "${Fixed:?Fixed not set in Infos.dat}"
: "${Posre:?Posre not set in Infos.dat}"
: "${ITP_FILES:?ITP not set in Infos.dat}"

# If Solvent is empty, default to Water
if [[ -z "${Solvent:-}" ]]; then
  Solvent="Water"
fi

echo "Solvent = $Solvent"
echo "Project = $Project"
echo "Resname = $Resname"
echo "Boxsize = $Boxsize"
echo "Ions = $Ions"
echo "Fixed = $Fixed"
echo "Posre = $Posre"
echo "Temperature = $Temperature"
echo "ITP files = $ITP_FILES"

mkdir -p asec
cd asec

# Copy needed files
cp ../${Project}_avg_dist.gro .
cp ../topol.top .
cp ../Infos.dat .



# Copy all ITPs mentioned in Infos.dat
for itp in $ITP_FILES; do
  cp ../"$itp" .
done

cp ../${Solvent}.itp .
cp $sd/merge_itp.py .
cp $sd/merge_3c.py .

# Merge all specified ITPs  into a single conf_${Resname}.itp
# First ITP in the list must be conf_${Resname}.itp
python3 merge_itp.py $ITP_FILES 

# New script always writes output.itp; rename to conf_${Resname}.itp
mv output.itp conf_${Resname}.itp

cp ${Project}_avg_dist.gro NA_avg_dist.gro
python3 merge_3c.py
mv topol_merged.top topol.top
mv NA_avg_dist_merged.gro ${Project}_avg_dist.gro

# ===== 5. Create index groups properly ===================================

echo "q" | gmx make_ndx -f ${Project}_avg_dist.gro -o index_temp.ndx 2>&1 | tee make_ndx_output.txt

NGROUPS=$(grep -c "^ *[0-9]" make_ndx_output.txt || echo 0)
echo "Number of default groups: $NGROUPS"

rm -f index_temp.ndx make_ndx_output.txt

gmx make_ndx -f ${Project}_avg_dist.gro -o index.ndx << EOF
r ${Fixed}
!r ${Solvent}
name ${NGROUPS} ${Fixed}_posres
name $((NGROUPS + 1)) non-${Solvent}
q
EOF

echo "Created groups: '${Fixed}_posres' and 'non-${Solvent}'"

if ! grep -q "\[ ${Fixed}_posres \]" index.ndx; then
  echo "ERROR: Failed to create group '${Fixed}_posres' in index.ndx" >&2
  echo "Available groups:"
  grep "^\[" index.ndx
  exit 1
fi

if ! grep -q "\[ non-${Solvent} \]" index.ndx; then
  echo "ERROR: Failed to create group 'non-${Solvent}' in index.ndx" >&2
  echo "Available groups:"
  grep "^\[" index.ndx
  exit 1
fi

gmx genrestr -f ${Project}_avg_dist.gro -n index.ndx -o posre_${Fixed}.itp -fc "$Posre" "$Posre" "$Posre" << EOF
${Fixed}_posres
EOF

####------ RUN ASEC MD -------####

cp $sd/npt_asec.mdp .
cp $sd/asec_mdrun.sh .

bash ~/bin/mBASH asec_mdrun -Nnod tornado -Npr 24 -Time 30

