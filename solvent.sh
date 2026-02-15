#!/usr/bin/env bash
set -euo pipefail

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
        Ions)     Ions="$value" ;;
        Fixed)    Fixed="$value" ;;
        Posre)    Posre="$value" ;;
    esac
done < Infos.dat

for var in Solvent Project Resname Boxsize Ions Fixed Posre; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: $var not set in Infos.dat" >&2
        exit 1
    fi
done

echo "Building solvated box for $Resname in $Solvent"

# 1. Build box with solvent
gmx solvate \
  -cp conf_${Resname}.gro \
  -cs ${Solvent}.gro \
  -box ${Boxsize} ${Boxsize} ${Boxsize} \
  -o solv.gro

# 2. Prepare topology
cp topol_${Solvent}.top topol.top

# 3. Generate ion placement TPR
gmx grompp \
  -f ions.mdp \
  -c solv.gro \
  -p topol.top \
  -o ions.tpr \
  -maxwarn 1

# 4. Add ions to neutralize
echo "SOL" | gmx genion \
  -s ions.tpr \
  -o solv_ions.gro \
  -p topol.top \
  -pname NA -nname CL \
  -neutral \
  -conc ${Ions}

# 5. Move ions near solute using Python script
if [[ -f swap_ion_solvent.py ]]; then
    python swap_ion_solvent.py solv_ions.gro solv_ions_swapped.gro \
        --solute "resname ${Resname}" \
        --solvent "resname ${Solvent}" \
        --ions "resname NA or resname CL" \
        --cutoff 5.0
    mv solv_ions_swapped.gro solv_ions.gro
fi

# 6. Create index groups
{
    echo "q"
} | gmx make_ndx -f solv_ions.gro -o index.ndx

# Add non-solvent group
{
    echo "! r ${Solvent}"
    echo "name 13 non-${Solvent}"
    echo "q"
} | gmx make_ndx -f solv_ions.gro -n index.ndx -o index.ndx

# 7. Generate position restraints
echo "${Fixed}" | gmx genrestr \
  -f solv_ions.gro \
  -n index.ndx \
  -o posre_${Fixed}.itp \
  -fc ${Posre} ${Posre} ${Posre}

# 8. Update topology to include position restraints
if ! grep -q "posre_${Fixed}.itp" topol.top; then
    sed -i "/\[ moleculetype \]/a\
#ifdef POSRES\n#include \"posre_${Fixed}.itp\"\n#endif" topol.top
fi

# 9. Energy minimization
gmx grompp \
  -f em.mdp \
  -c solv_ions.gro \
  -p topol.top \
  -n index.ndx \
  -o em.tpr \
  -maxwarn 1

gmx mdrun -v -deffnm em

# 10. Submit equilibration job
bash solvent_run.sh

echo "Solvent box setup complete. Check em.gro and wait for equilibration to finish."
