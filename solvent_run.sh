#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f Infos.dat ]]; then
    echo "ERROR: Infos.dat not found in current directory." >&2
    exit 1
fi

declare Solvent Project

while read -r key value rest; do
    case "$key" in
        Solvent)  Solvent="$value" ;;
        Project)  Project="$value" ;;
    esac
done < Infos.dat

for var in Solvent Project; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: $var not set in Infos.dat" >&2
        exit 1
    fi
done

echo "Running equilibration for $Project"

# NVT equilibration
gmx grompp \
  -f nvt.mdp \
  -c em.gro \
  -r em.gro \
  -p topol.top \
  -n index.ndx \
  -o nvt.tpr \
  -maxwarn 1

gmx mdrun -v -deffnm nvt

# NPT equilibration  
gmx grompp \
  -f npt.mdp \
  -c nvt.gro \
  -r nvt.gro \
  -t nvt.cpt \
  -p topol.top \
  -n index.ndx \
  -o npt.tpr \
  -maxwarn 1

gmx mdrun -v -deffnm npt

# Analysis: RMSD
echo "Backbone" | gmx rms \
  -s npt.tpr \
  -f npt.xtc \
  -o rmsd_npt.xvg \
  -tu ns

# Analysis: Density
echo "System" | gmx energy \
  -f npt.edr \
  -o density_npt.xvg <<< "Density"

echo "Equilibration complete. Output files: npt.gro, npt.cpt, rmsd_npt.xvg, density_npt.xvg"
