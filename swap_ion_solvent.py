#!/usr/bin/env python3
"""
Swap one ion (anywhere) with one solvent molecule located within a cutoff (Å) of a solute.

Example (UNK + THF + ions):
  python swap_ion_solvent.py em.gro out.gro \
      --solute "resname UNK" \
      --solvent "resname THF" \
      --ions "resname NA or resname CL" \
      --cutoff 5.0
"""

import argparse
import sys
try:
    import MDAnalysis as mda
except ImportError:
    print("ERROR: MDAnalysis not installed. Install with:")
    print("  conda install -c conda-forge mdanalysis")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Swap ion with solvent near solute"
    )
    parser.add_argument("infile", help="Input structure (GRO/PDB)")
    parser.add_argument("outfile", help="Output structure (GRO/PDB)")
    parser.add_argument("--solute", required=True,
                        help="Solute selection (e.g., 'resname UNK')")
    parser.add_argument("--solvent", required=True,
                        help="Solvent selection (e.g., 'resname THF')")
    parser.add_argument("--ions", required=True,
                        help="Ion selection (e.g., 'resname NA or resname CL')")
    parser.add_argument("--cutoff", type=float, default=5.0,
                        help="Distance cutoff (Å) for nearby solvent (default: 5.0)")
    args = parser.parse_args()

    u = mda.Universe(args.infile)
    solute = u.select_atoms(args.solute)
    if len(solute) == 0:
        print(f"WARNING: No atoms found for solute selection: {args.solute}")
        print("Writing input structure unchanged.")
        u.atoms.write(args.outfile)
        return

    # Find solvent within cutoff of solute
    solvent_all = u.select_atoms(args.solvent)
    nearby_solvent = solvent_all.select_atoms(f"around {args.cutoff} global group solute",
                                               solute=solute)
    
    # Get full residues for nearby solvent
    nearby_solvent_res = u.select_atoms(f"same residue as group nearby",
                                         nearby=nearby_solvent)
    
    # Get ions (anywhere in the system)
    ions = u.select_atoms(args.ions)

    if len(nearby_solvent_res) == 0 or len(ions) == 0:
        print(f"WARNING: No nearby solvent or ions. Nearby solvent: {len(nearby_solvent_res)}, Ions: {len(ions)}")
        print("Writing input structure unchanged.")
        u.atoms.write(args.outfile)
        return

    # Get residues
    solvent_residues = list(set(nearby_solvent_res.residues))
    ion_residues = list(set(ions.residues))

    # Swap as many as possible
    n_swaps = min(len(solvent_residues), len(ion_residues))
    print(f"Swapping {n_swaps} ion(s) with solvent near solute (cutoff={args.cutoff} Å)")

    for i in range(n_swaps):
        sol_res = solvent_residues[i]
        ion_res = ion_residues[i]
        
        # Compute centers of mass
        sol_com = sol_res.atoms.center_of_mass()
        ion_com = ion_res.atoms.center_of_mass()
        
        # Swap positions by translating each residue
        shift_sol = ion_com - sol_com
        shift_ion = sol_com - ion_com
        sol_res.atoms.positions += shift_sol
        ion_res.atoms.positions += shift_ion

    # Write updated structure
    u.atoms.write(args.outfile)
    print(f"Wrote {args.outfile} with {n_swaps} swapped ion(s)")

if __name__ == "__main__":
    main()
