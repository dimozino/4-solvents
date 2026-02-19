# General Protocol

1. Get the .xyz, ff.par and mol.psf for your system (solute)
2. Copy the Infos.dat and required files for solvent/ions (e.g. THF.itp THF.gro CL.itp) from ~/script_library/4-solvents
3. Modify Infos.dat
4. Copy and modify topol.top from ~/script_library/4-solvents
5. python ~/script_library/4-solvents/xyz_psf_to_pdb.py $x.xyz mol.psf $x.pdb
6. python ~/script_library/4-solvents/charmm2gmx_itp.py mol.psf ff.par conf_UNK.itp
7. ~/script_library/4-solvents/solvent.sh
8. Wait for end of MD, check equilibration
9. Modify and run extract_average.sh
ASEC starts here
10. Select the solute group, check it is written correctly in Infos.dat
11. Run asec_gensnapshots.sh
12. When finished, run asec_extract.sh
13. Finally, if several solvent boxes for ASEC are given: equalize_asecs.sh from the general folder (inside which folders for all solutes are located).
