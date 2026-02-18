#!/usr/bin/env python3
import sys
import math
from collections import defaultdict

# ---------- PSF PARSING ----------

def read_psf(psf_path):
    with open(psf_path, "r") as f:
        lines = f.readlines()

    atoms = []
    bonds = []
    angles = []
    dihedrals = []
    impropers = []

    i = 0
    n = len(lines)

    def skip_to(tag):
        nonlocal i
        while i < n:
            if tag in lines[i]:
                return
            i += 1

    # Atoms
    skip_to("!NATOM")
    natom = int(lines[i].split()[0])
    i += 1
    for _ in range(natom):
        fields = lines[i].split()
        idx = int(fields[0])
        segid = fields[1]
        resid = int(fields[2])
        resname = fields[3]
        atomname = fields[4]    # this is the label like C0, O5, etc.
        atomtype = fields[5]
        charge = float(fields[6])
        mass = float(fields[7])
        atoms.append(
            dict(
                index=idx,
                segid=segid,
                resid=resid,
                resname=resname,
                atomname=atomname,
                atomtype=atomtype,
                charge=charge,
                mass=mass,
            )
        )
        i += 1

    # Bonds
    skip_to("!NBOND")
    nbond = int(lines[i].split()[0])
    i += 1
    read_pairs = 0
    while read_pairs < nbond:
        nums = [int(x) for x in lines[i].split()]
        for k in range(0, len(nums), 2):
            a1, a2 = nums[k:k+2]
            bonds.append((a1, a2))
            read_pairs += 1
            if read_pairs >= nbond:
                break
        i += 1

    # Angles
    skip_to("!NTHETA")
    ntheta = int(lines[i].split()[0])
    i += 1
    read_triples = 0
    while read_triples < ntheta:
        nums = [int(x) for x in lines[i].split()]
        for k in range(0, len(nums), 3):
            a1, a2, a3 = nums[k:k+3]
            angles.append((a1, a2, a3))
            read_triples += 1
            if read_triples >= ntheta:
                break
        i += 1

    # Dihedrals
    skip_to("!NPHI")
    nphi = int(lines[i].split()[0])
    i += 1
    read_quads = 0
    while read_quads < nphi:
        nums = [int(x) for x in lines[i].split()]
        for k in range(0, len(nums), 4):
            a1, a2, a3, a4 = nums[k:k+4]
            dihedrals.append((a1, a2, a3, a4))
            read_quads += 1
            if read_quads >= nphi:
                break
        i += 1

    # Impropers
    skip_to("!NIMPHI")
    nimphi = int(lines[i].split()[0])
    i += 1
    read_improp = 0
    while read_improp < nimphi:
        nums = [int(x) for x in lines[i].split()]
        for k in range(0, len(nums), 4):
            a1, a2, a3, a4 = nums[k:k+4]
            impropers.append((a1, a2, a3, a4))
            read_improp += 1
            if read_improp >= nimphi:
                break
        i += 1

    return atoms, bonds, angles, dihedrals, impropers


# ---------- ff.par PARSING ----------

def read_par(par_path):
    with open(par_path, "r") as f:
        lines = f.readlines()

    section = None
    bond_params = {}        # (name1,name2) sorted -> (kb, r0)
    angle_params = {}       # (n1,n2,n3) -> (ktheta, theta0)
    dihed_params = defaultdict(list)  # (n1,n2,n3,n4) -> list of (k, n, delta)
    nb_params = {}          # name -> (epsilon, sigma)

    def norm_pair(a, b):
        return tuple(sorted((a, b)))

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        if line.startswith("BONDS"):
            section = "BONDS"
            i += 1
            continue
        if line.startswith("THETAS"):
            section = "THETAS"
            i += 1
            continue
        if line.startswith("PHI"):
            section = "PHI"
            i += 1
            continue
        if line.startswith("IMPHI"):
            section = "IMPHI"
            i += 1
            continue
        if line.startswith("NONBONDED"):
            section = "NONBONDED"
            # skip header line with CUTNB, etc.
            i += 1
            continue
        if line.startswith("END"):
            break

        fields = line.split()

        if section == "BONDS":
            # C0 C1 451.145 1.452
            if len(fields) >= 4:
                a, b, kb, r0 = fields[:4]
                kb = float(kb) * 836.8   # kcal/mol/Å^2 -> kJ/mol/Å^2
                r0 = float(r0) * 0.1     # Å -> nm
                bond_params[norm_pair(a, b)] = (kb, r0)

        elif section == "THETAS":
            # C0 C1 C3 44.689 128.173
            if len(fields) >= 5:
                a, b, c, ktheta, theta0 = fields[:5]
                ktheta = float(ktheta) * 8.368  # kcal/mol/rad^2 -> kJ
                theta0 = float(theta0)          # degrees
                angle_params[(a, b, c)] = (ktheta, theta0)

        elif section == "PHI":
            # C0 C1 C3 O5 -1.538 1 0.000
            if len(fields) >= 7:
                a, b, c, d = fields[:4]
                k = float(fields[4]) * 4.184    # kcal/mol -> kJ/mol
                n_mult = int(fields[5])
                delta = float(fields[6])        # degrees
                dihed_params[(a, b, c, d)].append((k, n_mult, delta))

        elif section == "NONBONDED":
            # C0 0.000 -0.080 1.964 0.000 -0.040 1.964
            #        eps     Rmin       eps14    Rmin14
            if len(fields) >= 7:
                name = fields[0]
                eps = float(fields[2]) * -4.184      # kcal/mol -> kJ/mol
                rmin = float(fields[3]) * 0.1       # Å -> nm
                # CHARMM: Rmin/2 often, but here the file appears to give
                # full Rmin; approximate sigma from Rmin: Rmin = 2^(1/6)*sigma
                sigma = rmin * 2.0 / (2.0 ** (1.0 / 6.0))
                nb_params[name] = (eps, sigma)

        # IMPHI section can be added similarly if needed.

        i += 1

    return bond_params, angle_params, dihed_params, nb_params


# ---------- ITP WRITER ----------

def write_itp(out_path, atoms, bonds, angles, dihedrals, impropers,
              bond_params, angle_params, dihed_params, nb_params):

    idx2atom = {a["index"]: a for a in atoms}
    resname = atoms[0]["resname"] if atoms else "UNK"

    with open(out_path, "w") as out:
        out.write(";\n; Generated from CHARMM PSF + ff.par\n;\n\n")

        # Atomtypes from NONBONDED (names like C0, O5...)
        out.write("[ atomtypes ]\n")
        out.write("; name  at.num   mass   charge  ptype    sigma      epsilon\n")
        # Guess element and mass from first char
        elem_mass = {
            "H": 1.008, "C": 12.011, "N": 14.007,
            "O": 15.999, "S": 32.060, "P": 30.974,
            "F": 18.998, "C": 12.011
        }
        for name in sorted(nb_params.keys()):
            eps, sigma = nb_params[name]
            elem = name[0]
            mass = elem_mass.get(elem, 12.011)
            out.write(f"{name:6s}  {elem:4s}  {mass:8.4f}   0.000    A  "
                      f"{sigma:11.6e}  {eps:11.6e}\n")
        out.write("\n")

        # Moleculetype
        out.write("[ moleculetype ]\n")
        out.write("; name  nrexcl\n")
        out.write(f"{resname:6s}  3\n\n")

        # Atoms block: use atomname as gromacs atom type,
        # but LJ parameters are defined by atomname via NONBONDED.
        out.write("[ atoms ]\n")
        out.write("; nr  type  resnr  residue  atom  cgnr   charge    mass\n")
        for a in atoms:
            nr = a["index"]
            aname = a["atomname"]
            resnr = a["resid"]
            charge = a["charge"]
            mass = a["mass"]
            cgnr = resnr
            out.write(f"{nr:5d}  {aname:6s}  {resnr:4d}  {resname:6s}  "
                      f"{aname:4s}  {cgnr:4d}  {charge:9.6f}  {mass:8.4f}\n")
        out.write("\n")

        # Bonds
        out.write("[ bonds ]\n")
        out.write("; ai  aj  funct     r0(nm)           kb\n")
        for i, j in bonds:
            n1 = idx2atom[i]["atomname"]
            n2 = idx2atom[j]["atomname"]
            key = tuple(sorted((n1, n2)))
            if key in bond_params:
                kb, r0 = bond_params[key]
            else:
                kb, r0 = 0.0, 0.1
            out.write(f"{i:5d}{j:6d}   1   {r0:9.5f}   {kb:12.3f}\n")
        out.write("\n")

        # Angles
        out.write("[ angles ]\n")
        out.write("; ai  aj  ak  funct   theta(deg)      ktheta\n")
        for i, j, k in angles:
            n1 = idx2atom[i]["atomname"]
            n2 = idx2atom[j]["atomname"]
            n3 = idx2atom[k]["atomname"]
            key = (n1, n2, n3)
            if key not in angle_params:
                key = (n3, n2, n1)
            if key in angle_params:
                ktheta, theta0 = angle_params[key]
            else:
                ktheta, theta0 = 0.0, 109.5
            out.write(f"{i:5d}{j:6d}{k:6d}   1   {theta0:9.3f}   {ktheta:12.3f}\n")
        out.write("\n")

        # Dihedrals (proper) – use function type 9 with multiplicities
        out.write("[ dihedrals ]\n")
        out.write("; proper dihedrals: function 9 (k, n, delta)\n")
        out.write("; ai  aj  ak  al  funct       phi0      k       mult\n")
        used = set()
        for i, j, k, l in dihedrals:
            if (i, j, k, l) in used:
                continue
            n1 = idx2atom[i]["atomname"]
            n2 = idx2atom[j]["atomname"]
            n3 = idx2atom[k]["atomname"]
            n4 = idx2atom[l]["atomname"]
            key = (n1, n2, n3, n4)
            if key not in dihed_params:
                key = (n4, n3, n2, n1)
            terms = dihed_params.get(key, [])
            if not terms:
                # no parameters, still write a zero torsion so list is complete
                out.write(f"{i:5d}{j:6d}{k:6d}{l:6d}   9   0.0   0.0   1\n")
            else:
                for kchi, n_mult, delta in terms:
                    out.write(f"{i:5d}{j:6d}{k:6d}{l:6d}   9   "
                              f"{delta:7.2f}   {kchi:7.3f}   {n_mult:d}\n")
            used.add((i, j, k, l))
        out.write("\n")

        # Impropers (if present, need parameters; here put placeholders)
        if impropers:
            out.write("; impropers (placeholders)\n")
            for i, j, k, l in impropers:
                out.write(f"{i:5d}{j:6d}{k:6d}{l:6d}   4   180.0   0.0\n")
            out.write("\n")


# ---------- MAIN ----------

def main():
    if len(sys.argv) != 4:
        print("Usage: python charmm2gmx_itp.py mol.psf ff.par output.itp")
        sys.exit(1)

    psf_path = sys.argv[1]
    par_path = sys.argv[2]
    out_path = sys.argv[3]

    atoms, bonds, angles, dihedrals, impropers = read_psf(psf_path)
    bond_params, angle_params, dihed_params, nb_params = read_par(par_path)
    write_itp(out_path, atoms, bonds, angles, dihedrals, impropers,
              bond_params, angle_params, dihed_params, nb_params)


if __name__ == "__main__":
    main()

