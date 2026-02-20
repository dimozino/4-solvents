#!/usr/bin/env bash
# check_convergence.sh: assess MD convergence from rmsd_npt.xvg and density_npt.xvg
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base
set -euo pipefail

RMSD_XVG="${1:-rmsd_npt.xvg}"
DENS_XVG="${2:-density_npt.xvg}"

python3 - << 'EOF'
import sys
import math

def read_xvg(path, col=1):
    t = []
    v = []
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(('#','@')):
                    continue
                parts = line.split()
                try:
                    t.append(float(parts[0]))
                    v.append(float(parts[col]))
                except (ValueError, IndexError):
                    continue
    except FileNotFoundError:
        print(f"ERROR: file not found: {path}", file=sys.stderr)
        sys.exit(1)
    if not t:
        print(f"ERROR: no data in {path}", file=sys.stderr)
        sys.exit(1)
    return t, v

def tail_stats(t, v, tail_frac=0.2):
    n = len(t)
    start = int((1.0 - tail_frac) * n)
    if start < 0: start = 0
    tt = t[start:]
    vv = v[start:]
    n_tail = len(vv)
    mean = sum(vv)/n_tail
    var = sum((x-mean)**2 for x in vv) / n_tail
    std = math.sqrt(var)
    # simple linear trend via endpoints
    slope = (vv[-1] - vv[0]) / (tt[-1] - tt[0]) if tt[-1] != tt[0] else 0.0
    return mean, std, slope

def check_file(label, path, col, tail_frac, std_thresh, slope_thresh):
    t, v = read_xvg(path, col=col)
    mean, std, slope = tail_stats(t, v, tail_frac=tail_frac)
    converged = (std <= std_thresh) and (abs(slope) <= slope_thresh)
    print(f"--- {label} ---")
    print(f"file          : {path}")
    print(f"tail fraction : {tail_frac:.2f}")
    print(f"mean (tail)   : {mean:.6f}")
    print(f"std  (tail)   : {std:.6f}")
    print(f"slope (tail)  : {slope:.6g} per ps")
    print(f"converged?    : {'YES' if converged else 'NO'}",
          f"(std <= {std_thresh}, |slope| <= {slope_thresh})")
    print()
    return converged

# User-adjustable parameters
tail_frac      = 0.20   # analyze last 20% of trajectory
rmsd_std_th    = 0.05   # e.g. 0.05 nm
rmsd_slope_th  = 1e-4   # nm per ps
dens_std_th    = 5.0    # kg/m^3 (or whatever units your density_npt.xvg has)
dens_slope_th  = 1e-2   # units per ps

rmsd_path = sys.argv[1] if len(sys.argv) > 1 else "rmsd_npt.xvg"
dens_path = sys.argv[2] if len(sys.argv) > 2 else "density_npt.xvg"

ok_rmsd = check_file("RMSD", rmsd_path, col=1,
                     tail_frac=tail_frac,
                     std_thresh=rmsd_std_th,
                     slope_thresh=rmsd_slope_th)

ok_dens = check_file("Density", dens_path, col=1,
                     tail_frac=tail_frac,
                     std_thresh=dens_std_th,
                     slope_thresh=dens_slope_th)

if ok_rmsd and ok_dens:
    print("Overall convergence: YES (RMSD and density tails are stable)")
    sys.exit(0)
else:
    print("Overall convergence: NO (see details above)")
    sys.exit(1)
EOF

