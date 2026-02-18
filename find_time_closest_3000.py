#!/usr/bin/env python3

import numpy as np
import sys

def read_xvg(filename, column=1):
    times = []
    values = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('#', '@')):
                continue
            parts = line.split()
            times.append(float(parts[0]))      # time (ps)
            values.append(float(parts[column]))
    return np.array(times), np.array(values)

def main(xvg_file, column=1, window_ns=3.0):
    times, dists = read_xvg(xvg_file, column=column)

    if times.size == 0:
        print("No data points found", file=sys.stderr)
        sys.exit(1)

    # Last window_ns ns (in ps)
    window_ps = window_ns * 1000.0
    t_max = times.max()
    t_min_window = t_max - window_ps

    # Use only frames in that window
    mask = times >= t_min_window
    if not np.any(mask):
        # If trajectory shorter than window: use all data
        mask = np.ones_like(times, dtype=bool)

    times_win = times[mask]
    dists_win = dists[mask]

    avg = dists_win.mean()

    # Find closest frame WITHIN the window
    idx_local = np.argmin(np.abs(dists_win - avg))
    snapshot_time = times_win[idx_local]

    # Print absolute snapshot time (ps)
    print(f"{snapshot_time:.6f}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} dist.xvg [column_index] [window_ns]", file=sys.stderr)
        sys.exit(1)

    xvg_file = sys.argv[1]
    column = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    window_ns = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0

    main(xvg_file, column, window_ns)

