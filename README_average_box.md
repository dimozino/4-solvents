# Distance-Based Snapshot Extraction Tools

This mini-toolkit identifies a frame in a GROMACS trajectory where a selected distance is closest to its time-averaged value over a specified time window, and then extracts the corresponding structure.

Scripts:

- `find_time_closest_3000.py` – parses a GROMACS `.xvg` time series and prints the time of the frame whose value is closest to the average over the last *N* ns. [file:307]
- `extract_average.sh` – computes a specific distance in `npt.xtc`, calls the Python helper to get the representative time, and extracts the corresponding snapshot as `<Project>_avg_dist.gro`. [file:306]

## 1. `find_time_closest_3000.py`

### Purpose

Given a GROMACS `.xvg` file containing a time series (e.g. a distance or RMSD), this script:

1. Reads a specified data column vs time.
2. Restricts analysis to the last `window_ns` nanoseconds of the trajectory.
3. Computes the average value in that window.
4. Locates the frame whose value is closest to this average.
5. Prints the corresponding **absolute time (ps)** to stdout. [file:307]

This time can then be passed to `gmx trjconv -dump` to extract the representative snapshot. [file:306][file:307]

### Usage

```bash
python3 find_time_closest_3000.py dist.xvg [column_index] [window_ns]
```

- `dist.xvg` – input `.xvg` file from GROMACS (headers starting with `#` or `@` are handled). [file:307]
- `column_index` (optional, default `1`) – zero-based index of the data column to analyze:
    - `0` – first data column after time,
    - `1` – second data column, etc.
For typical 2-column `time value` files, use `1`. [file:307]
- `window_ns` (optional, default `3.0`) – length of the analysis window in nanoseconds, counting backwards from the last frame. [file:307]

Example (use column 1, last 2 ns):

```bash
tavg=$(python3 find_time_closest_3000.py dist_UNK_CL.xvg 1 2.0)
echo "Representative time: $tavg ps"
```


### Behavior

Internally, the script: [file:307]

1. Calls `read_xvg(filename, column)`:
    - Skips lines starting with `#` or `@`.
    - Splits remaining lines on whitespace.
    - Reads:
        - `times.append(float(parts[^0]))`  \# time in ps
        - `values.append(float(parts[column]))`  \# selected observable
    - Returns NumPy arrays `times`, `values`. [file:307]
2. Computes the time window: [file:307]

```python
window_ps = window_ns * 1000.0
t_max = times.max()
t_min_window = t_max - window_ps
mask = times >= t_min_window
```

If no frames fall within the window (trajectory shorter than `window_ns`), it falls back to using all data points. [file:307]
3. Computes the average within the window:

```python
avg = dists_win.mean()
```

4. Finds the index of the value closest to `avg`:

```python
idx_local = np.argmin(np.abs(dists_win - avg))
snapshot_time = times_win[idx_local]
```

5. Prints `snapshot_time` to stdout with 6 decimal places. [file:307]

If no data points are found, it prints an error to stderr and exits with code 1. [file:307]

## 2. `extract_average.sh`

### Purpose

Automate the full workflow for a specific solute–ion distance in an NPT trajectory:

1. Compute the distance between:
    - a reference atom in the solute (here: `resid 1 and name C03`), and
    - all chloride ions (`resname CL and name CL`). [file:306]
2. Analyze the distance time series over the last part of the trajectory using `find_time_closest_3000.py`.
3. Extract the structure at the time where the distance is closest to its average in that window. [file:306]

The final output is `<Project>_avg_dist.gro`, where `Project` is taken from `Infos.dat`. [file:306]

### Usage

Run inside a project directory containing:

- `Infos.dat` (with a `Project` entry).
- `npt.tpr`, `npt.xtc` – equilibrated NPT setup and trajectory.

Then:

```bash
bash extract_average.sh
```


### Steps

1. **Read project name** from `Infos.dat`: [file:306]

```bash
Project=$(grep "Project" Infos.dat | awk '{ print $2 }')
```

2. **Set script directory** for helper scripts: [file:306]

```bash
sd=/home/ibsbnicigt/nikolaev_d/script_library/4-solvents/
```

3. **Load GROMACS and MPI environment** (cluster-specific module + env vars): [file:306]

```bash
module add gromacs/2019.4/gcc
export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"
```

Adjust this block to your local cluster configuration.
4. **Compute the distance time series** with `gmx distance`: [file:306]

```bash
gmx distance \
    -s npt.tpr \
    -f npt.xtc \
    -select 'resid 1 and name C03 plus resname CL and name CL' \
    -oall dist_UNK_CL.xvg
```

This selection measures the distance between atom `C03` of residue ID 1 and each `CL` atom (`resname CL and name CL`), writing all resulting series to `dist_UNK_CL.xvg`. [file:306]
5. **Obtain the helper Python script**: [file:306]

```bash
cp $sd/find_time_closest_3000.py .
```

6. **Activate Python environment** (Conda): [file:306]

```bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base
```

Ensure that this environment has NumPy installed; `find_time_closest_3000.py` has no other dependencies. [file:307]
7. **Find representative time** using `find_time_closest_3000.py`: [file:306]

```bash
tavg=$(python3 find_time_closest_3000.py dist_UNK_CL.xvg 1 2.0)
```

    - Input: `dist_UNK_CL.xvg`
    - Column index: `1` (typical for `time value`)
    - Window: last `2.0` ns of the trajectory [file:306][file:307]

The script prints the time (in ps) whose distance is closest to the average distance over that 2 ns window; this is captured in `tavg`. [file:307]
8. **Extract the snapshot at `tavg`** with `gmx trjconv`: [file:306]

```bash
gmx trjconv -s npt.tpr -f npt.xtc -dump "$tavg" -o ${Project}_avg_dist.gro << EOF
0
EOF
```

    - `-dump "$tavg"`: extract the frame at `tavg` ps.
    - Interactive selection (`0`) typically chooses the first group (often “System”); adjust as needed for your setup. [file:306]

The resulting structure is saved as `${Project}_avg_dist.gro`, e.g. `0a_avg_dist.gro`. [file:306]

## 3. Typical workflow

1. After NPT equilibration, ensure you have `npt.tpr` and `npt.xtc` in the working directory.
2. Configure `extract_average.sh` if:
    - your solute reference atom is not `resid 1 and name C03`,
    - your ion name is not `CL`,
    - or you wish to use a different time window (change `2.0` passed to `find_time_closest_3000.py`). [file:306][file:307]
3. Run:

```bash
bash extract_average.sh
```

4. Use `<Project>_avg_dist.gro` as a representative structure for e.g. QM calculations, further equilibration, or visualization. [file:306]
```
<span style="display:none">[^1][^2]</span>

<div align="center">⁂</div>

[^1]: extract_average.sh
[^2]: find_time_closest_3000.py```
