import pandas as pd
import matplotlib.pyplot as plt
import os

def _load_and_filter(csv_path, **filters):
    """
    Load the CSV, filter by Overlap, Blocksize, T_stop as given in filters,
    and only keep the 'final time and size' checkpoint.
    Returns a DataFrame with an extra column 'Time_s'.
    """
    df = pd.read_csv(csv_path)
    # always filter checkpoint
    mask = (df.Checkpoint == 'final time and size')
    # apply any additional filters passed in
    for col, val in filters.items():
        mask &= (df[col] == val)
    out = df[mask].copy()
    if out.empty:
        raise ValueError(f"No rows for {filters}")
    # convert µs → seconds
    out['Time_s'] = out['Mean'] * 1e-6
    return out

def plot_T_stop_vs_time(csv_path, overlap_val, blocksize_val):
    df = _load_and_filter(csv_path, Overlap=overlap_val, Blocksize=blocksize_val)

    # Convert T_stop to numeric where possible
    df['T_stop_numeric'] = pd.to_numeric(df['T_stop'], errors='coerce')

    # Separate references
    da_row  = df[df.Run == 'experimental_run_95']
    dab_row = df[df.Run == 'experimental_run_72']
    cascade = df[~df.Run.isin(['experimental_run_95', 'experimental_run_72'])]

    plt.figure()

    # CascadeDA as scatter
    plt.scatter(cascade['T_stop_numeric'], cascade['Time_s'], marker='o', label='CascadeDA')

    # DA and DAB as horizontal lines
    if not da_row.empty:
        da_val = da_row['Time_s'].mean()
        plt.axhline(da_val, color='black', linestyle='--', label='DA')

    if not dab_row.empty:
        dab_val = dab_row['Time_s'].mean()
        plt.axhline(dab_val, color='gray', linestyle='--', label='DAB')

    plt.xlabel('T_stop')
    plt.ylabel('Total Clustering Time (s)')
    plt.title(f'Influence of T_stop on Clustering Time\n'
              f'Overlap={overlap_val}, Blocksize={blocksize_val}')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    out = f'Influence_of_firstbestastop_o{overlap_val}_bs{blocksize_val}.png'
    plt.savefig(out)
    plt.close()
    print("✅", out)


def plot_overlap_vs_time(csv_path, T_stop_val, blocksize_val):
    df = _load_and_filter(csv_path, Blocksize=blocksize_val)

    # Use only rows relevant for this plot
    sub = df[df.T_stop.isin([T_stop_val, 'DA', 'DAB'])].copy()
    sub['Time_s'] = sub['Mean'] * 1e-6

    # Split by Run
    numeric = sub[sub.T_stop == T_stop_val]
    da_row = sub[sub.Run == 'experimental_run_95']
    dab_row = sub[sub.Run == 'experimental_run_72']

    if numeric.empty:
        raise ValueError(f"No runs with T_stop={T_stop_val}")
    
    plt.figure()

    # CascadeDA points
    plt.scatter(numeric.Overlap, numeric.Time_s, marker='o', label='CascadeDA')

    # Reference lines
    if not da_row.empty:
        da_val = da_row['Time_s'].mean()
        plt.axhline(da_val, color='black', linestyle='--', label='DA')

    if not dab_row.empty:
        dab_val = dab_row['Time_s'].mean()
        plt.axhline(dab_val, color='gray', linestyle='--', label='DAB')

    plt.xlabel('Overlap')
    plt.ylabel('Total Clustering Time (s)')
    plt.title(
        f'Influence of Overlap on Clustering Time\n'
        f'firstbestastop={T_stop_val}, Blocksize={blocksize_val}'
    )
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    fname = f'overlap_firstbestastop{T_stop_val}_bs{blocksize_val}.png'
    plt.savefig(fname)
    plt.close()
    print(f"✅ Saved {fname}")


def plot_blocksize_vs_time(csv_path, T_stop_val, overlap_val):
    df = _load_and_filter(csv_path, Overlap=overlap_val)

    sub = df[df.T_stop.isin([T_stop_val, 'DA', 'DAB'])].copy()
    sub['Time_s'] = sub['Mean'] * 1e-6

    numeric = sub[sub.T_stop == T_stop_val]
    da_row = sub[sub.Run == 'experimental_run_95']
    dab_row = sub[sub.Run == 'experimental_run_72']

    if numeric.empty:
        raise ValueError(f"No runs with T_stop={T_stop_val}")
    
    plt.figure()

    # CascadeDA points
    plt.scatter(numeric.Blocksize, numeric.Time_s, marker='o', label='CascadeDA')

    # Reference lines
    if not da_row.empty:
        da_val = da_row['Time_s'].mean()
        plt.axhline(da_val, color='black', linestyle='--', label='DA')

    if not dab_row.empty:
        dab_val = dab_row['Time_s'].mean()
        plt.axhline(dab_val, color='gray', linestyle='--', label='DAB')

    plt.xlabel('Blocksize')
    plt.ylabel('Total Clustering Time (s)')
    plt.title(
        f'Influence of Blocksize on Clustering Time\n'
        f'firstbestastop={T_stop_val}, Overlap={overlap_val}'
    )
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    fname = f'blocksize_sweep_b{T_stop_val}_o{overlap_val}.png'
    plt.savefig(fname)
    plt.close()
    print(f"✅ Saved {fname}")



if __name__ == "__main__":
    csv_file = '/t3home/frejalom/cmssw/CMSSW_15_0_0_pre2/src/usercode/download/summary_stats_time_60.csv'
        
    # 1) T_stop sweep
    plot_T_stop_vs_time  (csv_file, overlap_val=0,  blocksize_val=256)
    
    # 2) Overlap sweep (fix T_stop & blocksize)
    plot_overlap_vs_time(csv_file, T_stop_val='0.5', blocksize_val=256)
    
    # 3) Blocksize sweep (fix T_stop & overlap)
    plot_blocksize_vs_time(csv_file, T_stop_val='0.5', overlap_val=0)
