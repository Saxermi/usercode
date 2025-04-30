import pandas as pd
import matplotlib.pyplot as plt
import os

def _load_and_filter(csv_path, **filters):
    """
    Load the CSV, filter by Overlap, Blocksize, Beta1 as given in filters,
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

def plot_beta1_vs_time(csv_path, overlap_val, blocksize_val):
    df = _load_and_filter(csv_path, Overlap=overlap_val, Blocksize=blocksize_val)
    # split out the reference run for a black “DA” marker
    ref   = df[df.Run == 'experimental_run_66']
    other = df[df.Run != 'experimental_run_56']
    
    plt.figure()
    plt.scatter(other.Beta1, other.Time_s, marker='o', label='CascadeDA')
    plt.scatter(ref.Beta1,   ref.Time_s,   marker='o', color='black', label='DA')
    plt.xlabel('firstbestastop')
    plt.ylabel('Total Clustering Time (s)')
    plt.title(f'Influence of firstbestastop on Clustering Time\n'
              f'Overlap={overlap_val}, Blocksize={blocksize_val}')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    out = f'Influence_of_firstbestastop_o{overlap_val}_bs{blocksize_val}.png'
    plt.savefig(out)
    plt.close()
    print("✅", out)

def plot_overlap_vs_time(csv_path, beta1_val, blocksize_val):
    # 1) load all totals for this blocksize
    df = _load_and_filter(csv_path, Blocksize=blocksize_val)

    sub = df[df.Beta1.isin([beta1_val, 'DA'])].copy()

    # 3) convert µs→s
    sub['Time_s'] = sub['Mean'] * 1e-6

    # 4) split out
    numeric = sub[sub.Beta1 == beta1_val]
    da_rows = sub[sub.Beta1 == 'DA']

    if numeric.empty:
        raise ValueError(f"No runs with Beta1={beta1_val}")
    if da_rows.empty:
        raise ValueError("No DA reference rows found")

    # 5) aggregate DA down to one point
    da_x = da_rows['Overlap'].median()
    da_y = da_rows['Time_s'].mean()

    # 6) plot
    plt.figure()
    plt.scatter(numeric.Overlap, numeric.Time_s,
                marker='o', label='CascadeDA')
    plt.scatter(da_x, da_y,
                marker='o', color='black', label='DA')

    plt.xlabel('Overlap')
    plt.ylabel('Total Clustering Time (s)')
    plt.title(
        f'Influence of Overlap on Clustering Time\n'
        f'firstbestastop={beta1_val}, Blocksize={blocksize_val}'
    )
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    fname = f'overlap_firstbestastop{beta1_val}_bs{blocksize_val}.png'
    plt.savefig(fname)
    plt.close()
    print(f"✅ Saved {fname}")


def plot_blocksize_vs_time(csv_path, beta1_val, overlap_val):
    # 1) load all totals for this overlap
    df = _load_and_filter(csv_path, Overlap=overlap_val)

    sub = df[df.Beta1.isin([beta1_val, 'DA'])].copy()

    # 3) convert µs→s
    sub['Time_s'] = sub['Mean'] * 1e-6

    # 4) split out
    numeric = sub[sub.Beta1 == beta1_val]
    da_rows = sub[sub.Beta1 == 'DA']

    if numeric.empty:
        raise ValueError(f"No runs with Beta1={beta1_val}")
    if da_rows.empty:
        raise ValueError("No DA reference rows found")

    # 5) aggregate DA down to one point
    da_x = da_rows['Blocksize'].median()
    da_y = da_rows['Time_s'].mean()

    # 6) plot
    plt.figure()
    plt.scatter(numeric.Blocksize, numeric.Time_s,
                marker='o', label='CascadeDA')
    plt.scatter(da_x, da_y,
                marker='o', color='black', label='DA')

    plt.xlabel('Blocksize')
    plt.ylabel('Total Clustering Time (s)')
    plt.title(
        f'Influence of Blocksize on Clustering Time\n'
        f'firstbestastop={beta1_val}, Overlap={overlap_val}'
    )
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    fname = f'blocksize_sweep_b{beta1_val}_o{overlap_val}.png'
    plt.savefig(fname)
    plt.close()
    print(f"✅ Saved {fname}")


if __name__ == "__main__":
    csv_file = '/t3home/frejalom/cmssw/CMSSW_15_0_0_pre2/src/usercode/download/summary_stats_time_60.csv'
        
    # 1) Beta1 sweep
    plot_beta1_vs_time  (csv_file, overlap_val=0,  blocksize_val=256)
    
    # 2) Overlap sweep (fix Beta1 & blocksize)
    plot_overlap_vs_time(csv_file, beta1_val='0.5', blocksize_val=256)
    
    # 3) Blocksize sweep (fix Beta1 & overlap)
    plot_blocksize_vs_time(csv_file, beta1_val='0.5', overlap_val=0)
