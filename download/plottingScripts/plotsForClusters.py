import pandas as pd
import matplotlib.pyplot as plt
import os

def _load_and_filter_clusters(csv_path, **filters):

    df = pd.read_csv(csv_path)
    mask = (df.Checkpoint.str.strip() == 'final time and size')
    for col, val in filters.items():
        mask &= (df[col] == val)
    out = df[mask].copy()
    if out.empty:
        raise ValueError(f"No rows for {filters}")
    return out


def plot_beta1_vs_clusters(csv_path, overlap_val, blocksize_val):
    df = _load_and_filter_clusters(csv_path,
                                   Overlap=overlap_val,
                                   Blocksize=blocksize_val)

    ref   = df[df.Run == 'experimental_run_56']
    other = df[df.Run != 'experimental_run_56']

    plt.figure()
    plt.scatter(other.Beta1, other.Mean,
                marker='o', label='CascadeDA')
    plt.scatter(ref.Beta1,   ref.Mean,
                marker='o', color='black', label='DA')

    plt.xlabel('firstbestastop')
    plt.ylabel('Mean # clusters')
    plt.title(
        f'Influence of firstbestastop on #Clusters\n'
        f'Overlap={overlap_val}, Blocksize={blocksize_val}'
    )
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    out = f'clusters_firstbestastop_o{overlap_val}_bs{blocksize_val}.png'
    plt.savefig(out)
    plt.close()
    print("✅", out)


def plot_overlap_vs_clusters(csv_path, beta1_val, blocksize_val):
    df = _load_and_filter_clusters(csv_path, Blocksize=blocksize_val)

    sub = df[df.Beta1.isin([str(beta1_val), 'DA'])].copy()
    numeric = sub[sub.Beta1 == str(beta1_val)]
    da_rows = sub[sub.Beta1 == 'DA']

    # aggregate DA to one point
    da_x = da_rows['Overlap'].median()
    da_y = da_rows['Mean'].mean()

    plt.figure()
    plt.scatter(numeric.Overlap, numeric.Mean,
                marker='o', label='CascadeDA')
    plt.scatter(da_x, da_y,
                marker='o', color='black', label='DA')

    plt.xlabel('Overlap')
    plt.ylabel('Mean # clusters')
    plt.title(
        f'Influence of Overlap on #Clusters\n'
        f'firstbestastop={beta1_val}, Blocksize={blocksize_val}'
    )
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    out = f'clusters_overlap_firstbestastop{beta1_val}_bs{blocksize_val}.png'
    plt.savefig(out)
    plt.close()
    print("✅", out)


def plot_blocksize_vs_clusters(csv_path, beta1_val, overlap_val):
    df = _load_and_filter_clusters(csv_path, Overlap=overlap_val)

    sub = df[df.Beta1.isin([str(beta1_val), 'DA'])].copy()
    numeric = sub[sub.Beta1 == str(beta1_val)]
    da_rows = sub[sub.Beta1 == 'DA']

    da_x = da_rows['Blocksize'].median()
    da_y = da_rows['Mean'].mean()

    plt.figure()
    plt.scatter(numeric.Blocksize, numeric.Mean,
                marker='o', label='CascadeDA')
    plt.scatter(da_x, da_y,
                marker='o', color='black', label='DA')

    plt.xlabel('Blocksize')
    plt.ylabel('Mean # clusters')
    plt.title(
        f'Influence of Blocksize on #Clusters\n'
        f'firstbestastop={beta1_val}, Overlap={overlap_val}'
    )
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    out = f'clusters_blocksize_firstbestastop{beta1_val}_o{overlap_val}.png'
    plt.savefig(out)
    plt.close()
    print("✅", out)


if __name__ == "__main__":
    csv_clusters = '/t3home/frejalom/cmssw/CMSSW_15_0_0_pre2/src/usercode/download/summary_stats_clusters_50-54.csv'

    # 1) Beta1 sweep
    plot_beta1_vs_clusters(csv_clusters, overlap_val=0,  blocksize_val=256)

    # 2) Overlap sweep
    plot_overlap_vs_clusters(csv_clusters, beta1_val='0.5', blocksize_val=256)

    # 3) Blocksize sweep
    plot_blocksize_vs_clusters(csv_clusters, beta1_val='0.5', overlap_val=0)
