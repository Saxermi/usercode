import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Paths
base_path = "/t3home/frejalom/cmssw/CMSSW_15_0_0_pre2/src/usercode/download"
cluster_csv = os.path.join(base_path, "summary_stats_clusters_50-54.csv")
time_csv = os.path.join(base_path, "summary_stats_time_50-54.csv")
print('got paths')
# Load data
cluster_df = pd.read_csv(cluster_csv)
time_df = pd.read_csv(time_csv)
print('loaded data')

cluster_df["Beta1"] = cluster_df["Beta1"].astype(str)
time_df["Beta1"] = time_df["Beta1"].astype(str)

# Style setup
sns.set(style="whitegrid")
beta1_styles = {0.001: "o", 0.01: "s", 0.5: "D"}  # Marker by beta1
palette = sns.color_palette("viridis", 6)  # Color by block size

def scatter_plot_per_checkpoint(df, metric_name, ylabel, filename_prefix):
    import matplotlib.pyplot as plt
    import seaborn as sns

    checkpoints = df["Checkpoint"].unique()
    #checkpoints = df["Checkpoint"].unique()[:3]  # only plot first 3
    block_sizes = sorted(df["Blocksize"].unique())
    block_color_map = dict(zip(block_sizes, palette))

    for i, checkpoint in enumerate(checkpoints):
        print(f"🔄 Plotting checkpoint {i+1}/{len(checkpoints)}: {checkpoint}")

        # Subset once
        subset = df[df["Checkpoint"] == checkpoint][["Overlap", "Mean", "Blocksize", "Beta1"]].copy()
        if subset.empty:
            print(f"⚠️ Skipping empty checkpoint: {checkpoint}")
            continue

        # Map visual features
        subset["Color"] = subset.apply(
            lambda row: "black" if row["Beta1"] == "DA" else block_color_map.get(row["Blocksize"], "gray"), axis=1
        )
        subset["Marker"] = subset["Beta1"].map(beta1_styles).fillna("*")

        # Set up plot
        fig, ax = plt.subplots(figsize=(8, 6))

        # Plot each (beta1) marker separately for color control
        plotted_labels = set()
        for marker in subset["Marker"].unique():
            marker_subset = subset[subset["Marker"] == marker]
            for bs in marker_subset["Blocksize"].unique():
                bs_subset = marker_subset[marker_subset["Blocksize"] == bs]
                label = f"DA Reference" if bs_subset["Beta1"].iloc[0] == "DA" else f"BS{bs}, β={bs_subset['Beta1'].iloc[0]}"
                if label in plotted_labels:
                    label = None
                else:
                    plotted_labels.add(label)

                ax.scatter(
                    bs_subset["Overlap"].values,
                    bs_subset["Mean"].values,
                    color=block_color_map[bs],
                    marker=marker,
                    s=50,
                    label=label
                )

        ax.set_title(f"{metric_name} — {checkpoint}")
        ax.set_xlabel("Overlap")
        ax.set_ylabel(ylabel)
        ax.legend(
            loc="center left",
            bbox_to_anchor=(1.05, 0.5),
            fontsize="small",
            borderaxespad=0.
        )
        plt.tight_layout()

        safe_checkpoint = checkpoint.replace(" ", "_").replace(";", "_").replace("/", "_")
        out_path = os.path.join("/t3home/frejalom/cmssw/CMSSW_14_2_0_pre4/src/usercode/download/plotter2Results", f"{filename_prefix}_{safe_checkpoint}.png")
        plt.savefig(out_path)
        print(f"✅ Saved: {out_path}")
        plt.close()

# Cluster size scatter plots per checkpoint
scatter_plot_per_checkpoint(cluster_df, "Number of Clusters", "Cluster Size", "clusters")

# Time scatter plots per checkpoint
scatter_plot_per_checkpoint(time_df, "Time (µs)", "Time (µs)", "time")
