import os
import math
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from datetime import datetime

# Step 1: Find all downloaded CSV files

def find_csv_files(exp_run_dir, overlap, blocksize):
    # Build the pattern with a wildcard for the subset directory
    pattern = os.path.join(exp_run_dir, '*', str(overlap), str(blocksize), '*.csv')
    csv_files = glob.glob(pattern)
    return csv_files

# Example usage:
exp_run_dir = '/work/msaxer/ba/experimental_run_29'
overlap = 50
blocksize = 512

csv_files = find_csv_files(exp_run_dir, overlap, blocksize)

if not csv_files:
    print("❌ No CSV files found. Make sure you've downloaded the data before running this script.")
    exit(1)

# Initialize lists to collect data
all_data = []
order_of_checkpoints = []

# Read all CSV files and extract relevant data
for file in csv_files:
    try:
        df = pd.read_csv(file, delimiter=';', skiprows=2)  # Skip first two lines (descriptive text)
        df = df.iloc[:, :3]  # Keep only relevant columns: Checkpoint, Time it took, Number of clusters
        df.columns = ["Checkpoint", "Time (µs)", "Number of Clusters"]  # Standardize column names
        all_data.append(df)

        # Record order of checkpoints if this is the first file being read
        if not order_of_checkpoints:
            order_of_checkpoints = df["Checkpoint"].tolist()

    except Exception as e:
        print(f"⚠️ Warning: Could not process file {file}. Error: {e}")

# Ensure we have valid data before proceeding
if not all_data:
    print("❌ No valid CSV files were processed. Please check the data format.")
    exit(1)

# Concatenate all data from current run (ensuring no Datenübernahme aus vorherigen Runden)
all_data_df = pd.concat(all_data)

# Ensure the order of checkpoints is maintained based on the first file read
all_data_df["Checkpoint"] = pd.Categorical(
    all_data_df["Checkpoint"], categories=order_of_checkpoints, ordered=True
)
all_data_df = all_data_df.sort_values("Checkpoint")

# Get the earliest date from filenames
try:
    earliest_date = min([datetime.strptime(os.path.basename(f)[6:20], "%Y%m%d_%H%M%S") for f in csv_files])
except ValueError:
    print("❌ Error parsing timestamps from filenames. Ensure the format is correct.")
    exit(1)

# Number of CSV files used for averaging
n_files = len(csv_files)

# Prepare data for boxplots per checkpoint (freshly computed for this run)
checkpoints = order_of_checkpoints
time_data = [all_data_df.loc[all_data_df["Checkpoint"] == cp, "Time (µs)"].values for cp in checkpoints]
cluster_data = [all_data_df.loc[all_data_df["Checkpoint"] == cp, "Number of Clusters"].values for cp in checkpoints]

# Determine figure width based on longest checkpoint label
max_label_length = max(len(str(cp)) for cp in checkpoints)
fig_width = max(10, len(checkpoints) * 1.5 + max_label_length * 0.1)

# Dynamically adjust scale for "Number of Clusters" in steps of 25
max_cluster = all_data_df["Number of Clusters"].max()
max_cluster_tick = math.ceil(max_cluster / 25) * 25

# Define moderate colors
time_color = "steelblue"
size_color = "seagreen"

# Get current creation time text for footer
created_text = f"Plot created at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"

# --------------------------
# Plot 1: Boxplot for Time (µs)
# --------------------------
plt.figure(figsize=(fig_width, 6))
bp_time = plt.boxplot(time_data, positions=range(1, len(checkpoints) + 1), patch_artist=True)
# Set box colors for time
for box in bp_time['boxes']:
    box.set(facecolor=time_color)
plt.xticks(range(1, len(checkpoints) + 1), checkpoints, rotation=45, ha='right')
plt.ylabel("Time (µs)")
plt.title(f"Time taken for each step (n={n_files}, Earliest: {earliest_date.strftime('%Y-%m-%d %H:%M:%S')})")

# Compute statistics freshly for this plot
time_stats = []
for cp, data in zip(checkpoints, time_data):
    if len(data) > 0:
        mean_val = data.mean()
        median_val = pd.Series(data).median()
        std_val = data.std()
        time_stats.append(f"{cp}: mean={mean_val:.1f}, median={median_val:.1f}, std={std_val:.1f}")
stats_text = " | ".join(time_stats) + "    " + created_text

plt.figtext(0.5, 0.01, stats_text, wrap=True, horizontalalignment='center', fontsize=8)
plt.tight_layout(rect=[0, 0.03, 1, 1])
plt.savefig(f"time_per_step_overlap_{overlap}_blocksize{blocksize}.png")
print(f"✅ Saved: time_per_step_overlap_{overlap}_blocksize{blocksize}.png")
plt.close()

# --------------------------
# Plot 2: Boxplot for Number of Clusters
# --------------------------
plt.figure(figsize=(fig_width, 6))
bp_cluster = plt.boxplot(cluster_data, positions=range(1, len(checkpoints) + 1), patch_artist=True)
# Set box colors for clusters
for box in bp_cluster['boxes']:
    box.set(facecolor=size_color)
plt.xticks(range(1, len(checkpoints) + 1), checkpoints, rotation=45, ha='right')
plt.ylabel("Size after step")
plt.title(f"Size after each step (n={n_files}, Earliest: {earliest_date.strftime('%Y-%m-%d %H:%M:%S')})")

# Set y-axis for clusters in steps of 25
plt.ylim(0, max_cluster_tick)
plt.yticks(np.arange(0, max_cluster_tick + 1, 25))

# Compute statistics freshly for this plot
cluster_stats = []
for cp, data in zip(checkpoints, cluster_data):
    if len(data) > 0:
        mean_val = data.mean()
        median_val = pd.Series(data).median()
        std_val = data.std()
        cluster_stats.append(f"{cp}: mean={mean_val:.1f}, median={median_val:.1f}, std={std_val:.1f}")
stats_text = " | ".join(cluster_stats) + "    " + created_text

plt.figtext(0.5, 0.01, stats_text, wrap=True, horizontalalignment='center', fontsize=8)
plt.tight_layout(rect=[0, 0.03, 1, 1])
plt.savefig(f"size_per_step_overlap_{overlap}_blocksize{blocksize}.png")
print(f"✅ Saved: size_per_step_overlap_{overlap}_blocksize{blocksize}.png")
plt.close()

# --------------------------
# Plot 3: Combined overlay of Time (µs) and Number of Clusters using boxplots
# --------------------------
fig, ax1 = plt.subplots(figsize=(fig_width, 6))

# Offset positions so the two boxplots don't overlap completely
positions_time = [i - 0.15 for i in range(1, len(checkpoints) + 1)]
positions_cluster = [i + 0.15 for i in range(1, len(checkpoints) + 1)]

# Plot boxplot for Time on ax1
bp_time_overlay = ax1.boxplot(time_data, positions=positions_time, widths=0.3, patch_artist=True)
for box in bp_time_overlay['boxes']:
    box.set(facecolor=time_color)
ax1.set_ylabel("Time (µs)", color=time_color)
ax1.tick_params(axis="y", labelcolor=time_color)

# Plot boxplot for Number of Clusters on twin axis
ax2 = ax1.twinx()
bp_cluster_overlay = ax2.boxplot(cluster_data, positions=positions_cluster, widths=0.3, patch_artist=True)
for box in bp_cluster_overlay['boxes']:
    box.set(facecolor=size_color)
ax2.set_ylabel("Size after step", color=size_color)
ax2.tick_params(axis="y", labelcolor=size_color)

# Set dynamic scale for clusters on ax2
ax2.set_ylim(0, max_cluster_tick)
ax2.set_yticks(np.arange(0, max_cluster_tick + 1, 25))

ax1.set_xticks(range(1, len(checkpoints) + 1))
ax1.set_xticklabels(checkpoints, rotation=45, ha='right')
plt.title(f"Time and Size after each step (n={n_files}, Earliest: {earliest_date.strftime('%Y-%m-%d %H:%M:%S')})")

# Compute statistics freshly for the overlay plot
time_stats_overlay = []
for cp, data in zip(checkpoints, time_data):
    if len(data) > 0:
        mean_val = data.mean()
        median_val = pd.Series(data).median()
        std_val = data.std()
        time_stats_overlay.append(f"{cp}: mean={mean_val:.1f}, median={median_val:.1f}, std={std_val:.1f}")

cluster_stats_overlay = []
for cp, data in zip(checkpoints, cluster_data):
    if len(data) > 0:
        mean_val = data.mean()
        median_val = pd.Series(data).median()
        std_val = data.std()
        cluster_stats_overlay.append(f"{cp}: mean={mean_val:.1f}, median={median_val:.1f}, std={std_val:.1f}")

stats_text_overlay = ("Time (µs): " + " | ".join(time_stats_overlay) +
                      "\nSize: " + " | ".join(cluster_stats_overlay) + "    " + created_text)
plt.figtext(0.5, 0.01, stats_text_overlay, wrap=True, horizontalalignment='center', fontsize=8)

# Add legend using custom patches
legend_handles = [Patch(facecolor=time_color, label="Time (µs)"),
                  Patch(facecolor=size_color, label="Size after step")]
ax1.legend(handles=legend_handles, loc="upper right")

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig(f"time_and_size_overlay_overlap_{overlap}_blocksize{blocksize}.png")
print(f"✅ Saved: time_and_size_overlay_overlap_{overlap}_blocksize{blocksize}.png")
plt.close()

print("✅ All plots saved successfully!")
