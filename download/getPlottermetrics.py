import os
import re
import pandas as pd

# Base directory for all runs
base_path = "/work/frejalom/ba"

# Define mapping of run folders to beta1
runs = {
    "experimental_run_62": "DAB"
}

# Collect results separately by metric
time_stats = []
cluster_stats = []

# Loop through all runs
for run_name, beta1 in runs.items():
    run_path = os.path.join(base_path, run_name)

    # Step 1: loop over event-type directories
    for event_type in os.listdir(run_path):
        event_path = os.path.join(run_path, event_type)
        if not os.path.isdir(event_path):
            continue

        # Step 2: within each event type, loop over overlap/blocksize folders
        for folder in os.listdir(event_path):
            folder_path = os.path.join(event_path, folder)
            if not os.path.isdir(folder_path):
                continue

            # parse overlap & blocksize from folder name
            m = re.match(r"^o(\d+)bs(\d+)$", folder)
            if m:
                overlap   = int(m.group(1))
                blocksize = int(m.group(2))
            else:
                m2 = re.match(r"^(\d+)_(\d+)$", folder)
                if m2:
                    overlap   = int(m2.group(1))
                    blocksize = int(m2.group(2))
                else:
                    print(f"⚠️ Skipping non-config folder: {folder}")
                    continue

            # find all the 'daten*.csv' files
            csv_files = [
                f for f in os.listdir(folder_path)
                if f.endswith(".csv") and f.startswith("daten")
            ]

            for file in csv_files:
                file_path = os.path.join(folder_path, file)
                try:
                    # load and trim to the 3 relevant columns
                    df = pd.read_csv(file_path, delimiter=";", skiprows=2)
                    df = df.iloc[:, :3]
                    df.columns = ["Checkpoint", "Time (µs)", "Number of Clusters"]

                    # for each metric, compute per-file mean/median/std
                    for metric, target_list in [
                        ("Time (µs)", time_stats),
                        ("Number of Clusters", cluster_stats)
                    ]:
                        grouped = (
                            df
                            .groupby("Checkpoint")[metric]
                            .agg(["mean", "median", "std"])
                            .reset_index()
                        )
                        for _, row in grouped.iterrows():
                            target_list.append({
                                "Run":       run_name,
                                "Beta1":     beta1,
                                "EventType": event_type,
                                "Overlap":   overlap,
                                "Blocksize": blocksize,
                                "Checkpoint": row["Checkpoint"],
                                "Mean":      row["mean"],
                                "Median":    row["median"],
                                "Std":       row["std"],
                            })
                except Exception as e:
                    print(f"❌ Error processing {file_path}: {e}")

            print(f"✅ Processed all files in {event_type}/{folder}")

# Build DataFrames of all the per‐file stats
time_df    = pd.DataFrame(time_stats)
cluster_df = pd.DataFrame(cluster_stats)

# Function to compute true mean/median/std across files per unique configuration
def summarize(df, value_col="Mean"):
    summary = (
        df
        .groupby(
            ["Run", "Beta1", "EventType", "Overlap", "Blocksize", "Checkpoint"]
        )[value_col]
        .agg(["mean", "median", "std"])
        .reset_index()
    ).rename(columns={
        "mean":   "Mean",
        "median": "Median",
        "std":    "Std",
    })
    return summary

time_summary    = summarize(time_df,    value_col="Mean")
cluster_summary = summarize(cluster_df, value_col="Mean")

# Save to CSV
output_dir       = "/t3home/frejalom/cmssw/CMSSW_15_0_0_pre2/src/usercode/download"
time_csv_path    = os.path.join(output_dir, "summary_stats_time_old_dab.csv")
cluster_csv_path = os.path.join(output_dir, "summary_stats_clusters_old_dab.csv")

time_summary.to_csv(time_csv_path,    index=False)
cluster_summary.to_csv(cluster_csv_path, index=False)

print(f"✅ Saved time summary (with Std) to:    {time_csv_path}")
print(f"✅ Saved cluster summary (with Std) to: {cluster_csv_path}")
