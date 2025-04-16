import os
import pandas as pd

# Base directory for all runs
base_path = "/work/frejalom/ba"

# Define mapping of run folders to beta1
runs = {
    "experimental_run_50": 0.5,
    "experimental_run_51": 0.01,
    "experimental_run_55": 0.001,
    "experimental_run_54": 0.1,
    "experimental_run_53": "DA"
}

# Collect results separately by metric
time_stats = []
cluster_stats = []

# Loop through all runs
for run_name, beta1 in runs.items():
    run_path = os.path.join(base_path, run_name)

    for folder in os.listdir(run_path):
        folder_path = os.path.join(run_path, folder)
        if not os.path.isdir(folder_path) or not folder.startswith("o"):
            continue

        try:
            overlap = int(folder.split("bs")[0][1:])
            blocksize = int(folder.split("bs")[1])
        except:
            print(f"⚠️ Skipping malformed folder name: {folder}")
            continue

        csv_files = [f for f in os.listdir(folder_path) if f.endswith(".csv") and f.startswith('daten')]
        for file in csv_files:
            file_path = os.path.join(folder_path, file)
            try:
                df = pd.read_csv(file_path, delimiter=';', skiprows=2)
                df = df.iloc[:, :3]
                df.columns = ["Checkpoint", "Time (µs)", "Number of Clusters"]

                for metric, target_list in [("Time (µs)", time_stats), ("Number of Clusters", cluster_stats)]:
                    grouped = df.groupby("Checkpoint")[metric].agg(["mean", "median", "std"]).reset_index()
                    for _, row in grouped.iterrows():
                        target_list.append({
                            "Run": run_name,
                            "Beta1": beta1,
                            "Overlap": overlap,
                            "Blocksize": blocksize,
                            "Checkpoint": row["Checkpoint"],
                            "Mean": row["mean"],
                            "Median": row["median"],
                            "Std": row["std"]
                        })
            except Exception as e:
                print(f"❌ Error processing {file_path}: {e}")
        print(f"✅ Processed all files in {folder}")

# Create DataFrames
time_df = pd.DataFrame(time_stats)
cluster_df = pd.DataFrame(cluster_stats)

# Group to remove duplicates (average same config across files)
time_df = time_df.groupby(["Run", "Beta1", "Overlap", "Blocksize", "Checkpoint"]).agg({
    "Mean": "mean",
    "Median": "mean",
    "Std": "mean"
}).reset_index()

cluster_df = cluster_df.groupby(["Run", "Beta1", "Overlap", "Blocksize", "Checkpoint"]).agg({
    "Mean": "mean",
    "Median": "mean",
    "Std": "mean"
}).reset_index()

# Save to CSVs
output_dir = "/t3home/frejalom/cmssw/CMSSW_14_2_0_pre4/src/usercode/download"
time_csv_path = os.path.join(output_dir, "summary_stats_time.csv")
cluster_csv_path = os.path.join(output_dir, "summary_stats_clusters.csv")

time_df.to_csv(time_csv_path, index=False)
cluster_df.to_csv(cluster_csv_path, index=False)

print(f"✅ Saved cleaned time stats to: {time_csv_path}")
print(f"✅ Saved cleaned cluster stats to: {cluster_csv_path}")
