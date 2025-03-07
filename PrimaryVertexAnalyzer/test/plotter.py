import os
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# Step 1: Find all downloaded CSV files
csv_files = [f for f in os.listdir() if f.endswith('.csv')]

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

# Concatenate all data and compute averages
all_data_df = pd.concat(all_data)
averaged_df = all_data_df.groupby("Checkpoint", as_index=False).mean()

# Ensure the order of checkpoints is maintained based on the first file read
averaged_df["Checkpoint"] = pd.Categorical(
    averaged_df["Checkpoint"], categories=order_of_checkpoints, ordered=True
)
averaged_df = averaged_df.sort_values(
    by="Checkpoint", key=lambda x: x.map({k: i for i, k in enumerate(order_of_checkpoints)})
)

# Get the earliest date from filenames
try:
    earliest_date = min([datetime.strptime(f[6:20], "%Y%m%d_%H%M%S") for f in csv_files])
except ValueError:
    print("❌ Error parsing timestamps from filenames. Ensure the format is correct.")
    exit(1)

# Number of CSV files used for averaging
n_files = len(csv_files)

# Plot 1: Time taken for each step (bar chart)
plt.figure(figsize=(10, 5))
plt.barh(averaged_df["Checkpoint"], averaged_df["Time (µs)"], color="blue")
plt.xlabel("Time (µs)")
plt.ylabel("Processing Steps")
plt.title(f"Time taken for each step (n={n_files}, Earliest: {earliest_date.strftime('%Y-%m-%d %H:%M:%S')})")
plt.gca().invert_yaxis()
plt.savefig("time_per_step.png")
print("✅ Saved: time_per_step.png")

# Plot 2: Size after each step (bar chart)
plt.figure(figsize=(10, 5))
plt.barh(averaged_df["Checkpoint"], averaged_df["Number of Clusters"], color="green")
plt.xlabel("Size after step")
plt.ylabel("Processing Steps")
plt.title(f"Size after each step (n={n_files}, Earliest: {earliest_date.strftime('%Y-%m-%d %H:%M:%S')})")
plt.gca().invert_yaxis()
plt.savefig("size_per_step.png")
print("✅ Saved: size_per_step.png")

# Plot 3: Overlaying both (time and size)
fig, ax1 = plt.subplots(figsize=(10, 5))

ax1.set_xlabel("Processing Steps")
ax1.set_ylabel("Time (µs)", color="blue")
ax1.bar(averaged_df["Checkpoint"], averaged_df["Time (µs)"], color="blue", alpha=0.6, label="Time (µs)")
ax1.tick_params(axis="y", labelcolor="blue")

ax2 = ax1.twinx()
ax2.set_ylabel("Size after step", color="green")
ax2.plot(averaged_df["Checkpoint"], averaged_df["Number of Clusters"], color="green", marker="o", linestyle="dashed", label="Size after step")
ax2.tick_params(axis="y", labelcolor="green")

plt.title(f"Time and Size after each step (n={n_files}, Earliest: {earliest_date.strftime('%Y-%m-%d %H:%M:%S')})")
fig.autofmt_xdate(rotation=30)
plt.savefig("time_and_size_overlay.png")
print("✅ Saved: time_and_size_overlay.png")

print("✅ All plots saved successfully!")
