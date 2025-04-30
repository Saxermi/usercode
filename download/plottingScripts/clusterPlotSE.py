import pandas as pd
import matplotlib.pyplot as plt

# Load and filter the data
df = pd.read_csv('/t3home/frejalom/cmssw/CMSSW_15_0_0_pre2/src/usercode/download/summary_stats_clusters_60.csv')
df_final = df[
    (df['Checkpoint'] == 'final time and size') & 
    (df['Overlap'] == 0) & 
    (df['Blocksize'] == 256)
].copy()

# Define groups for special runs
dab_runs = ['experimental_run_63', 'experimental_run_64']
df_final['Group'] = df_final['Run'].apply(
    lambda x: 'DAB' if x in dab_runs else ('DA' if x == 'experimental_run_66' else 'CascadeDA')
)

# Order Beta1 categories
categories = ['0.01', '0.1', '0.25', '0.5', 'DA']
df_final['Beta1'] = pd.Categorical(df_final['Beta1'], categories=categories, ordered=True)

# Set up color mapping for event types and marker mapping for groups
event_types = df_final['EventType'].unique()
cmap = plt.get_cmap('tab10')
color_map = {etype: cmap(i) for i, etype in enumerate(event_types)}
marker_map = {'CascadeDA': 'o', 'DA': 's', 'DAB': '^'}

# Create the plot
fig, ax = plt.subplots(figsize=(10, 6))
for _, row in df_final.iterrows():
    ax.plot(
        row['Beta1'], row['Mean'],
        marker=marker_map[row['Group']],
        color=color_map[row['EventType']],
        markersize=10,
        linestyle=''    # no connecting line
    )

# Build legends
# Legend for event types
for etype in event_types:
    ax.plot([], [], marker='o', linestyle='', color=color_map[etype], label=etype)
# Legend for groups
for grp, marker in marker_map.items():
    ax.plot([], [], marker=marker, linestyle='', color='black', label=grp)

# Labels, title, and styling
ax.set_xlabel('firstbestastop', fontsize=12)
ax.set_ylabel('Mean Number of Clusters', fontsize=12)
ax.set_title('Inlfuence of firstbestastop on Number of Clusters', fontsize=14)
ax.set_xticks(categories)
ax.legend(title='SE Types & Algorithm', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
fname = f'clusterVsFirstbetastopNoErrors.png'
plt.savefig(fname)
plt.close()
print(f"✅ Saved {fname}")
