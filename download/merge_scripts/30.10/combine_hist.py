import ROOT
import re
import os

# List of histograms to retrieve and plot
histogram_names = [
    "PUBlockBordersvsZdeltayprofile",
    "PUBlockBordersvsEfficencyprofile",
    "PUPurityVsZaxisprofile",
    "SEEfficiencyVsZaxisProfile",
    "PUEfficiencyVsZaxisProfile",
    "SEResolutionNormalizedBlockprofile",
    "SEBlockBordersvsPurityprofile",
    "SEBlockBordersvsEfficencyprofile",
    "PUBlockBordersvsFakeVertProfi",
    "SEPurityVsZaxisProfile",
]

# List of ROOT files with full path
root_files = [
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_HiggsGluonFusion_01_512_20241024_052146.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_HiggsGluonFusion_01_n512_20241025_181727.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_TTbar_01_512_20241024_063454.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_TTbar_01_n512_20241025_180817.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_ZMM_01_512_20241024_064931.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_ZMM_01_n512_20241025_184510.root",
]

# Colors assigned to each subset if grouping by dataset
subset_colors = {
    "HiggsGluonFusion": ROOT.kRed,
    "TTbar": ROOT.kBlue,
    "ZMM": ROOT.kGreen+2,
}

# Line styles to differentiate files within each subset if grouping by dataset
line_styles = [1, 2, 3, 4]  # Solid, Dashed, Dotted, Dash-Dotted

# Distinct colors for individual files if not grouping by dataset
distinct_colors = [
    ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kOrange+2,
    ROOT.kCyan+1, ROOT.kYellow+2, ROOT.kPink+1, ROOT.kAzure+1, ROOT.kTeal+2,
    ROOT.kViolet, ROOT.kSpring+2
]

# Dictionary to hold errors for each file
error_log = {}

# Set global ROOT styles
ROOT.gStyle.SetOptStat(0)  # Disable statistics box
ROOT.gStyle.SetTextFont(42)  # Set text font to match CMS style

# Create a canvas
canvas = ROOT.TCanvas("canvas", "", 800, 700)

# Regex pattern to match and extract parts of the filename
pattern = re.compile(r"pvSubset_([A-Za-z]+)_(\d+)_(n?)(\d+)_\d{8}_\d{6}\.root")

# Axis label dictionary for each histogram
axis_labels = {
    "PUBlockBordersvsZdeltayprofile": ("Distance to nearest block (mm)", "Delta Y (mm)"),
   # "PUBlockBordersvsEfficencyprofile": ("Distance to nearest block (mm)", "Efficiency (%)"),
   # "PUPurityVsZaxisprofile": ("Z-axis position (mm)", "Purity (%)"),
   # "SEEfficiencyVsZaxisProfile": ("Z-axis position (mm)", "Efficiency (%)"),
   # "PUEfficiencyVsZaxisProfile": ("Z-axis position (mm)", "Efficiency (%)"),
   # "SEResolutionNormalizedBlockprofile": ("Normalized Block", "Resolution"),
   # "SEBlockBordersvsPurityprofile": ("Distance to nearest block (mm)", "Purity (%)"),
   # "SEBlockBordersvsEfficencyprofile": ("Distance to nearest block (mm)", "Efficiency (%)"),
   # "PUBlockBordersvsFakeVertProfi": ("Distance to nearest block (mm)", "Fake Vertex Rate (%)"),
    #"SEPurityVsZaxisProfile": ("Z-axis position (mm)", "Purity (%)"),
}

# Flag to toggle grouping by dataset
group_by_dataset = False  # Set to False to use distinct colors for each dataset

# Loop through each histogram name
for hist_name in histogram_names:
    # Prepare a legend and list for each histogram type
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    histograms_with_legend = []

    # Track line style usage within each subset
    subset_line_style_index = {}

    # Loop through the files and draw histograms
    first_histogram = True
    for i, file_path in enumerate(root_files):
        # Initialize or retrieve the error list for each file
        error_log.setdefault(file_path, [])

        # Open ROOT file
        root_file = ROOT.TFile.Open(file_path)

        # Check if the file opened successfully
        if root_file and not root_file.IsZombie():
            print(f"Successfully opened: {file_path}")
        else:
            error_log[file_path].append("Failed to open file")
            continue  # Skip this file and move to the next

        # Extract filename and match against the regex
        filename = os.path.basename(file_path)
        match = pattern.search(filename)
        if match:
            # Extract relevant parts from the filename
            subset_type = match.group(1)  # E.g., "HiggsGluonFusion", "ZMM"
            subset_num = match.group(2)   # E.g., "01"
            is_negative = match.group(3)  # E.g., "n" if present, otherwise ""
            blocksize = match.group(4)    # E.g., "512"

            # Determine block size with appropriate sign
            blocksize = int(blocksize)
            if is_negative:
                blocksize = -blocksize

            # Construct legend entry
            legend_entry = f"{subset_type} {subset_num} blocksize {blocksize}"
            print(f"Generated legend entry: {legend_entry}")
        else:
            legend_entry = filename  # Use the filename as the legend entry
            print(f"Filename {filename} does not match pattern. Using filename as legend entry.")

        # Navigate to the 'offlinePrimaryVertices/efficiency' directory
        efficiency_dir = root_file.Get("offlinePrimaryVertices/efficiency")
        if not efficiency_dir:
            error_log[file_path].append(
                "'offlinePrimaryVertices/efficiency' directory not found"
            )
            root_file.Close()
            continue

        # Retrieve the histogram for the current name
        hist = efficiency_dir.Get(hist_name)
        if hist:
            hist.SetDirectory(0)

            # Check flag and set color/line style accordingly
            if group_by_dataset:
                # Assign color and line style by subset type
                color = subset_colors.get(subset_type, ROOT.kBlack)
                subset_line_style_index.setdefault(subset_type, 0)
                line_style = line_styles[subset_line_style_index[subset_type] % len(line_styles)]
                subset_line_style_index[subset_type] += 1
            else:
                # Use distinct colors for each file
                color = distinct_colors[i % len(distinct_colors)]
                line_style = 1  # Uniform line style for all files

            hist.SetMarkerColor(color)
            hist.SetLineStyle(line_style)
            hist.SetLineColor(color)

            # Apply axis labels from dictionary if available
            if hist_name in axis_labels:
                x_label, y_label = axis_labels[hist_name]
                hist.GetXaxis().SetTitle(x_label)
                hist.GetYaxis().SetTitle(y_label)

            hist.GetXaxis().SetTitleSize(0.04)
            hist.GetYaxis().SetTitleSize(0.04)

            # Add to list of histograms for legend
            histograms_with_legend.append((hist, legend_entry))

            # Draw histogram
            draw_option = "E" if first_histogram else "E SAME"
            hist.Draw(draw_option)
            first_histogram = False
        else:
            error_log[file_path].append(f"Histogram '{hist_name}' not found")
        root_file.Close()

    # Finalize the legend
    for hist, legend_entry in histograms_with_legend:
        legend.AddEntry(hist, legend_entry, "lep")
    legend.Draw()

    # Add CMS simulation labels
    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextSize(0.040)
    label.DrawLatex(0.100, 0.920, "#bf{CMS Simulation}")
    label.DrawLatex(0.550, 0.920, "#sqrt{s} = 14 TeV")

    # Update and save canvas
    canvas.Update()
    output_file = f"/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/git reps/usercode/download/merge_scripts/30.10/imgout/combined_{hist_name}_overlap_histograms.png"
    canvas.SaveAs(output_file)

# Print formatted error report at the end
print("\nError Summary:")
for file_path, errors in error_log.items():
    error_text = "None" if not errors else "\n  - " + "\n  - ".join(errors)
    print(f"{file_path}:\n  Errors: {error_text}")

