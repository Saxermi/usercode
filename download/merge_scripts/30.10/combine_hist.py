import ROOT
import re
import os

# List of histograms to retrieve and plot
histogram_names = [
    # "PUBlockBordersvsZdeltayprofile",
    "PUBlockBordersvsEfficencyprofile",
    #   "PUPurityVsZaxisprofile",
    #   "SEEfficiencyVsZaxisProfile",
    "PUEfficiencyVsZaxisProfile",
    #  "SEResolutionNormalizedBlockprofile",
    #  "SEBlockBordersvsPurityprofile",
    #  "SEBlockBordersvsEfficencyprofile",
    #  "PUBlockBordersvsFakeVertProfi",
    # "SEPurityVsZaxisProfile",
]

# List of ROOT files with full path
root_files = [
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_TTbar_01_512_20241024_063454.root",
    # "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_TTbar_01_512_20241024_063454A.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_TTbar_01_n512_20241025_180817.root",
    # "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_HiggsGluonFusion_01_n512_20241025_181727.root",
    # Add other ROOT files as needed
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_HiggsGluonFusion_01_512_20241024_052146.root",
    # "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_ZMM_01_512_20241024_064931.root",
]

# Colors assigned to each subset if grouping by dataset
subset_colors = {
    "HiggsGluonFusion": ROOT.kRed,
    "TTbar": ROOT.kBlue,
    "ZMM": ROOT.kGreen + 2,
}

# Line styles to differentiate files within each subset if grouping by dataset
line_styles = [1, 2, 3, 4]  # Solid, Dashed, Dotted, Dash-Dotted

# Distinct colors for individual files if not grouping by dataset
distinct_colors = [
    ROOT.kRed,
    ROOT.kBlue,
    ROOT.kGreen + 2,
    ROOT.kMagenta,
    ROOT.kOrange + 2,
    ROOT.kCyan + 1,
    ROOT.kYellow + 2,
    ROOT.kPink + 1,
    ROOT.kAzure + 1,
    ROOT.kTeal + 2,
    ROOT.kViolet,
    ROOT.kSpring + 2,
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
    "PUBlockBordersvsZdeltayprofile": (
        "Distance to nearest block (mm)",
        "Delta Y (mm)",
    ),
    # Add other axis labels as needed
}

# Flag to toggle grouping by dataset
group_by_dataset = False  # Set to False to use distinct colors for each dataset


# Function to find histogram recursively
def find_histogram_in_directory(directory, hist_name):

    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        if (
            obj.InheritsFrom("TH1")
            or obj.InheritsFrom("TProfile")
            or obj.InheritsFrom("TH2")
        ):
            if obj.GetName() == hist_name:
                return obj
        elif obj.InheritsFrom("TDirectory"):
            found_hist = find_histogram_in_directory(obj, hist_name)
            if found_hist:
                return found_hist
    return None


histograms_with_legend = []

# Loop through each histogram name
for hist_name in histogram_names:
    # Prepare a legend and list for each histogram type
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

    # Track line style usage within each subset
    subset_line_style_index = {}

    # Clear the canvas before drawing new histograms
    canvas.Clear()

    # Loop through the files and draw histograms
    first_histogram = True
    root_file = []
    for i, file_path in enumerate(root_files):
        # Initialize or retrieve the error list for each file
        error_log.setdefault(file_path, [])

        # Open ROOT file
        root_file = ROOT.TFile.Open(file_path)

        root_file.cd("testVertices_test/efficiency")

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
            subset_type = match.group(1)
            subset_num = match.group(2)
            is_negative = match.group(3)
            blocksize = match.group(4)
            blocksize = int(blocksize)
            if is_negative:
                blocksize = -blocksize
            legend_entry = f"{subset_type} {subset_num} blocksize {blocksize}"
            print(f"Generated legend entry: {legend_entry}")
        else:
            legend_entry = filename
            print(
                f"Filename {filename} does not match pattern. Using filename as legend entry."
            )

        # Retrieve the histogram for the current name
        # hist = find_histogram_in_directory(root_file, hist_name)
        # Retrieve the histogram for the current name
        hist = find_histogram_in_directory(
            root_file.GetDirectory("testVertices_test/efficiency"), hist_name
        )
        print(f"Retrieved histogram for {hist_name}: {hist}")

        if hist:
            hist.SetDirectory(0)

            # Clone the histogram with a unique name to prevent overwriting
            unique_hist_name = f"{hist.GetName()}_{os.path.basename(file_path)}"
            cloned_hist = hist.Clone(unique_hist_name)
            cloned_hist.SetDirectory(0)
            root_file.Close()

            # Check if cloning was successful
            if not cloned_hist:
                error_log[file_path].append(f"Failed to clone histogram '{hist_name}'")
                root_file.Close()
                continue

            # Check flag and set color/line style accordingly
            if group_by_dataset:
                color = subset_colors.get(subset_type, ROOT.kBlack)
                subset_line_style_index.setdefault(subset_type, 0)
                line_style = line_styles[
                    subset_line_style_index[subset_type] % len(line_styles)
                ]
                subset_line_style_index[subset_type] += 1
            else:
                color = distinct_colors[i % len(distinct_colors)]
                line_style = 1

            cloned_hist.SetMarkerColor(color)
            cloned_hist.SetLineStyle(line_style)
            cloned_hist.SetLineColor(color)
            print(
                f"Set color {color} and line style {line_style} for histogram '{unique_hist_name}'"
            )

            # Apply axis labels from dictionary if available
            # if hist_name in axis_labels:
            #     x_label, y_label = axis_labels[hist_name]
            #     cloned_hist.GetXaxis().SetTitle(x_label)
            #     cloned_hist.GetYaxis().SetTitle(y_label)
            #     print(
            #         f"Set axis labels for '{hist_name}': X='{x_label}', Y='{y_label}'"
            #     )

            cloned_hist.GetXaxis().SetTitleSize(0.04)
            cloned_hist.GetYaxis().SetTitleSize(0.04)

            # Add to list of histograms for legend
            histograms_with_legend.append((cloned_hist, legend_entry))

            # Draw histogr    ROOT.gSystem.Sleep(5000)

            # draw_option = "E SAME" if not first_histogram else "E"
            if first_histogram:
                cloned_hist.Draw("E")
            else:
                cloned_hist.Draw("SAME E")

            # print(f"Drawing histogram '{unique_hist_name}' with option '{draw_option}'")
            # cloned_hist.Draw(draw_option)
            canvas.Update()
            ROOT.gSystem.Sleep(5000)

            first_histogram = False
        else:
            error_log[file_path].append(f"Histogram '{hist_name}' not found in file")
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
    output_file = f"/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/git reps/usercode/download/merge_scripts/30.10/imgout/combined_{hist_name}_overlap_histograms.pdf"
    canvas.SaveAs(output_file)
    print(f"Saved histogram as '{output_file}'")

# Print formatted error report at the end
print("\nError Summary:")
for file_path, errors in error_log.items():
    error_text = "None" if not errors else "\n  - " + "\n  - ".join(errors)
    print(f"{file_path}:\n  Errors: {error_text}")
