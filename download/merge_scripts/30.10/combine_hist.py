import ROOT
import re
import os

# List of ROOT files with full path
root_files = [
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_HiggsGluonFusion_01_512_20241024_052146.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_HiggsGluonFusion_01_n512_20241025_181727.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_TTbar_01_512_20241024_063454.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_TTbar_01_n512_20241025_180817.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_ZMM_01_512_20241024_064931.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/30.10/experimental_run_6/pvSubset_ZMM_01_n512_20241025_184510.root",
]

# Colors to use for each file's error bars
colors = [
    ROOT.kRed,
    ROOT.kBlue,
    ROOT.kGreen,
    ROOT.kMagenta,
    ROOT.kCyan,
    ROOT.kOrange,
    ROOT.kViolet,
    ROOT.kYellow,
    ROOT.kPink,
]

# Set global ROOT styles
ROOT.gStyle.SetOptStat(0)  # Disable statistics box
ROOT.gStyle.SetTextFont(42)  # Set text font to match CMS style

# Create a canvas
canvas = ROOT.TCanvas("canvas", "", 800, 700)

# List to store histograms and their legend entries
histograms_with_legend = []

# Variable to track if the first histogram is drawn
first_histogram = True

# Regex pattern to match and extract parts of the filename
pattern = re.compile(r"pvSubset_([A-Za-z]+)_(\d+)_(n?)(\d+)_\d{8}_\d{6}\.root")

# Loop through the files and draw histograms
for i, file_path in enumerate(root_files):
    print(f"Opening ROOT file: {file_path}")

    # Open ROOT file
    root_file = ROOT.TFile.Open(file_path)

    # Check if the file opened successfully
    if root_file and not root_file.IsZombie():
        print(f"Successfully opened: {file_path}")
    else:
        print(f"Failed to open file: {file_path}")
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

    # Check if the directory exists
    if not efficiency_dir:
        print(
            f"Error: 'offlinePrimaryVertices/efficiency' directory not found in {file_path}!"
        )
        # Close the ROOT file to free resources
        root_file.Close()
        continue  # Skip this file and move to the next

    # Retrieve the TProfile histogram
    hist = efficiency_dir.Get("PUBlockBordersvsPurityprofile")

    # Check if the histogram was retrieved successfully
    if hist:
        print(f"Successfully retrieved histogram from {file_path}")

        # Prevent the histogram from being deleted when the file is closed
        hist.SetDirectory(0)

        # Set the color for the markers and error bars
        hist.SetMarkerColor(colors[i % len(colors)])
        hist.SetLineStyle(0)
        hist.SetLineColor(colors[i % len(colors)])

        # Set X and Y axis labels with appropriate size
        hist.GetXaxis().SetTitle("Distance to nearest block (mm)")
        hist.GetYaxis().SetTitle("Purity (%)")
        hist.GetXaxis().SetTitleSize(0.04)
        hist.GetYaxis().SetTitleSize(0.04)

        # Add the histogram and its legend entry to the list
        histograms_with_legend.append((hist, legend_entry))

        # Draw the TProfile histogram with error bars only
        draw_option = "E" if first_histogram else "E SAME"
        if first_histogram:
            first_histogram = False
        hist.Draw(draw_option)  # Draw on the canvas
    else:
        print(
            f"Error: TProfile histogram 'PUBlockBordersvsPurityprofile' not found in {file_path}"
        )
        # Close the ROOT file to free resources
        root_file.Close()
        continue  # Skip this file and move to the next

    # Close the ROOT file to free resources
    root_file.Close()

# Add a legend
print("Adding a legend")
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
for hist, legend_entry in histograms_with_legend:
    legend.AddEntry(hist, legend_entry, "lep")  # Use stored entry text

legend.Draw()

# Add labels for CMS style
label = ROOT.TLatex()
label.SetNDC(True)
label.SetTextSize(0.040)
label.DrawLatex(0.100, 0.920, "#bf{CMS Simulation}")
label.DrawLatex(0.550, 0.920, "#sqrt{s} = 14 TeV")

# Update the canvas to show everything
canvas.Update()

# Save the plot as a .png file to the specified directory
output_file = "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/git reps/usercode/download/merge_scripts/30.10/combined_PUBlockBordesvsPurity_overlap_histograms.png"
canvas.SaveAs(output_file)

