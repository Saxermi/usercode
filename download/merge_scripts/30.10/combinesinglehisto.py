import ROOT
import os
import re

# Define the input root files
root_files = [
    "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_16/TTbar_overlap_40_blocksize_512.root",
]

# List of histograms to combine
list_histos = ["PUBlockBordersvsPurityprofile1", "PURandomBlockBordersvsPurityprofile1", "PUDeterBlockBordersvsPurityprofile1"]

# Save folder
savefolder = "test"

# Ensure the save directory exists
output_path = f"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/combinedHistograms/{savefolder}"
os.makedirs(output_path, exist_ok=True)

# Distinct colors for histograms
distinct_colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2]

# Set global ROOT styles
ROOT.gStyle.SetOptStat(0)  # Disable statistics box
ROOT.gStyle.SetTextFont(42)  # Set text font to match CMS style

# Function to find histogram recursively
def find_histogram_in_directory(directory, hist_name):
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TH1") or obj.InheritsFrom("TProfile") or obj.InheritsFrom("TH2"):
            if obj.GetName() == hist_name:
                return obj
        elif obj.InheritsFrom("TDirectory"):
            found_hist = find_histogram_in_directory(obj, hist_name)
            if found_hist:
                return found_hist
    return None

# Process each ROOT file
for file_path in root_files:
    print(f"Processing file: {file_path}")
    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print(f"Failed to open file: {file_path}")
        continue

    directory = root_file.GetDirectory("testVertices_test/efficiency")
    if not directory:
        print(f"Directory 'testVertices_test/efficiency' not found in file: {file_path}")
        root_file.Close()
        continue

    # Create a canvas
    canvas = ROOT.TCanvas("canvas", "", 800, 700)

    # Add a legend
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
    legend.SetTextFont(42)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    first_histogram = True
    histograms_with_legend = []

    for i, hist_name in enumerate(list_histos):
        hist = find_histogram_in_directory(directory, hist_name)
        if hist:
            hist.SetDirectory(0)

            # Clone the histogram with a unique name
            unique_hist_name = f"{hist.GetName()}_{os.path.basename(file_path)}"
            cloned_hist = hist.Clone(unique_hist_name)
            cloned_hist.SetDirectory(0)

            # Style the histogram
            color = distinct_colors[i % len(distinct_colors)]
            cloned_hist.SetLineColor(color)
            cloned_hist.SetLineStyle(1)
            cloned_hist.SetMarkerColor(color)

            # Add to legend
            legend.AddEntry(cloned_hist, hist_name, "l")

            # Draw the histogram
            if first_histogram:
                cloned_hist.SetMinimum(0)
                cloned_hist.SetMaximum(100)
                cloned_hist.Draw("E")
                first_histogram = False
            else:
                cloned_hist.Draw("SAME E")

            histograms_with_legend.append(cloned_hist)
        else:
            print(f"Histogram '{hist_name}' not found in file: {file_path}")

    # Draw legend
    legend.Draw()

    # Save the canvas
    base_filename = os.path.basename(file_path).replace(".root", "")
    output_file = f"{output_path}/{base_filename}_combined_histograms.pdf"
    canvas.SaveAs(output_file)
    print(f"Saved combined histogram as: {output_file}")

    root_file.Close()
