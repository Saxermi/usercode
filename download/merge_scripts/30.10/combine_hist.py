import ROOT
import re
import os

# List of histograms to retrieve and plot
histogram_names = [
    "SERecoVsSimZPosition",
    "PURecoVsSimZPosition",
    "SETracksPurity",
    "SETracksEfficiency",
    "PUTracksPurity",
    "PUTracksEfficiency",
    "SEResolution",
    "PUResolution",
    "reco_vs_true_z_position_hist_categorial_c1",
    "reco_vs_true_z_position_hist_categorial_c2",
    "reco_vs_true_z_position_hist_categorial_c3",
    "True_3D_point_to_plane_distance",
    "SE_reco_index_hist",
    "SE_reco_index_histHR",
    "PUTracksPurityBlock",
    "SETracksPurityBlock",
    "NVertexVSCPUTime",
    "SERecoVsSimZPositionBlock",
    "PURecoVsSimZPositionBlock",
    "PUBlockBordersvsPurityprofile",
    "PUBlockBordersvsPurity",
    "PUBlockBordersvsEfficencyprofile",
    "PUBlockBordersvsEfficency",
    "PUBlockBordersvsZdeltayprofile",
    "PUBlockBordersvsZdelta",
    "PUPurityVsZaxisprofile",
    "PUPurityVsZaxis",
    "SEEfficiencyVsZaxisProfile",
    "SEEfficiencyVsZaxis",
    "PUEfficiencyVsZaxisProfile",
    "PUEfficiencyVsZaxis",
    "SEResolutionNormalized",
    "SEResolutionNormalizedBlockprofile",
    "SEResolutionNormalizedBlock",
    "PUResolutionNormalized",
    "SEResolutionVsTrackPurity",
    "SEBlockBordersvsPurityprofile",
    "SEBlockBordersvsPurity",
    "SEPurityVsZaxisProfile",
    "SEPurityVsZaxis",
    "SEBlockBordersvsEfficencyprofile",
    "SEBlockBordersvsEfficency",
    "PUBlockBordersvsFakeVertProfi",
    "PUBlockBorder",
    "SEBlockBorder",
    "BlockSizes",
    "BlockNumber",
    "SESimulatedVertices",
    "PUSimulatedVertices",
    "SEReconVertices",
    "PUReconVertices",
    "FakeVertices"
]

# List of ROOT files with full path
root_files = [
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_30_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_30_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_40_blocksize_256.root",
    "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_40_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_50_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_50_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_50_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_30_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_30_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_40_blocksize_256.root",
    "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_40_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_50_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_50_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_50_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_30_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_30_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_40_blocksize_256.root",
    "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_40_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_50_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_50_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_50_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_30_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_30_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_40_blocksize_256.root",
    "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_40_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_50_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_50_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_50_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_6/HiggsGluonFusion_overlap_n40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_6/SToMuMu_overlap_n40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_6/TTbar_overlap_n40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_6/ZMM_overlap_n40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_7/experimental_run_7/HiggsGluonFusion_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_7/experimental_run_7/SToMuMu_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_7/experimental_run_7/TTbar_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_7/experimental_run_7/ZMM_overlap_40_blocksize_512.root"

]

# Colors assigned to each subset if grouping by dataset
subset_colors = {
    "HiggsGluonFusion": ROOT.kRed,
    "TTbar": ROOT.kBlue,
    "ZMM": ROOT.kGreen + 2,
    "SToMuMu": ROOT.kOrange,
    "experimental_run_7": ROOT.kGray  # Set experimental_run_7 files to gray color
}
savefolder = "comparisonSE"
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

# Loop through each histogram name
for hist_name in histogram_names:
    # Prepare a legend and list for each histogram type
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

    # Track line style usage within each subset
    subset_line_style_index = {}
    histograms_with_legend = []

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
        elif "experimental_run_7" in file_path:
            legend_entry = f"{filename.split('_')[0]} no blocks"
        else:
            legend_entry = filename
            print(
                f"Filename {filename} does not match pattern. Using filename as legend entry."
            )

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
            if "experimental_run_7" in file_path:
                color = ROOT.kGray
                line_style = 1
            elif group_by_dataset:
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
            cloned_hist.GetXaxis().SetTitleSize(0.04)
            cloned_hist.GetYaxis().SetTitleSize(0.04)

            # Add to list of histograms for legend
            histograms_with_legend.append((cloned_hist, legend_entry))

            # Draw histogram
            if first_histogram:
                cloned_hist.Draw("E")
            else:
                cloned_hist.Draw("SAME E")

            canvas.Update()
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
    output_file = f"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/combinedHistograms/{savefolder}/combined_{hist_name}_blocksizes_histograms.pdf"
    canvas.SaveAs(output_file)
    print(f"Saved histogram as '{output_file}'")

# Additional section to divide SESimulatedVertices by SEReconVertices, and by (SEReconVertices + FakeVertices)
# and also for PUSimulatedVertices by PUReconVertices and (PUReconVertices + FakeVertices)
for file_path in root_files:
    root_file = ROOT.TFile.Open(file_path)
    if root_file and not root_file.IsZombie():
        directory = root_file.GetDirectory("testVertices_test/efficiency")
        if directory:
            # SE Histograms
            simulated_hist = find_histogram_in_directory(directory, "SESimulatedVertices")
            recon_hist = find_histogram_in_directory(directory, "SEReconVertices")
            fake_hist = find_histogram_in_directory(directory, "FakeVertices")

            if simulated_hist and recon_hist:
                # SEReconVertices / SESimulatedVertices
                ratio_hist_6 = recon_hist.Clone(f"ratio_SERecon_vs_SESimulated_{os.path.basename(file_path)}")
                ratio_hist_6.Divide(simulated_hist)
                ratio_hist_6.SetDirectory(0)
                ratio_hist_6.SetTitle("SEReconVertices / SESimulatedVertices")
                ratio_hist_6.SetLineColor(ROOT.kGray)
                ratio_hist_6.Draw("E")

                # Create legend for SERecon / SESimulated
                se_legend_recon = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
                se_legend_recon.AddEntry(ratio_hist_6, "SEReconVertices / SESimulatedVertices", "lep")
                se_legend_recon.Draw()
                canvas.Update()
                output_file = f"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/combinedHistograms/{savefolder}/ratio_SERecon_vs_SESimulated_{os.path.basename(file_path)}.pdf"
                canvas.SaveAs(output_file)
                print(f"Saved ratio histogram as '{output_file}'")

                # FakeVertices / SESimulatedVertices
                ratio_hist_7 = fake_hist.Clone(f"ratio_Fake_vs_SESimulated_{os.path.basename(file_path)}")
                ratio_hist_7.Divide(simulated_hist)
                ratio_hist_7.SetDirectory(0)
                ratio_hist_7.SetTitle("FakeVertices / SESimulatedVertices")
                ratio_hist_7.SetLineColor(ROOT.kRed)
                ratio_hist_7.Draw("E SAME")

                # Create legend for FakeVertices ratio
                fake_legend_se = ROOT.TLegend(0.7, 0.5, 0.9, 0.7)
                fake_legend_se.AddEntry(ratio_hist_7, "FakeVertices / SESimulatedVertices", "lep")
                fake_legend_se.Draw()
                canvas.Update()
                output_file = f"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/combinedHistograms/{savefolder}/ratio_Fake_vs_SESimulated_{os.path.basename(file_path)}.pdf"
                canvas.SaveAs(output_file)
                print(f"Saved ratio histogram as '{output_file}'")
                canvas.Update()
                output_file = f"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/combinedHistograms/{savefolder}/ratio_SESimulated_vs_SERecon_{os.path.basename(file_path)}.pdf"
                canvas.SaveAs(output_file)
                print(f"Saved ratio histogram as '{output_file}'")


            # PU Histograms
            pu_simulated_hist = find_histogram_in_directory(directory, "PUSimulatedVertices")
            pu_recon_hist = find_histogram_in_directory(directory, "PUReconVertices")

            if pu_simulated_hist and pu_recon_hist:
                ratio_hist_3 = pu_simulated_hist.Clone(f"ratio_PUSimulated_vs_PURecon_{os.path.basename(file_path)}")
                ratio_hist_3.Divide(pu_recon_hist)
                ratio_hist_3.SetDirectory(0)
                ratio_hist_3.SetTitle("PUSimulatedVertices / PUReconVertices")
                ratio_hist_3.Draw("E")

                # Create legend
                pu_legend_3 = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
                pu_legend_3.AddEntry(ratio_hist_3, "PUSimulatedVertices / PUReconVertices", "lep")
                pu_legend_3.Draw()
                canvas.Update()
                output_file = f"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/combinedHistograms/{savefolder}/ratio_PUSimulated_vs_PURecon_{os.path.basename(file_path)}.pdf"
                canvas.SaveAs(output_file)
                print(f"Saved ratio histogram as '{output_file}'")

            if pu_simulated_hist and pu_recon_hist and fake_hist:
                # PUReconVertices / PUSimulatedVertices
                ratio_hist_4 = pu_recon_hist.Clone(f"ratio_PURecon_vs_PUSimulated_{os.path.basename(file_path)}")
                ratio_hist_4.Divide(pu_simulated_hist)
                ratio_hist_4.SetDirectory(0)
                ratio_hist_4.SetTitle("PUReconVertices / PUSimulatedVertices")
                ratio_hist_4.SetLineColor(ROOT.kGray)
                ratio_hist_4.Draw("E")

                # Create legend
                pu_legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
                pu_legend.AddEntry(ratio_hist_4, "PUReconVertices / PUSimulatedVertices", "lep")
                pu_legend.Draw()
                canvas.Update()
                output_file = f"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/combinedHistograms/{savefolder}/ratio_PURecon_vs_PUSimulated_{os.path.basename(file_path)}.pdf"
                canvas.SaveAs(output_file)
                print(f"Saved ratio histogram as '{output_file}'")

                # FakeVertices / PUSimulatedVertices
                ratio_hist_5 = fake_hist.Clone(f"ratio_Fake_vs_PUSimulated_{os.path.basename(file_path)}")
                ratio_hist_5.Divide(pu_simulated_hist)
                ratio_hist_5.SetDirectory(0)
                ratio_hist_5.SetTitle("FakeVertices / PUSimulatedVertices")
                ratio_hist_5.SetLineColor(ROOT.kRed)
                ratio_hist_5.Draw("E SAME")

                # Create legend for FakeVertices ratio
                fake_legend = ROOT.TLegend(0.7, 0.5, 0.9, 0.7)
                fake_legend.AddEntry(ratio_hist_5, "FakeVertices / PUSimulatedVertices", "lep")
                fake_legend.Draw()
                canvas.Update()
                output_file = f"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/combinedHistograms/{savefolder}/ratio_Fake_vs_PUSimulated_{os.path.basename(file_path)}.pdf"
                canvas.SaveAs(output_file)
                print(f"Saved ratio histogram as '{output_file}'")

                # Create legend
                pu_legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
                pu_legend.AddEntry(ratio_hist_4, "PUSimulatedVertices / (PUReconVertices + FakeVertices)", "lep")
                pu_legend.Draw()
                canvas.Update()
                output_file = f"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/combinedHistograms/{savefolder}/ratio_PUSimulated_vs_PUReconPlusFake_{os.path.basename(file_path)}.pdf"
                canvas.SaveAs(output_file)
                print(f"Saved ratio histogram as '{output_file}'")

        root_file.Close()

# Print formatted error report at the end
print("\nError Summary:")
for file_path, errors in error_log.items():
    error_text = "None" if not errors else "\n  - " + "\n  - ".join(errors)
    print(f"{file_path}:\n  Errors: {error_text}")
