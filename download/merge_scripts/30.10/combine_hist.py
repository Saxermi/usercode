import ROOT
import re
import os

# List of histograms to retrieve and plot
histogram_names = [
    # "SERecoVsSimZPosition",
    # "PURecoVsSimZPosition",
    # "SETracksPurity",
    # "SETracksEfficiency",
    # "PUTracksPurity",
    # "PUTracksEfficiency",
    "SEResidual",
    "PUResidual",
    # "reco_vs_true_z_position_hist_categorial_c1",
    # "reco_vs_true_z_position_hist_categorial_c2",
    # "reco_vs_true_z_position_hist_categorial_c3",
    # "True_3D_point_to_plane_distance",
    # "SE_reco_index_hist",
    # "SE_reco_index_histHR",
    # "PUTracksPurityBlock",
    # "SETracksPurityBlock",
    # "NVertexVSCPUTime",
    # "SERecoVsSimZPositionBlock",
    # "PURecoVsSimZPositionBlock",
    # "PUBlockBordersvsPurityprofile",
    # "PURandomBlockBordersvsPurityprofile",
    # "PUBlockBordersvsPurity",
    # "PUBlockBordersvsEfficencyprofile",
    # "PUBlockBordersvsEfficency",
    # "PURandomBlockBordersvsEfficency",
    # "PUDeterBlockBordersvsPurityprofile",
    # "PUDeterBlockBordersvsPurity",
    # "PUDeterBlockBordersvsEfficencyprofile",
    # "PUDeterBlockBordersvsEfficency",
    "PUBlockBordersvsPurityprofile1",
    # "PUBlockBordersvsPurity1",
    # "PURandomBlockBordersvsPurity1",
    "PURandomBlockBordersvsPurityprofile1",
    # "PURandomBlockBordersvsPurity5",
    "PURandomBlockBordersvsPurityprofile5",
    "PUBlockBordersvsEfficencyprofile1",
    # "PURandomBlockBordersvsEfficency1",
    "PURandomBlockBordersvsEfficencyprofile1",
    # "PURandomBlockBordersvsEfficency5",
    # "PUBlockBordersvsEfficency5",
    "PUBlockBordersvsPurityprofile5",
    # "PUBlockBordersvsPurity5",
    "PUBlockBordersvsEfficencyprofile5",
    "PUDeterBlockBordersvsPurityprofile1",
    # "PUDeterBlockBordersvsPurity1",
    "PUDeterBlockBordersvsPurityprofile5",
    # "PUDeterBlockBordersvsPurity5",
    "PUDeterBlockBordersvsEfficencyprofile1",
    # "PUDeterBlockBordersvsEfficency1",
    "PUDeterBlockBordersvsEfficencyprofile5",
    # "PUDeterBlockBordersvsEfficency5",
    # "PUBlockBordersvsZdeltayprofile",
    # "PUBlockBordersvsZdelta",
    # "PUPurityVsZaxisprofile",
    # "PUPurityVsZaxis",
    # "PUPurityVsNumTracks",
    # "SEEfficiencyVsZaxisProfile",
    # "SEEfficiencyVsZaxis",
    # "SEEfficiencyVsNumTracks",
    # "PUEfficiencyVsZaxisProfile",
    # "PUEfficiencyVsZaxis",
    # "PUEfficiencyVsNumTracks",
    # "SEResidualNormalized",
    # "SEResidualNormalizedBlockprofile",
    # "SEResidualNormalizedBlock",
    # "PUResidualNormalized",
    # "SEResidualVsTrackPurity",
    # "SEResidualVsTrackPurityprofile",
    # "SEBlockBordersvsPurityprofile",
    # "SEBlockBordersvsPurity",
    # "SEBlockBordersvsEfficencyprofile",
    # "SEBlockBordersvsEfficency",
    # "SERandomBlockBordersvsPurityprofile",
    # "SERandomBlockBordersvsPurity",
    # "SEPurityVsNumTracks",
    # "SERandomBlockBordersvsEfficencyprofile",
    # "SERandomBlockBordersvsEfficency",
    # "SEDeterBlockBordersvsPurityprofile",
    # "SEDeterBlockBordersvsPurity",
    # "SEDeterBlockBordersvsEfficencyprofile",
    # "SEDeterBlockBordersvsEfficency",
    # "SEPurityVsZaxisProfile",
    # "SEPurityVsZaxis",
    "SEBlockBordersvsPurityprofile1",
    # "SEBlockBordersvsPurity1",
    "SEBlockBordersvsPurityprofile5",
    # "SEBlockBordersvsPurity5",
    "SEBlockBordersvsEfficencyprofile1",
    # "SEBlockBordersvsEfficency1",
    "SEBlockBordersvsEfficencyprofile5",
    # "SEBlockBordersvsEfficency5",
    "SERandomBlockBordersvsPurityprofile1",
    # "SERandomBlockBordersvsPurity1",
    "SERandomBlockBordersvsPurityprofile5",
    # "SERandomBlockBordersvsPurity5",
    "SERandomBlockBordersvsEfficencyprofile1",
    # "SERandomBlockBordersvsEfficency1",
    "SERandomBlockBordersvsEfficencyprofile5",
    # "SERandomBlockBordersvsEfficency5",
    "SEDeterBlockBordersvsPurityprofile1",
    # "SEDeterBlockBordersvsPurity1",
    "SEDeterBlockBordersvsPurityprofile5",
    # "SEDeterBlockBordersvsPurity5",
    "SEDeterBlockBordersvsEfficencyprofile1",
    # "SEDeterBlockBordersvsEfficency1",
    "SEDeterBlockBordersvsEfficencyprofile5",
    # "SEDeterBlockBordersvsEfficency5",
    # "PUBlockBordersvsFakeVertProfi",
    # "PUBlockBordersvsFakeVertProfi1",
    # "PUBlockBordersvsFakeVertProfi5",
    # "PUBlockBorder",
    # "PUBlockBorder1",
    # "SEBlockBorder",
    "BlockSizes",
    "BlockNumber",
    # "SESimulatedVertices",
    # "PUSimulatedVertices",
    # "SEReconVertices",
    # "PUReconVertices",
    # "FakeVertices",
    # "PUSimVertexTrackDist",
    "PUSimVertexTrackDistLog",
    # "PUReconVertexTrackDist",
    "PUReconVertexTrackDistLog",
    # "PUFakeVertexTrackDist",
    "PUFakeVertexTrackDistLog",
    # "SESimVertexTrackDist",
    "SESimVertexTrackDistLog",
    # "SEReconVertexTrackDist",
    "SEReconVertexTrackDistLog",
    "PUSimNumTracksZPos",
    "PUReconNumTracksZPos",
    "PUFakeNumTracksZPos",
    "SESimNumTracksZPos",
    "SEReconNumTracksZPos",
    # "PUSimNumTracksBlock",
    # "PUReconNumTracksBlock",
    # "PUFakeNumTracksBlock",
    # "SESimNumTracksBlock",
    # "SEReconNumTracksBlock",
    "PUSimNumTracksBlock1",
    "PUReconNumTracksBlock1",
    "PUFakeNumTracksBlock1",
    "SESimNumTracksBlock1",
    "SEReconNumTracksBlock1",
]



# List of ROOT files with full path
root_files = [
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_30_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_30_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_40_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_40_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_50_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_50_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/HiggsGluonFusion_overlap_50_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_30_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_30_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_40_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_40_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_50_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_50_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/SToMuMu_overlap_50_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_30_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_30_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_40_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_40_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_50_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_50_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/TTbar_overlap_50_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_30_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_30_blocksize_1024.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_40_blocksize_256.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_5/ZMM_overlap_40_blocksize_512.root",
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
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_7/experimental_run_7/ZMM_overlap_40_blocksize_512.root",
    #"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_15/TTbar_overlap_0_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_15/TTbar_overlap_0_blocksize_1024.root",
    #"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_15/TTbar_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_15/TTbar_overlap_40_blocksize_1024.root",
    #"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_15/TTbar_overlap_n40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_15/TTbar_overlap_n40_blocksize_1024.root",
    #"/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_16/TTbar_overlap_40_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_31_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_32_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_33_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_34_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_35_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_36_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_37_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_38_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_39_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/rootFileStorage/experimental_run_17/TTbar_overlap_30_blocksize_512.root",
    # "/t3home/frejalom/cmssw/CMSSW_14_1_0_pre7/src/usercode/PrimaryVertexAnalyzer/test/pv.root"
#"/t3home/msaxer/cmssw/CMSSW_14_1_0_pre7/src/usercode/PrimaryVertexAnalyzer/test/pv.root"
"/work/msaxer/ba/experimental_run_1/Subset_TTbar_01/40/512/pvSubset_TTbar_01_OL0.4_BS512_20250228_212321.root"
]

# Colors assigned to each subset if grouping by dataset
subset_colors = {
    "HiggsGluonFusion": ROOT.kRed,
    "TTbar": ROOT.kBlue,
    "ZMM": ROOT.kGreen + 2,
    "SToMuMu": ROOT.kOrange,
    "experimental_run_7": ROOT.kGray,  # Set experimental_run_7 files to gray color
    "experimental_run_16": ROOT.kGray,
}
savefolder = "quicktest"
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
pattern = re.compile(r"([A-Za-z]+)_overlap_(n?\d+)_blocksize_(\d+)\.root")

# # Axis label dictionary for each histogram
# axis_labels = {
#     "PUBlockBordersvsZdeltayprofile": (
#         "Distance to nearest block (mm)",
#         "Delta Y (mm)",
#     ),
#     # Add other axis labels as needed
# }

# Add ratio histograms to compare across files
ratio_names = [
    "SEReconVertices/SESimulatedVertices",
    "PUReconVertices/PUSimulatedVertices",
    "FakeVertices/PUSimulatedVertices",
]

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
for hist_name in histogram_names + ratio_names:  # Include ratios in the same loop
    # Prepare a legend and list for each histogram type
    legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.3)

    # Track line style usage within each subset
    subset_line_style_index = {}
    histograms_with_legend = []

    # Clear the canvas before drawing new histograms
    canvas.Clear()

    # Loop through the files and draw histograms
    first_histogram = True
    for i, file_path in enumerate(root_files):
        # Initialize or retrieve the error list for each file
        error_log.setdefault(file_path, [])

        # Open ROOT file
        root_file = ROOT.TFile.Open(file_path)
        if not root_file or root_file.IsZombie():
            error_log[file_path].append("Failed to open file")
            continue  # Skip this file and move to the next

        # Navigate to the correct directory within the ROOT file
        directory = root_file.GetDirectory("testVertices_test/efficiency")
        if not directory:
            error_log[file_path].append(f"Directory not found in file: {file_path}")
            root_file.Close()
            continue

        # Extract filename and match against the regex
        filename = os.path.basename(file_path)
        match = pattern.search(filename)
        if ("experimental_run_7" in file_path) or ("experimental_run_16" in file_path):
            legend_entry = f"{filename.split('_')[0]} no blocks"
        elif match:
            subset_type = match.group(1)  # Process name (e.g., HiggsGluonFusion, TTbar)
            overlap = match.group(2)  # Overlap (e.g., n40, 40)
            blocksize = int(match.group(3))  # Block size (e.g., 512)

            legend_entry = f"{subset_type} overlap {overlap} blocksize {blocksize}"
        else:
            legend_entry = filename


        # Handle regular histograms
        if hist_name not in ratio_names:
            hist = find_histogram_in_directory(directory, hist_name)
        # Handle ratio histograms
        else:
            if hist_name == "SEReconVertices/SESimulatedVertices":
                numerator = find_histogram_in_directory(directory, "SEReconVertices")
                denominator = find_histogram_in_directory(directory, "SESimulatedVertices")
            elif hist_name == "PUReconVertices/PUSimulatedVertices":
                numerator = find_histogram_in_directory(directory, "PUReconVertices")
                denominator = find_histogram_in_directory(directory, "PUSimulatedVertices")
            elif hist_name == "FakeVertices/PUSimulatedVertices":
                numerator = find_histogram_in_directory(directory, "FakeVertices")
                denominator = find_histogram_in_directory(directory, "PUSimulatedVertices")
            else:
                continue

            if numerator and denominator:
                numerator.Rebin(5)
                denominator.Rebin(5)
                hist = numerator.Clone(f"{hist_name}_{os.path.basename(file_path)}")
                hist.Divide(numerator, denominator, 1.0, 1.0, "B")
            else:
                error_log[file_path].append(f"Missing numerator/denominator for {hist_name}")
                root_file.Close()
                continue

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
                continue

            # Check flag and set color/line style accordingly
            if ("experimental_run_7" in file_path) or ("experimental_run_16" in file_path):
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

            # Add to list of histograms for legend
            histograms_with_legend.append((cloned_hist, legend_entry))

            # Draw histogram
            if first_histogram:
                if ("profile" in hist_name):
                    cloned_hist.SetMinimum(0)  # Set Y-axis minimum
                    cloned_hist.SetMaximum(100)  # Set Y-axis maximum
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

    # Update and save canvas
    safe_hist_name = hist_name.replace(" ", "_").replace("/", "_")
    output_file = f"/t3home/msaxer/{savefolder}/combined_{safe_hist_name}_blocksizes_histograms.pdf"
    canvas.SaveAs(output_file)
    print(f"Saved histogram as '{output_file}'")

