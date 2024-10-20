import os
import ROOT
import re

# List of ROOT files with full path
root_files = [
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/TTBar1_15_overlap_10.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/TTBar1_15_overlap_20.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/TTBar1_15_overlap_30.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/TTBar1_15_overlap_40.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/TTBar1_15_overlap_50.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/TTBar1_15_overlap_60.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/TTBar1_15_overlap_70.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/TTBar1_15_overlap_80.root",
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/TTBar1_15_overlap_90.root",
]

# List of relevant histogram names to be processed
relevant_hist_names = [
    "PUBlockBordersvsPurityprofile",
    "PUBlockBordersvsEfficencyprofile",
    "PUBlockBordersvsZdeltayprofile",
    "NVertexVSCPUTime",
]

# Dictionary to store statistics for PUBlockBordersvsPurityprofile
PUBlockBordersvsPurityprofile_stats = []

# Use this pattern to extract overlap values from file paths
overlap_pattern = re.compile(r"overlap_\d+")

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
print("Setting global ROOT styles")
ROOT.gStyle.SetOptStat(0)  # Disable statistics box
ROOT.gStyle.SetTextFont(42)  # Set text font to match CMS style

# Dictionary to store histograms information
hist_dict = {}


# Function to add CMS style labels to the canvas
def add_cms_labels(canvas, text1="CMS Simulation", text2="#sqrt{s} = 14 TeV"):
    print("Adding CMS simulation labels")
    label = ROOT.TLatex()
    label.SetNDC(True)  # Set to normalized device coordinates
    label.SetTextSize(0.040)
    label.DrawLatex(0.100, 0.920, f"#bf{{{text1}}}")
    label.DrawLatex(0.550, 0.920, f"{text2}")
    canvas.Update()


# Check if 'images' directory exists, create it if not
def ensure_images_directory():
    images_dir = "images"
    if not os.path.exists(images_dir):
        print(f"Creating directory: {images_dir}")
        os.makedirs(images_dir)
    return images_dir


# Function to save the canvas
def save_canvas(canvas, file_name):
    images_dir = ensure_images_directory()
    output_path = os.path.join(images_dir, file_name)
    print(f"Saving canvas as {output_path}")
    canvas.SaveAs(output_path)


def extract_run_name(file_path):
    # Normalize the path to ensure consistent separators
    path_parts = os.path.normpath(file_path).split(os.sep)
    # Look for the part of the path that contains 'experimental_run_'
    for part in path_parts:
        if "experimental_run_" in part:
            return part
    # Return a default value if not found
    return "run_unknown"


# Loop over each ROOT file and store histograms into a dictionary
for i, file_path in enumerate(root_files):
    print(f"Opening ROOT file: {file_path}")
    root_file = ROOT.TFile.Open(file_path)

    if not root_file or root_file.IsZombie():
        print(f"Failed to open file: {file_path}")
        continue

    # Get histograms from 'offlinePrimaryVertices/efficiency' directory
    efficiency_dir = root_file.Get("offlinePrimaryVertices/efficiency")
    if not efficiency_dir:
        print(f"Directory 'offlinePrimaryVertices/efficiency' not found in {file_path}")
        continue

    # Loop through all histograms in the directory
    for key in efficiency_dir.GetListOfKeys():
        obj = key.ReadObj()

        # Ensure we process 1D, 2D, or TProfile histograms
        if (
            isinstance(obj, ROOT.TH1)
            or isinstance(obj, ROOT.TH2)
            or isinstance(obj, ROOT.TProfile)
        ):
            hist_name = obj.GetName()
            hist_type = type(obj).__name__

            print(f"Found histogram: {hist_name}, Type: {hist_type}")

            # Detach the histogram from the file so it can be used after closing the file
            obj.SetDirectory(0)

            # Store histogram information in the dictionary
            if hist_name not in hist_dict:
                hist_dict[hist_name] = []

            # Append to the list of histograms for this name
            hist_dict[hist_name].append(
                {
                    "file": file_path,
                    "type": hist_type,
                    "histogram": obj,
                    "color": colors[i % len(colors)],  # Assign color to histogram
                }
            )

    # Close the ROOT file
    root_file.Close()

# Create a canvas to draw on
canvas = ROOT.TCanvas("canvas", "", 800, 700)

# Processing each diagram (histogram) in the dictionary
for hist_name, hist_info_list in hist_dict.items():
    # Check if the diagram is in the list of relevant histograms
    if hist_name in relevant_hist_names:
        print(f"Processing relevant histograms for: {hist_name}")

        # Create a legend for the histograms
        legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

        first_hist = True  # Flag to determine drawing options

        # Process each histogram with the same name but from different files
        for info in hist_info_list:
            # Get the histogram object
            hist = info["histogram"]

            # Set the color for the histogram
            hist.SetMarkerColor(info["color"])
            hist.SetLineColor(info["color"])

            # Extract the overlap_* part from the file path
            match = overlap_pattern.search(info["file"])
            if match:
                overlap_value = match.group()
            else:
                overlap_value = "unknown"
            # Draw the first histogram and then overlay the rest
            draw_option = "E" if first_hist else "E SAME"
            first_hist = False

            # Histogram type and name check for specific logic
            if info["type"] == "TProfile":
                if hist_name == "PUBlockBordersvsPurityprofile":
                    print(f"Applying logic for {hist_name}")
                    hist.GetXaxis().SetTitle("Distance to nearest block (mm)")
                    hist.GetYaxis().SetTitle("Purity (%)")
                    hist.GetXaxis().SetTitleSize(0.04)
                    hist.GetYaxis().SetTitleSize(0.04)

                    # Access the number of entries (total number of data points in the histogram)
                    entries = hist.GetEntries()

                    # Compute the mean y-value manually
                    total_sum_y = 0
                    total_entries = 0
                    for bin in range(1, hist.GetNbinsX() + 1):
                        bin_entries = hist.GetBinEntries(bin)
                        bin_content = hist.GetBinContent(bin)
                        total_sum_y += bin_entries * bin_content
                        total_entries += bin_entries

                    mean_y = total_sum_y / total_entries if total_entries > 0 else 0

                    # Access the mean of the x-axis
                    mean_x = hist.GetMean()  # Mean along x-axis

                    # Access the standard deviation of the x-axis
                    std_dev_x = hist.GetStdDev()  # Std dev along x-axis

                    # Note: Computing std_dev_y is more complex and might not be directly needed

                    # Store the values in a dictionary and append to the list
                    PUBlockBordersvsPurityprofile_stats.append(
                        {
                            "overlap": overlap_value,
                            "hist_name": hist_name,
                            "entries": entries,
                            "mean_x": mean_x,
                            "mean_y": mean_y,
                            "std_dev_x": std_dev_x,
                            # "std_dev_y": std_dev_y,  # Optional
                        }
                    )
                    print(hist)
                    print(overlap_value)
                    print(f"Entries: {entries}")
                    print(f"Mean X: {mean_x}")
                    print(f"Mean Y: {mean_y}")
                    print(f"Std Dev X: {std_dev_x}")
                    # print(f"Std Dev Y: {std_dev_y}")

                elif hist_name == "PUBlockBordersvsEfficencyprofile":
                    print(f"Applying specific logic for {hist_name}")
                    hist.GetXaxis().SetTitle("Distance to nearest block (mm)")
                    hist.GetYaxis().SetTitle("Efficiency (%)")
                    hist.GetXaxis().SetTitleSize(0.04)
                    hist.GetYaxis().SetTitleSize(0.04)
                    # Add specific settings for this histogram if needed

                elif hist_name == "PUBlockBordersvsZdeltayprofile":
                    print(f"Applying specific logic for {hist_name}")
                    hist.GetXaxis().SetTitle("Distance to nearest block (mm)")
                    hist.GetYaxis().SetTitle("Delta z axis (mm)")
                    hist.GetXaxis().SetTitleSize(0.04)
                    hist.GetYaxis().SetTitleSize(0.04)

            elif info["type"] == "TH2":
                if hist_name == "NVertexVSCPUTime":
                    print(f"Processing 2D histogram: {hist_name}")

                    # Create a TGraph to store the mean points from all histograms
                    graph = ROOT.TGraph()

                    # Initialize the legend
                    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

                    point_index = 0  # To index points in the graph
                    for info in hist_info_list:
                        # Get the histogram object
                        hist = info["histogram"]

                        # Extract the mean X and Y values (cluster center)
                        mean_x = hist.GetMean(1)  # Mean along X (Number of vertices)
                        mean_y = hist.GetMean(2)  # Mean along Y (CPU time)

                        print(
                            f"File: {info['file']} - Mean X: {mean_x}, Mean Y: {mean_y}"
                        )

                        # Set point in the graph for each file (mean cluster positions)
                        graph.SetPoint(point_index, mean_x, mean_y)

                        # Set marker color and style for each histogram
                        graph.SetMarkerColor(
                            info["color"]
                        )  # Set the color for each point
                        graph.SetMarkerStyle(
                            20 + point_index
                        )  # Different marker style for each file

                        point_index += 1

                        # Add entry to the legend with the file name (without .root extension)
                        legend.AddEntry(
                            graph,
                            os.path.basename(info["file"]).replace(".root", ""),
                            "lep",
                        )

                    # Now, draw the TGraph with the mean points on the canvas
                    graph.Draw("AP")  # "A" to draw axis, "P" to draw the points
                    graph.GetXaxis().SetTitle("Number of vertices reconstructed")
                    graph.GetYaxis().SetTitle("CPU Time used (s)")

                    # Draw the legend with file names corresponding to points
                    legend.Draw()

                    # Apply the CMS style labels to the canvas
                    add_cms_labels(canvas)

                    # Extract the run name from one of the file paths
                    run_name = extract_run_name(hist_info_list[0]["file"])

                    # Save the canvas with a filename indicating the run name
                    save_canvas(canvas, f"{hist_name}_{run_name}_mean.png")

                    # Clear the canvas for the next plot
                    canvas.Clear()

            else:
                print(f"Processing 1D histogram: {hist_name}")
                hist.Draw(draw_option)  # Default draw for 1D histograms

            # Draw the histogram on the canvas
            hist.Draw(draw_option)

            # Add entry to the legend
            legend.AddEntry(
                hist, os.path.basename(info["file"]).replace(".root", ""), "lep"
            )

        # Draw the legend
        legend.Draw()

        # Apply the CMS style labels to the canvas
        add_cms_labels(canvas)

        # Extract the run name from one of the file paths
        run_name = extract_run_name(hist_info_list[0]["file"])

        # Save the canvas with a filename indicating the run name
        save_canvas(canvas, f"{hist_name}_{run_name}.png")

        # Clear the canvas for the next plot
        canvas.Clear()

# Create a TGraph to plot the mean_y values against the overlap values
graph = ROOT.TGraph()

# Extract the overlap and mean_y values from hist_stats and populate the TGraph
point_index = 0  # To index points in the graph

for stat in PUBlockBordersvsPurityprofile_stats:
    # Extract the overlap number (convert overlap string like 'overlap_90' to integer)
    overlap_value = int(
        stat["overlap"].split("_")[-1]
    )  # Get the number after 'overlap_'
    mean_y = stat["mean_y"]  # Extract the mean_y value

    # Set point in the graph for each file
    graph.SetPoint(point_index, overlap_value, mean_y)

    point_index += 1

# Create a new canvas for the plot
canvas = ROOT.TCanvas("canvas_mean_y_vs_overlap", "", 800, 600)

# Draw the TGraph
graph.Draw("APL")  # "A" for axes, "P" for points, "L" for line

# Set graph title and axis labels
graph.SetTitle("Mean Purity vs Overlap")
graph.GetXaxis().SetTitle("Overlap (%)")
graph.GetYaxis().SetTitle("Mean Purity (%)")

# Optionally, you can style the graph, like marker type and line color
graph.SetMarkerStyle(20)  # Marker style (full circle)
graph.SetMarkerColor(ROOT.kBlue)  # Marker color
graph.SetLineColor(ROOT.kBlue)  # Line color

# Add CMS style labels to the canvas
add_cms_labels(canvas)

# Save the canvas as an image
save_canvas(canvas, "MeanPurity_vs_Overlap.png")
