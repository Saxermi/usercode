import os
import ROOT

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
    "AnotherHistogram1",
    "AnotherHistogram2",
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
                elif hist_name == "AnotherHistogram1":
                    print(f"Applying specific logic for {hist_name}")
                    # Add specific settings for AnotherHistogram1 if needed
                elif hist_name == "AnotherHistogram2":
                    print(f"Applying specific logic for {hist_name}")
                    # Add specific settings for AnotherHistogram2 if needed
            elif info["type"] == "TH2":
                print(f"Processing 2D histogram: {hist_name}")
                hist.Draw("COLZ")  # Draw 2D histograms with color palette
            else:
                print(f"Processing 1D histogram: {hist_name}")
                hist.Draw(draw_option)  # Default draw for 1D histograms

            # Draw the histogram on the canvas
            hist.Draw(draw_option)

            # Add entry to the legend
            # Add entry to the legend without the .root extension
            legend.AddEntry(
                hist, os.path.basename(info["file"]).replace(".root", ""), "lep"
            )

        # Draw the legend
        legend.Draw()

        # Apply the CMS style labels to the canvas
        add_cms_labels(canvas)

        # Save the canvas as a PNG file
        save_canvas(canvas, f"{hist_name}.png")

        # Clear the canvas for the next plot
        canvas.Clear()
