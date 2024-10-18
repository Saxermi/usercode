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

# Create a canvas
print("Creating a canvas for drawing")
canvas = ROOT.TCanvas("canvas", "", 800, 700)

# List to store TProfile histograms
histograms = []

# Variable to track if the first histogram is drawn
first_histogram = True

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

    # Navigate to the 'offlinePrimaryVertices/efficiency' directory
    print(f"Checking for directory 'offlinePrimaryVertices/efficiency' in {file_path}")
    efficiency_dir = root_file.Get("offlinePrimaryVertices/efficiency")

    # Check if the directory exists
    if not efficiency_dir:
        print(
            f"Error: 'offlinePrimaryVertices/efficiency' directory not found in {file_path}!"
        )
        continue  # Skip this file and move to the next

    # Retrieve the TProfile histogram
    print(f"Checking for histogram 'PUBlockBordersvsPurityprofile' in {file_path}")
    hist = efficiency_dir.Get("PUBlockBordersvsPurityprofile")
    # Check if the object is a TProfile
    if isinstance(hist, ROOT.TProfile):
        print(f"'{file_path}' contains a TProfile histogram.")
    else:
        print(f"'{file_path}' does not contain a TProfile histogram.")

    # Check if the histogram was retrieved successfully
    if hist:
        print(f"Successfully retrieved histogram from {file_path}")

        # Prevent the histogram from being deleted when the file is closed
        hist.SetDirectory(0)  # Add this line to detach the histogram from the file

        # Set the color for the markers and error bars
        hist.SetMarkerColor(colors[i])
        hist.SetLineStyle(0)  # Hide the line, only show markers and error bars
        hist.SetLineColor(colors[i])  # Ensure error bars have the correct color

        # Set X and Y axis labels with appropriate size
        hist.GetXaxis().SetTitle("Distance to nearest block (mm)")
        hist.GetYaxis().SetTitle("Purity (%)")
        hist.GetXaxis().SetTitleSize(0.04)
        hist.GetYaxis().SetTitleSize(0.04)

        # Add the histogram to the list for future use
        histograms.append(hist)

        # Draw the TProfile histogram with error bars only
        draw_option = "E" if first_histogram else "E SAME"
        if first_histogram:
            print(
                f"Drawing first TProfile histogram (error bars only) from {file_path}"
            )
            first_histogram = False
        else:
            print(
                f"Drawing additional TProfile histogram (error bars only) from {file_path}"
            )
        hist.Draw(draw_option)  # Draw on the canvas
    else:
        print(
            f"Error: TProfile histogram 'PUBlockBordersvsPurityprofile' not found in {file_path}"
        )
        continue  # Skip this file and move to the next

    # Close the ROOT file to free resources (optional)
    root_file.Close()

# Add a legend
print("Adding a legend")
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
for i, hist in enumerate(histograms):
    if hist:  # Ensure hist is not None
        overlap_value = 10 * (i + 1)
        print(f"Adding entry for Overlap {overlap_value} to legend")
        legend.AddEntry(
            hist, f"Overlap {overlap_value}", "lep"
        )  # "lep" for marker and error bars
    else:
        print(f"Skipping entry for Overlap {10 * (i+1)} as the histogram is None")

legend.Draw()

# Add labels for CMS style
print("Adding CMS simulation labels")
label = ROOT.TLatex()
label.SetNDC(True)  # Set to normalized device coordinates
label.SetTextSize(0.040)
label.DrawLatex(0.100, 0.920, "#bf{CMS Simulation}")
label.DrawLatex(0.550, 0.920, "#sqrt{s} = 14 TeV")

# Update the canvas to show everything
print("Updating canvas")
canvas.Update()

# Save the plot as a .png file
output_file = "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/18.10/experimental_run_1/combined_PUBlockBordesvsPurity_overlap_histograms_experimental.png"
print(f"Saving canvas as {output_file}")
canvas.SaveAs(output_file)

# Create an application to keep the canvas open
print("Starting ROOT application to keep canvas open")
app = ROOT.TApplication("app", 0, [])
canvas.Draw()

# Run the application and wait for user to close it
app.Run()
