import ROOT


# Open the ROOT file located in the working directory
file = ROOT.TFile.Open(
    "/home/sam/ownCloud - Michael Saxer (zhaw.ch)@drive.switch.ch/PA/Root_archiv/12.10/pv.root"
)


# Function to recursively list keys in directories
def list_keys(directory, indent=""):
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        print(f"{indent}{key.GetName()} ({obj.ClassName()})")
        # If the object is a directory, call the function recursively
        if obj.InheritsFrom("TDirectory"):
            list_keys(obj, indent + "  ")


# List all keys starting from the root directory
list_keys(file)
# Navigate to the correct directory
efficiency_dir = file.Get("offlinePrimaryVertices/efficiency")
if not efficiency_dir:
    print("Error: 'offlinePrimaryVertices/efficiency' directory not found!")
    exit()

# Retrieve the histogram by name
hist = efficiency_dir.Get("PURecoVsSimZPositionBlock")
if not hist:
    print("Error: Histogram 'PURecoVsSimZPositionBlock' not found!")
    exit()

# Create a canvas to draw the histogram
canvas = ROOT.TCanvas(
    "canvas", "Distribution of Z-axis Positions with Shaded Border Area", 800, 600
)

# Draw the 2D histogram
hist.Draw("COLZ")

# Define and draw the borders
upper_border = ROOT.TLine(1.0, 4.5, 5.0, 4.5)  # Upper border (Theoretical, y = 4.5)
upper_border.SetLineColor(ROOT.kRed)
upper_border.SetLineStyle(2)  # Dashed
upper_border.Draw()

lower_border = ROOT.TLine(1.0, 1.5, 5.0, 1.5)  # Lower border (Theoretical, y = 1.5)
lower_border.SetLineColor(ROOT.kGreen)
lower_border.SetLineStyle(2)  # Dashed
lower_border.Draw()

upper_border_measured = ROOT.TLine(
    4.5, 1.0, 4.5, 5.0
)  # Upper border (Measured, x = 4.5)
upper_border_measured.SetLineColor(ROOT.kMagenta)
upper_border_measured.SetLineStyle(2)  # Dashed
upper_border_measured.Draw()

lower_border_measured = ROOT.TLine(
    2.0, 1.0, 2.0, 5.0
)  # Lower border (Measured, x = 2.0)
lower_border_measured.SetLineColor(ROOT.kOrange)
lower_border_measured.SetLineStyle(2)  # Dashed
lower_border_measured.Draw()

# Define and draw the shaded area
shaded_area = ROOT.TBox(2.0, 1.5, 4.5, 4.5)  # Box representing the shaded area
# shaded_area.SetFillColorAlpha(ROOT.kRed, 0.005)  # Semi-transparent light red color
shaded_area.SetLineColor(ROOT.kRed)  # Set the line color for the outline
shaded_area.SetLineStyle(2)  # Optional: Dashed or different style
shaded_area.SetFillStyle(0)
shaded_area.Draw("same")

# Create and draw the legend
legend = ROOT.TLegend(0.1, 0.7, 0.3, 0.9)
legend.AddEntry(hist, "Data Points", "p")
legend.AddEntry(upper_border, "Upper Border (Theoretical)", "l")
legend.AddEntry(lower_border, "Lower Border (Theoretical)", "l")
legend.AddEntry(upper_border_measured, "Upper Border (Measured)", "l")
legend.AddEntry(lower_border_measured, "Lower Border (Measured)", "l")
legend.AddEntry(shaded_area, "Overlap Area", "f")
legend.Draw()

# Save the canvas to a file
canvas.SaveAs("PURecoVsSimZPositionBlock_with_borders.png")

# Keep the canvas open for review
input("Press Enter to exit...")
