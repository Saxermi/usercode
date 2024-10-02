import ROOT

# Create a new ROOT file
f = ROOT.TFile("testfile.root", "RECREATE")

# Create a histogram
h = ROOT.TH1F("h", "test histogram", 100, -4, 4)

# Fill the histogram with random Gaussian-distributed data
h.FillRandom("gaus")

# Write the histogram to the ROOT file
h.Write()

# Create a canvas to draw the histogram
c = ROOT.TCanvas("c", "Canvas for histogram", 800, 600)

# Draw the histogram on the canvas
h.Draw()

# Save the canvas as a PNG file
c.SaveAs("histogram.png")

# Close the ROOT file
f.Close()
