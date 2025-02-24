import ROOT

def fit_histogram():
    # Open the ROOT file
    file = ROOT.TFile.Open("CutVarD0_2024_triggered.root", "READ")
    if not file or file.IsZombie():
        print("Error: Could not open file.")
        return
    
    # Get the histogram
    hist = file.Get("hRawFracPrompt")
    if not hist:
        print("Error: Histogram 'hRawFracPrompt' not found in file.")
        file.Close()
        return
    
    # Define the fitting function (pol0: constant function)
    # Customize histogram appearance
    hist.SetTitle("Raw Fraction of Prompt D^{0}")  # Set title
    hist.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")  # X-axis label
    hist.GetYaxis().SetTitle("Prompt fraction")  # Y-axis label
    hist.GetXaxis().SetRangeUser(1, 16)  # Set X-axis range
    hist.GetYaxis().SetRangeUser(0.9, 1.05)   # Set Y-axis range
    
    # Define the fitting function (pol0: constant function)
    fit_func = ROOT.TF1("fit_func", "pol0", 1.0, 16.0)  # Fit within specified range
    fit_func.SetLineColor(ROOT.kBlue)  # Set fit line color to red
    fit_func.SetLineWidth(2)  # Increase line width
    fit_func.SetLineStyle(1)  # Dashed line

    # Fit the histogram
    fit_result = hist.Fit(fit_func, "S")  # "S" option returns a TFitResultPtr
    
    # Print the fit results
    if fit_result:
        fit_result.Print("V")  # Print verbose fit result details
    
    # Draw the histogram and the fit
    canvas = ROOT.TCanvas("canvas", "Fit Result", 800, 600)
    hist.Draw("E")
    fit_func.Draw("same")

    # Add fit results as text on the plot
    fit_param = fit_func.GetParameter(0)  # Get the parameter of pol0
    fit_error = fit_func.GetParError(0)  # Get the error of the parameter
    text = ROOT.TLatex()
    text.SetNDC(True)  # Use normalized coordinates (0-1)
    text.SetTextSize(0.04)
    text.DrawLatex(0.15, 0.85, f"Fit result: {fit_param:.4f} #pm {fit_error:.4f}")

    # Save the plot as an image
    img_filename = "raw_prompt_frac_pol0.png"
    root_filename = "raw_prompt_frac_pol0.root"
    canvas.SaveAs(img_filename)
    canvas.SaveAs(root_filename)

    # Close the file
    file.Close()

# Run the function in the notebook
test = fit_histogram()
