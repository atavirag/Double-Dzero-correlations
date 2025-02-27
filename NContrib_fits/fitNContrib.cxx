/// Quick macro to calculate the luminosity of the triggered dataset

#include <iostream>
#include <functional>

#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TKey.h"


using std::cout;
using std::cerr;
using std::endl;

bool verbose = true;

void fitNContrib() {

    //gSystem->Load("BifurcatedCB_cxx.so");  // Load the compiled library

    TFile *file = TFile::Open("~/MyMacros/Correlations_github/AnalysisResults_newModels_2024af.root", "read");
    TFile *fileMc = TFile::Open("~/MyMacros/Correlations_github/AnalysisResults_k3_PVContrib.root", "read"); // Using old models + short PV range
    TFile *fout = TFile::Open("fitNContrib.root", "RECREATE");

    if (!file) {
        cerr << ">> ERROR File not well readout" << endl;
        return;
    }

    if (!fileMc) {
        cerr << ">> ERROR MC File not well readout" << endl;
        return;
    }

    // Retrieve histogram
    TDirectory *dir = dynamic_cast<TDirectory*>(file->Get("hf-correlator-d-meson-pairs"));
    TDirectory *dirMc = dynamic_cast<TDirectory*>(fileMc->Get("hf-correlator-d-meson-pairs"));
    if (!dir || !dirMc) {
        cerr << ">> ERROR: Directory 'hf-correlator-d-meson-pairs' not found" << endl;
        return;
    }


    TH1F *hNContrib = dynamic_cast<TH1F*>(dir->Get("hPVContrib"));
    TH1F *hNContribMc = dynamic_cast<TH1F*>(dirMc->Get("hNContribMcGen"));
    if (!hNContrib) {
        cerr << ">> ERROR: Histogram 'hPVContrib' not found" << endl;
        return;
    }
    if (!hNContribMc) {
        cerr << ">> ERROR: Histogram 'hNContribMcGen' not found" << endl;
        return;
    }

    hNContrib->Rebin(1); // Rebin by a factor of 5 (merges 5 bins into 1)
    hNContribMc->Rebin(1); // Rebin by a factor of 5 (merges 5 bins into 1)
    hNContrib->GetXaxis()->SetRangeUser(1.0, 180.0);
    hNContribMc->GetXaxis()->SetRangeUser(1.0, 180.0);

    int nBinsX = hNContrib->GetNbinsX();
    std::cout << "Number of bins in X: " << nBinsX << std::endl;

    // Define the fitting function
    double mean = hNContrib->GetMean();
    double rms = hNContrib->GetRMS();

    double meanMc = hNContribMc->GetMean();
    double rmsMc = hNContribMc->GetRMS();

    // Define the variable corresponding to the x-axis (the PV contributors)
    RooRealVar pvContrib("pvContrib", "PV Contributors", 1.0, 180.0);
    RooRealVar mu("mu", "Mean", mean-rms, mean - 4*rms, mean + rms);
    RooRealVar sigma("sigma", "Sigma", rms, 0.0, 4 * rms);

    RooRealVar muMc("muMc", "Mean", meanMc-rmsMc, meanMc - 2*rmsMc, meanMc + rmsMc);
    RooRealVar sigmaMc("sigmaMc", "Sigma", rmsMc, 0.0, 2 * rmsMc);
    RooGaussian gauss("gauss", "Gaussian Core", pvContrib, muMc, sigmaMc);

    RooRealVar lambda("lambda", "Exponential decay", 0.05, -2, 2);
    RooExponential expo("expo", "Exponential Tail", pvContrib, lambda);

    RooRealVar fraction("fraction", "Fraction of Gaussian", 0.8, 0.0, 1.0);
    RooAddPdf modelMc("model", "Gaussian + Exponential", RooArgList(gauss, expo), fraction);

    RooRealVar tail("tail", "Tail parameter", -0.5, -2, 2);
    RooRealVar tailMc("tailMc", "Tail parameter", -0.5, -2, 2);

    RooNovosibirsk novo("novo", "Novosibirsk Function", pvContrib, mu, sigma, tail);
    RooNovosibirsk novoMc("novoMc", "Novosibirsk Function", pvContrib, muMc, sigmaMc, tailMc);


    RooDataHist data("data", "Histogram Data", RooArgList(pvContrib), hNContrib);
    RooDataHist dataMc("dataMc", "Histogram Data MC", RooArgList(pvContrib), hNContribMc);
    // Fit to data
    RooFitResult* fitResult = novo.fitTo(data, RooFit::Extended(true), RooFit::Save());
    //RooFitResult* fitResult = bifCB->fitTo(data, RooFit::Extended(true), RooFit::Save());
    //RooFitResult* fitResult = doubleGauss.fitTo(data, RooFit::Save());
    fitResult->Print("v");

    RooFitResult* fitResultMc = modelMc.fitTo(dataMc, RooFit::Extended(true), RooFit::Save());
    fitResultMc->Print("v");

    RooPlot* frame = pvContrib.frame();
    data.plotOn(frame);
    novo.plotOn(frame);

    RooPlot* frameMc = pvContrib.frame();
    dataMc.plotOn(frameMc);
    novoMc.plotOn(frameMc);

    // Create canvas and plot
    TCanvas *c = new TCanvas("c", "NContrib Fit", 800, 600);
    frame->Draw();

    // Save results
    fout->cd();
    frame->Write();
    c->SaveAs("hNContrib_fit_Novosibirsk.png");

    TCanvas *cMc = new TCanvas("cMc", "NContrib Fit MC", 800, 600);
    frameMc->Draw();

    // Save results
    fout->cd();
    frameMc->Write();
    cMc->SaveAs("hNContrib_MC_fit_Novosibirsk.png");

    // Compute residuals
    TCanvas *cRes = new TCanvas("cRes", "NContrib Residuals", 800, 600);
    RooHist* residuals = frame->residHist();
    RooPlot* residFrame = pvContrib.frame();
    residFrame->addPlotable(residuals, "P");
    residFrame->SetTitle("Fit Residuals for PV Contributors");  // Change the title
    residFrame->Draw();
    residFrame->Write();
    cRes->SaveAs("hNContrib_residuals.png");


    // Finally, divide the two functions
    //RooFormulaVar weightFunc("weightFunc", "@0/@1", RooArgList(novo, novoMc));
    RooRealVar normData("normData", "Normalization Data", novo.createIntegral(pvContrib)->getVal());
    RooRealVar normMc("normMc", "Normalization MC", novoMc.createIntegral(pvContrib)->getVal());

    std::cout << "Normalization Data: " << normData.getVal() << std::endl;
    std::cout << "Normalization MC: " << normMc.getVal() << std::endl;


    RooFormulaVar weightFunc("weightFunc", "(@0/@1)*(@3/@2)", 
                             RooArgList(novo, novoMc, normData, normMc));

    // Create a frame for plotting
    RooPlot* frameWeight = pvContrib.frame(RooFit::Title("Weight Function"));

    // Sample the function and plot it
    weightFunc.plotOn(frameWeight, RooFit::LineColor(kRed), RooFit::LineWidth(2));
    frameWeight->SetMinimum(0.001);

    // Retrieve the RooCurve from the frame
    RooCurve* curve = (RooCurve*)frameWeight->getObject(0); // The first object is usually the curve

    // Draw the frame
    TCanvas c1("c1", "Weight Function", 800, 600);
    c1.SetLogy();
    frameWeight->Draw();
    c1.SaveAs("weightFunction_rooPlot.png");  // Save the plot if needed

    TCanvas c2("c2", "Weight Function", 800, 600);
    c2.SetLogy();
    hNContrib->Scale(1.0 / hNContrib->Integral());
    hNContribMc->Scale(1.0 / hNContribMc->Integral());

    TH1F* hWeight = (TH1F*)hNContrib->Clone("hWeight");
    hWeight->Divide(hNContribMc);  // Compute bin-by-bin ratio

    hWeight->Draw();
    c2.SaveAs("weightFunction_histos.png");  // Save the plot if needed
    hWeight->Write();


    // Cleanup
    fout->Close();
    file->Close();
    delete c;
    delete cMc;
    delete cRes;
}


    /* RooRealVar alphaL("alphaL", "Left Alpha", 1.5, 0.1, 10.0);
    RooRealVar nL("nL", "Left n", 5.0, 0.1, 30.0);
    RooRealVar alphaR("alphaR", "Right Alpha", 2, 0.1, 15.0);
    RooRealVar nR("nR", "Right n", 10.0, 0.1, 30.0);

    nL.setVal(8.0);  // Fix nL
    nL.setConstant(kTRUE);
    nR.setVal(8.0);  // Fix nR
    nR.setConstant(kTRUE); */

    // Bifurcated CB (not good enough) -- If needed, load BifurcatedCB class
    //RooAbsPdf *bifCB = new BifurcatedCB("bifCB", "Bifurcated CB", pvContrib, mu, sigma, alphaL, nL, alphaR, nR);

    // Double gaussian (not good enough)
    /* RooRealVar frac("frac", "Fraction", 0.5, 0.0, 1.0);
    RooRealVar sigma1("sigma1", "Sigma", rms, 0.0, 2 * rms);
    RooRealVar sigma2("sigma2", "Sigma", rms, 0.0, 4 * rms);
    RooGaussian gauss1("gauss1", "Core Gaussian", pvContrib, mu, sigma1);
    RooGaussian gauss2("gauss2", "Wide Gaussian", pvContrib, mu, sigma2);
    RooAddPdf doubleGauss("doubleGauss", "Double Gaussian", RooArgList(gauss1, gauss2), RooArgList(frac)); */