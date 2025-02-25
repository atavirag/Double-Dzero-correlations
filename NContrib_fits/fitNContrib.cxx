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

    hNContrib->Rebin(5); // Rebin by a factor of 5 (merges 5 bins into 1)
    hNContribMc->Rebin(5); // Rebin by a factor of 5 (merges 5 bins into 1)
    hNContrib->GetXaxis()->SetRangeUser(10.0, 180.0);
    hNContribMc->GetXaxis()->SetRangeUser(10.0, 120.0);

    // Define the fitting function
    double mean = hNContrib->GetMean();
    double rms = hNContrib->GetRMS();

    double meanMc = hNContribMc->GetMean();
    double rmsMc = hNContribMc->GetRMS();

    // Define the variable corresponding to the x-axis (the PV contributors)
    RooRealVar pvContrib("pvContrib", "PV Contributors", 10.0, 180.0);
    RooRealVar mu("mu", "Mean", mean-rms, mean - 4*rms, mean + rms);
    RooRealVar sigma("sigma", "Sigma", rms, 0.0, 4 * rms);

    RooRealVar pvContribMc("pvContribMc", "PV MC Contributors", 10.0, 120.0); // not needed when I get the right PV MC range
    RooRealVar muMc("muMc", "Mean", meanMc-rmsMc, meanMc - 2*rmsMc, meanMc + rmsMc);
    RooRealVar sigmaMc("sigmaMc", "Sigma", rmsMc, 0.0, 2 * rmsMc);

    RooRealVar tail("tail", "Tail parameter", -0.5, -2, 2);
    RooRealVar tailMc("tailMc", "Tail parameter", -0.5, -2, 2);

    RooNovosibirsk novo("novo", "Novosibirsk Function", pvContrib, mu, sigma, tail);
    RooNovosibirsk novoMc("novoMc", "Novosibirsk Function", pvContribMc, muMc, sigmaMc, tailMc);


    RooDataHist data("data", "Histogram Data", RooArgList(pvContrib), hNContrib);
    RooDataHist dataMc("dataMc", "Histogram Data MC", RooArgList(pvContribMc), hNContribMc);
    // Fit to data
    RooFitResult* fitResult = novo.fitTo(data, RooFit::Extended(true), RooFit::Save());
    //RooFitResult* fitResult = bifCB->fitTo(data, RooFit::Extended(true), RooFit::Save());
    //RooFitResult* fitResult = doubleGauss.fitTo(data, RooFit::Save());
    fitResult->Print("v");

    RooFitResult* fitResultMc = novoMc.fitTo(dataMc, RooFit::Extended(true), RooFit::Save());
    fitResultMc->Print("v");

    RooPlot* frame = pvContrib.frame();
    data.plotOn(frame);
    novo.plotOn(frame);

    RooPlot* frameMc = pvContribMc.frame();
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