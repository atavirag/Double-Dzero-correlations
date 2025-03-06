#include "TCanvas.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TH2.h"

void computeCrossSection() {
    // Calculate the cross sections
    // xsec = (N_corr * f_prompt^2) / (BR^2 * L_int)
    float const BR = 0.03950;
    float const BRErr = 0.0003;

    // Get lumi from computeLumiTriggered.cxx
    TFile *fLumi = new TFile("Integrated_Lumi.root", "READ");
    if (!fLumi || fLumi->IsZombie()) {
        std::cerr << "Error: Could not open ROOT file with integrated luminosity!" << std::endl;
        return;
    }

    // Retrieve the histogram
    TH1F *hCorrLumi = (TH1F *)fLumi->Get("hCorrLumi");
    if (!hCorrLumi) {
        std::cerr << "Error: Histogram hCorrLumi not found in file!" << std::endl;
        fLumi->Close();
        return;
    }

    // Get the integrated luminosity value
    double corrLumi = hCorrLumi->GetBinContent(1);
    std::cout << "Corrected Integrated Luminosity: " << corrLumi << " pb⁻¹" << std::endl;

    // Get prompt fraction value
    TFile *fPromptFrac = new TFile("cut_variation/raw_prompt_frac_pol0.root", "READ");
    if (!fPromptFrac || fPromptFrac->IsZombie()) {
        std::cerr << "Error: Could not open ROOT file with prompt fraction!" << std::endl;
        return;
    }
    // Retrieve the histogram
    TH1F *hPromptFrac = (TH1F *)fPromptFrac->Get("hFitParam");
    if (!hPromptFrac) {
        std::cerr << "Error: Histogram hPromptFrac not found in file!" << std::endl;
        fPromptFrac->Close();
        return;
    }

    // Get the prompt fraction value
    double promptFraction = hPromptFrac->GetBinContent(1);
    std::cout << "Prompt fraction: " << promptFraction << std::endl;

    // Get the corrected yields
    TFile *fYield_LS = new TFile("~/MyMacros/corelations_with_phi/correlations_LS_newBDTs.root", "READ");
    TFile *fYield_OS = new TFile("~/MyMacros/corelations_with_phi/correlations_OS_newBDTs.root", "READ");
    if (!fYield_OS || fYield_OS->IsZombie()) {
        std::cerr << "Error: Could not open ROOT file with yields OS!" << std::endl;
        return;
    }

    TH1F *hYield_LS = (TH1F *)fYield_LS->Get("hYieldsCorr");
    TH1F *hYield_OS = (TH1F *)fYield_OS->Get("hYieldsCorr");

    if (!hYield_OS) {
        std::cerr << "Error: Histogram hYield_OS not found in file!" << std::endl;
        fYield_OS->Close();
        return;
    }

    double nCorr_LS = hYield_LS->GetBinContent(1);
    double nCorr_OS = hYield_OS->GetBinContent(1);
    std::cout << "N_corr LS: " << nCorr_LS << std::endl;
    std::cout << "N_corr OS: " << nCorr_OS << std::endl;

    // Calculate the cross sections
    double xsec_LS =  (nCorr_LS * promptFraction * promptFraction) / (BR * BR * corrLumi);
    std::cout << "        Like-sign cross section: " << xsec_LS << endl;

    double xsec_OS =  (nCorr_OS * promptFraction * promptFraction) / (BR * BR * corrLumi);
    std::cout << "        Opposite-sign cross section: " << xsec_OS << endl;
}