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
    TFile *fYield_LS = new TFile("~/MyMacros/corelations_with_phi/correlations_LS_oldBDTs.root", "READ");
    TFile *fYield_OS = new TFile("~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs.root", "READ");
    if (!fYield_OS || fYield_OS->IsZombie()) {
        std::cerr << "Error: Could not open ROOT file with yields OS!" << std::endl;
        return;
    }

    TFile *fileEffInt = TFile::Open("/home/andrea/MyMacros/Common/AccEffPreselD0ToKPi_no_ambiguous_bothBDTs_integrated.root", "read");
    TH1F *hEffInt = (TH1F *)fileEffInt->Get("hAccEffPreselD0ToKPi_no_ambiguous_bothBDTs_integratedAll");
    double integratedEff = hEffInt->GetBinContent(1);
    double integratedEffErr = hEffInt->GetBinError(1);

    TH1F *hYieldCorr_LS = (TH1F *)fYield_LS->Get("hYieldsCorr");
    TH1F *hYieldCorr_OS = (TH1F *)fYield_OS->Get("hYieldsCorr");

    TH1F *hYield_LS = (TH1F *)fYield_LS->Get("hYields");
    TH1F *hYield_OS = (TH1F *)fYield_OS->Get("hYields");

    if (!hYield_OS) {
        std::cerr << "Error: Histogram hYield_OS not found in file!" << std::endl;
        fYield_OS->Close();
        return;
    }

    double nCorr_LS = hYieldCorr_LS->GetBinContent(1);
    double nCorr_OS = hYieldCorr_OS->GetBinContent(1);
    double dNcorr_LS = hYieldCorr_LS->GetBinError(1);
    double dNcorr_OS = hYieldCorr_OS->GetBinError(1);

    double n_LS = hYield_LS->GetBinContent(1);
    double n_OS = hYield_OS->GetBinContent(1);
    double dN_LS = hYield_LS->GetBinError(1);
    double dN_OS = hYield_OS->GetBinError(1);

    double dCorrLumi = 5.83;

    std::cout << "N_corr LS: " << nCorr_LS << "+- " << hYieldCorr_LS->GetBinError(1) << std::endl;
    std::cout << "N_corr OS: " << nCorr_OS << "+- " <<  hYieldCorr_OS->GetBinError(1) << std::endl;
    std::cout << "N_corr OS / N_corr LS: " << "+- " <<  nCorr_OS /nCorr_LS  << std::endl;

    std::cout << "N LS: " << n_LS << "+- " << hYield_LS->GetBinError(1) << std::endl;
    std::cout << "N OS: " << n_OS << "+- " <<  hYield_OS->GetBinError(1) << std::endl;
    std::cout << "N OS / N LS: " << n_OS /n_LS  << std::endl;

    double relUnc_nCorr_LS = dNcorr_LS / nCorr_LS;
    double relUnc_nCorr_OS = dNcorr_OS / nCorr_OS;
    double relUnc_n_LS = dN_LS / n_LS;
    double relUnc_n_OS = dN_OS / n_OS;

    double relUnc_corrLumi = dCorrLumi / corrLumi;
    double relUnc_promptFraction = 0.0003 / promptFraction;
    double relUnc_BR = 0.0003 / BR;

    double relUnc_xsec_LS = sqrt(pow(relUnc_nCorr_LS, 2) +
                              4 * pow(relUnc_promptFraction, 2) +
                              4 * pow(relUnc_BR, 2) +
                              pow(relUnc_corrLumi, 2));

    double relUnc_xsec_OS = sqrt(pow(relUnc_nCorr_OS, 2) +
                              4 * pow(relUnc_promptFraction, 2) +
                              4 * pow(relUnc_BR, 2) +
                              pow(relUnc_corrLumi, 2));

    // Calculate the cross sections
    double xsec_LS =  (nCorr_LS * promptFraction * promptFraction) / (BR * BR * corrLumi);
    double dXsec_LS = xsec_LS * relUnc_xsec_LS;
    std::cout << "        Like-sign cross section: " << xsec_LS << "+- " << dXsec_LS << endl;

    double xsec_OS =  (nCorr_OS * promptFraction * promptFraction) / (BR * BR * corrLumi);
    double dXsec_OS = xsec_OS * relUnc_xsec_OS;
    std::cout << "        Opposite-sign cross section: " << xsec_OS << "+- " << dXsec_OS << endl;

    double relUnc_Eff = integratedEffErr / integratedEff;
    double relUnc_xsec_int_LS = sqrt(pow(relUnc_n_LS, 2) +
                              4 * pow(relUnc_promptFraction, 2) +
                              4 * pow(relUnc_BR, 2) +
                              4 * pow(relUnc_Eff, 2) +
                              pow(relUnc_corrLumi, 2));

    double relUnc_xsec_int_OS = sqrt(pow(relUnc_n_OS, 2) +
                              4 * pow(relUnc_promptFraction, 2) +
                              4 * pow(relUnc_BR, 2) +
                              4 * pow(relUnc_Eff, 2) +
                              pow(relUnc_corrLumi, 2));

    // Calculate the cross sections with the integrated efficiency
    double xsec_int_LS =  (n_LS * promptFraction * promptFraction) / (integratedEff * integratedEff * BR * BR * corrLumi);
    double dXsec_int_LS = xsec_int_LS * relUnc_xsec_int_LS;
    std::cout << "        Like-sign cross section with integrated efficiency: " << xsec_int_LS << "+- " << dXsec_int_LS << endl;

    double xsec_int_OS =  (n_OS * promptFraction * promptFraction) / (integratedEff * integratedEff * BR * BR * corrLumi);
    double dXsec_int_OS = xsec_int_OS * relUnc_xsec_int_OS;
    std::cout << "        Opposite-sign cross section with intetgrated efficiency: " << xsec_int_OS << "+- " << dXsec_int_OS << endl;
}