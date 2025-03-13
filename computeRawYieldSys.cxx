#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <string>

void computeRawYieldSys() {
    TFile* fmain = TFile::Open("~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_main.root", "READ");
    if (!fmain || fmain->IsZombie()) {
        std::cerr << "Error: Could not open fmain" << std::endl;
        return;
    }

    TH1D* hYields_main = dynamic_cast<TH1D*>(fmain->Get("hYields"));
    if (!hYields_main) {
        std::cerr << "Error: Histogram 'hYields' not found in fmain" << std::endl;
        fmain->Close();
        return;
    }
    // Get main yields
    double rawYield_main = hYields_main->GetBinContent(1);
    double rawYieldErr_main = hYields_main->GetBinError(1);

    std::vector<const char*> sysLabels = {"Exp(pol2) Range 1.74-2.05", "Exp(pol2) Range 1.73-2.04",
                                          "pol1 Range 1.73-2.04",  "pol2 Range 1.73-2.04", "pol2 Range 1.74-2.05",
                                          "expo Range 1.74-2.05", "expo Range 1.73-2.04"};
    const std::vector<std::string>& filenames = {"~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_range205.root",
                                                 "~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_range173.root",
                                                 "~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_range173_pol1.root",
                                                 "~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_range173_pol2.root",
                                                 "~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_range205_pol2.root",
                                                 "~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_range205_expo.root",
                                                 "~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_range173_expo.root"};
    std::vector<double> rawYields = {};
    std::vector<double> rawYieldsErr = {};
    for (const auto& filename : filenames) {
        TFile* file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open " << filename << std::endl;
            continue;
        }

        TH1D* hYields = dynamic_cast<TH1D*>(file->Get("hYields"));
        if (!hYields) {
            std::cerr << "Error: Histogram 'hYields' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        int binIndex = 1; // SgnSgn yield
        double binContent = hYields->GetBinContent(binIndex);
        rawYields.push_back(binContent);
        double binError = hYields->GetBinError(binIndex);
        rawYieldsErr.push_back(binError);

        std::cout << "File: " << filename << " | SgnSgn Content: " << binContent 
                  << " | SgnSgn Error: " << binError << std::endl;

        file->Close();
    }

    // Compute the ratio wrt the main value
    std::vector<double> ratioYields = {};
    std::vector<double> ratioYieldsErr = {};
    for (int i = 0; i < rawYields.size(); i++) {
        ratioYields.push_back(rawYields[i]/rawYield_main);
        ratioYieldsErr.push_back(sqrt(pow(rawYieldsErr[i]/rawYields[i],2)+pow(rawYieldErr_main/rawYield_main,2)));
    }

    // plot systematics
    TH1D *hRawYieldSys = new TH1D("hRawYieldSys", "Raw Yield Sys.", rawYields.size(), 0, rawYields.size());
    double squaredSum = 0.;
    for (int i = 0; i < rawYields.size(); i++) {
        if (sysLabels.size() != rawYields.size()) {
            std::cerr << "Mismatch between systematics labels and raw yields size!" << std::endl;
        }
        hRawYieldSys->GetXaxis()->SetBinLabel(i+1, sysLabels[i]);
        hRawYieldSys->SetBinContent(i+1, ratioYields[i]);
        hRawYieldSys->SetBinError(i+1, ratioYieldsErr[i]);
        double deviation = ratioYields[i] - 1;
        squaredSum += deviation*deviation;
    }
    double finalYieldSys = sqrt(squaredSum);
    std::cout << "Raw Yield systematics (quadrature sum): " << finalYieldSys << std::endl;

    hRawYieldSys->SetMarkerColor(kBlue);
    hRawYieldSys->SetMarkerStyle(21);
    hRawYieldSys->SetLineColor(kBlue);
    hRawYieldSys->GetYaxis()->SetRangeUser(0.8, 1.2);

    TLine *line = new TLine(0, 1, rawYields.size(), 1);

    TLegend *legend = new TLegend(0.15, 0.7, 0.3, 0.8);
    legend->AddEntry(hRawYieldSys, "trial/main result", "PE");
    TCanvas *c = new TCanvas("c", "Raw Yield Systematics", 800, 600);
    gStyle->SetOptStat(0);
    hRawYieldSys->Draw("PE");
    line->Draw("same");
    legend->Draw();
    c->SaveAs("rawYieldSys_OS.png");
    c->SaveAs("rawYieldSys_OS.root");
    c->Close();
    fmain->Close();
}

