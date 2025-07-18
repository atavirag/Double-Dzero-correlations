#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <string>

bool isStatLinear=true; 
void calculateSys(std::vector<std::string> filenames, std::vector<std::string> sysLabels, std::string outname, TFile* fout);
void GetStatistics(TH1D* histo, double &average, double &averageStat, double &averageSyst);

void computeRawYieldSys() {

    std::vector<std::string> sysLabels = {};
    std::vector<std::string> filenamesOS = {};
    std::vector<std::string> filenamesLS = {};

    std::vector<std::pair<double, double>> massRanges = {
        {1.74, 2.04}, {1.73, 2.03}, {1.75, 2.05}
    };
    std::vector<std::string> bkgFuncs = {"poly2", "expo", "expPoly2"};

    // Add poly2 (1.74, 2.05)
    std::string filenameLS_extra = "systematics/rawYield_sys/correlations_LS_poly2_1.74_2.05_effs.root";
    std::string filenameOS_extra = "systematics/rawYield_sys/correlations_OS_poly2_1.74_2.05_effs.root";

    std::string label_extra = "poly2 ( 1.74, 2.05)";

    filenamesOS.push_back(filenameOS_extra);
    filenamesLS.push_back(filenameLS_extra);

    sysLabels.push_back(label_extra);
    
    for (const auto& bkgFunc : bkgFuncs) {
        for (const auto& massLim : massRanges) {
            std::string massLimStrFirst = Form("%.2f", massLim.first);
            std::string massLimStrSecond = Form("%.2f", massLim.second);

            std::string filenameLS = "systematics/rawYield_sys/correlations_LS_" + bkgFunc + "_" + 
                                        massLimStrFirst + "_" + massLimStrSecond + "_effs.root";
            std::string filenameOS = "systematics/rawYield_sys/correlations_OS_" + bkgFunc + "_" + 
                                        massLimStrFirst + "_" + massLimStrSecond + "_effs.root";

            std::string label = bkgFunc + " (" + massLimStrFirst + ", " + massLimStrSecond + ")";

            if (bkgFunc.compare("poly2") == 0 && massLim.first == 1.74 && massLim.second == 2.04) {
                continue;
            }
            if (bkgFunc.compare("expPoly2") == 0 && massLim.first == 1.74 && massLim.second == 2.04) {
                continue;
            }
            filenamesOS.push_back(filenameOS);
            filenamesLS.push_back(filenameLS);

            cout << label << endl;
            sysLabels.push_back(label);
        }

    }
    // Add expPoly2 (1.73, 2.04)
    std::string filenameLS_extra2 = "systematics/rawYield_sys/correlations_LS_expPoly2_1.72_2.05_effs.root";
    std::string filenameOS_extra2 = "systematics/rawYield_sys/correlations_OS_expPoly2_1.72_2.05_effs.root";

    std::string label_extra2 = "expPoly2 ( 1.72, 2.05)";

    filenamesOS.push_back(filenameOS_extra2);
    filenamesLS.push_back(filenameLS_extra2);

    sysLabels.push_back(label_extra2);

    TFile* fout = TFile::Open("corrYieldSystematics.root", "RECREATE");
    calculateSys(filenamesOS, sysLabels, "CorrYieldSys_OS", fout);
    calculateSys(filenamesLS, sysLabels, "CorrYieldSys_LS", fout);
}

void calculateSys(std::vector<std::string> filenames, std::vector<std::string> sysLabels, std::string outname, TFile* fout) {

    double average = 0., averageStat = 0., averageSyst = 0.;

    std::vector<double> rawYields = {};
    std::vector<double> rawYieldsErr = {};
    std::vector<double> chi2_ndfs = {};
    std::vector<double> widths = {};
    std::vector<double> widthsErr = {};
    for (const auto& filename : filenames) {
        TFile* file = TFile::Open(filename.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Could not open " << filename << std::endl;
            continue;
        }

        TH1D* hYields = dynamic_cast<TH1D*>(file->Get("hYieldsCorr"));
        if (!hYields) {
            std::cerr << "Error: Histogram 'hYields' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        TH1F* hChi2_ndf = dynamic_cast<TH1F*>(file->Get("hChi2_ndf_Corr"));
        if (!hChi2_ndf) {
            std::cerr << "Error: Histogram 'hChi2_ndf_' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        TH1F* hWidth = dynamic_cast<TH1F*>(file->Get("hWidth"));
        if (!hWidth) {
            std::cerr << "Error: Histogram 'hWidth' not found in " << filename << std::endl;
            file->Close();
            continue;
        }

        int binIndex = 1; // SgnSgn yield
        double binContent = hYields->GetBinContent(binIndex);
        rawYields.push_back(binContent);
        double binError = hYields->GetBinError(binIndex);
        rawYieldsErr.push_back(binError);

        double chi2_ndf = hChi2_ndf->GetBinContent(binIndex);
        chi2_ndfs.push_back(chi2_ndf);

        double width = hWidth->GetBinContent(binIndex);
        widths.push_back(width);
        double widthError = hWidth->GetBinError(binIndex);
        widthsErr.push_back(widthError);

        cout << "binContent: " << binContent << ", chi2_ndf: " << chi2_ndf << ", width: " << width << endl;

        file->Close();
    }

    // plot systematics
    TH1D *hRawYieldSys = new TH1D("hRawYieldSys", "Raw Yield Systematics", rawYields.size(), 0, rawYields.size());
    TH1D *hWidthSys = new TH1D("hWidthSys", "Widths", rawYields.size(), 0, rawYields.size());
    TH1D *hChi2Sys = new TH1D("hChi2Sys", "#chi^{2}/ndf", rawYields.size(), 0, rawYields.size());
    double squaredSum = 0.;
    for (int i = 0; i < rawYields.size(); i++) {
        if (sysLabels.size() != rawYields.size()) {
            std::cerr << "Mismatch between systematics labels and raw yields size!" << std::endl;
        }
        hRawYieldSys->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());
        hWidthSys->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());
        hChi2Sys->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());
        if (rawYieldsErr[i] == 0) continue;
        hRawYieldSys->SetBinContent(i+1, rawYields[i]);
        hRawYieldSys->SetBinError(i+1, rawYieldsErr[i]);

        hWidthSys->SetBinContent(i+1, widths[i]);
        hWidthSys->SetBinError(i+1, widthsErr[i]);

        hChi2Sys->SetBinContent(i+1, chi2_ndfs[i]);
    }
    GetStatistics(hRawYieldSys, average, averageStat, averageSyst);
    std::cout << "Average: " << average << ", average stat: " << averageStat << ", average syst: " << averageSyst << endl;

    //hRawYieldSys->SetMarkerColor(kBlue);
    hRawYieldSys->SetMarkerStyle(21);
    //hRawYieldSys->SetLineColor(kBlue);
    hRawYieldSys->GetYaxis()->SetRangeUser(average/3, average*5/3);

    // Create TPaveText
    TPaveText *pave = new TPaveText(0.15, 0.7, 0.8, 0.85, "NDC");  // Normalized coordinates
    pave->SetFillColor(0);  // Transparent background
    pave->SetTextAlign(12); // Align left
    pave->SetTextSize(0.03);
    pave->AddText(Form("#it{N}_{DD}: %.0f #pm %.0f (stat., %.2f%%) #pm %.0f (syst., %.2f%%)", average, averageStat, (averageStat/average)*100, averageSyst, (averageSyst/average)*100));

    TLine *line = new TLine(0, average, rawYields.size(), average);
    line->SetLineWidth(2);
    line->SetLineColor(kRed);

    TLine *statsUpper = new TLine(0, average + averageStat, rawYields.size(), average + averageStat);
    statsUpper->SetLineWidth(2);
    statsUpper->SetLineColor(kGray+1);
    statsUpper->SetLineStyle(kDashed);

    TLine *statsLower = new TLine(0, average - averageStat, rawYields.size(), average - averageStat);
    statsLower->SetLineWidth(2);
    statsLower->SetLineColor(kGray+1);
    statsLower->SetLineStyle(kDashed);

    TCanvas *c = new TCanvas("c", "Raw Yield Systematics", 800, 600);

    gStyle->SetOptStat(0);

    TBox *colorBox = new TBox(0.01, average-averageSyst, rawYields.size(), average + averageSyst);  // (x1, y1, x2, y2)
    colorBox->SetFillColor(kRed-10);               // Set box color

    hRawYieldSys->Draw("PE");
    colorBox->Draw("same");
    line->Draw("same");
    statsLower->Draw("same");
    pave->Draw("same");
    statsUpper->Draw("same");
    hRawYieldSys->Draw("samePE");
    //legend->Draw();
    std::string outPng = outname + ".png";
    c->SaveAs(outPng.c_str());

    // Add plot with the chi2_ndf
    TCanvas *c_chi2 = new TCanvas("c_chi2", "#chi^{2}/ndf", 800, 400);

    hChi2Sys->SetMarkerStyle(21);
    hChi2Sys->SetMarkerSize(1.5);
    hChi2Sys->GetYaxis()->SetRangeUser(hChi2Sys->GetBinContent(1)-0.6, hChi2Sys->GetBinContent(1)+0.6);
    hChi2Sys->Draw("P");
    std::string outChiPng = "Chi2_over_ndf_" + outname + ".png";
    c_chi2->SaveAs(outChiPng.c_str());

    // Add plot with the widths
    TCanvas *c_width = new TCanvas("c_width", "Widths", 800, 400);
    hWidthSys->SetMarkerStyle(21);
    hWidthSys->SetMarkerSize(1.3);
    hWidthSys->GetYaxis()->SetRangeUser(hWidthSys->GetBinContent(1)-0.003, hWidthSys->GetBinContent(1)+0.003);
    //hWidthSys->SetLineSize(2);
    hWidthSys->SetLineColor(kBlack);

    hWidthSys->Draw("PE");
    std::string outWidthPng = "Width_" + outname + ".png";
    c_width->SaveAs(outWidthPng.c_str());

    fout->cd();
    c->Write();
    c_chi2->Write();
    c_width->Write();

    c->Close();
    c_chi2->Close();
    c_width->Close();

    delete hRawYieldSys;
    delete hWidthSys;
    delete hChi2Sys;
}

void GetStatistics(TH1D* histo, double &average, double &averageStat, double &averageSyst) {
    double av=0., avSt=0., rms=0.;
    int nbins = histo->GetNbinsX();
    int nFits=0;

    for (int i=1; i<nbins; i++){
        double value = histo->GetBinContent(i);
        double unc = histo->GetBinError(i);
        if(value>0){
            av += value;
            if(isStatLinear) {
                avSt += unc;
            } else{
                avSt += unc*unc;
            }
            nFits++;
        }
        cout <<" i="<<i<<" val="<<value;
    }
    cout<<endl;
    av /= nFits;
    if(isStatLinear) {
        avSt /= nFits;
    } else  avSt = TMath::Sqrt( avSt )/nFits;

    for (int i=1; i<nbins; i++){
        double value = histo->GetBinContent(i);
        if(value>0) rms += (value - av)*(value - av);
    }

    rms = TMath::Sqrt( rms /(nFits-1) );
    average=av;
    averageStat=avSt;
    averageSyst=rms;
}