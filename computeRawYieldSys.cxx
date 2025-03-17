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

bool isStatLinear=true; 
void calculateSys(std::vector<std::string> filenames, std::vector<std::string> sysLabels, std::string outname);
void GetStatistics(TH1D* histo, double &average, double &averageStat, double &averageSyst);

void computeRawYieldSys() {

    std::vector<std::string> sysLabels = {};
    std::vector<std::string> filenamesOS = {};
    std::vector<std::string> filenamesLS = {};

    std::vector<std::pair<double, double>> massRanges = {
        {1.74, 2.05}, {1.73, 2.05}, {1.72, 2.04}
    };
    std::vector<std::string> bkgFuncs = {"poly2", "expPoly1", "expo", "expPoly2"};
    
    for (const auto& bkgFunc : bkgFuncs) {
        for (const auto& massLim : massRanges) {
            std::string massLimStrFirst = Form("%.2f", massLim.first);
            std::string massLimStrSecond = Form("%.2f", massLim.second);

            std::string filenameLS = "rawYield_sys/correlations_LS_" + bkgFunc + "_" + 
                                        massLimStrFirst + "_" + massLimStrSecond + ".root";
            std::string filenameOS = "rawYield_sys/correlations_OS_" + bkgFunc + "_" + 
                                        massLimStrFirst + "_" + massLimStrSecond + ".root";

            std::string label = bkgFunc + " (" + massLimStrFirst + ", " + massLimStrSecond + ")";
            if (filenameOS.compare(0, filenameOS.size(), "rawYield_sys/correlations_OS_poly2_1.72_2.03.root") == 0) {
                continue;  // this one didn't fit properly
            }

            filenamesOS.push_back(filenameOS);
            filenamesLS.push_back(filenameLS);

            cout << label << endl;
            sysLabels.push_back(label);
        }
    }
    calculateSys(filenamesOS, sysLabels, "rawYieldSys_OS");
    calculateSys(filenamesLS, sysLabels, "rawYieldSys_LS");
}

void calculateSys(std::vector<std::string> filenames, std::vector<std::string> sysLabels, std::string outname) {
    TFile* fmain;
    if (outname.compare(outname.size() - 2, outname.size(), "OS") == 0) {
        fmain = TFile::Open("~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_main.root", "READ");
    } else {
        fmain = TFile::Open("~/MyMacros/corelations_with_phi/correlations_LS_oldBDTs_main.root", "READ");
    }
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
    double average = 0., averageStat = 0., averageSyst = 0.;
    double rawYield_main = hYields_main->GetBinContent(1);
    double rawYieldErr_main = hYields_main->GetBinError(1);

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

        //std::cout << "File: " << filename << " | SgnSgn Content: " << binContent 
        //          << " | SgnSgn Error: " << binError << std::endl;

        file->Close();
    }
    
    // Compute the ratio wrt the main value
    std::vector<double> ratioYields = {};
    std::vector<double> ratioYieldsErr = {};
    for (int i = 0; i < rawYields.size(); i++) {
        ratioYields.push_back(rawYields[i]/rawYield_main);
        cout << rawYields[i]/rawYield_main << endl;
        ratioYieldsErr.push_back(sqrt(pow(rawYieldsErr[i]/rawYields[i],2)+pow(rawYieldErr_main/rawYield_main,2)));
    }

    // plot systematics
    TH1D *hRawYieldSys = new TH1D("hRawYieldSys", "Raw Yield Sys.", rawYields.size(), 0, rawYields.size());
    double squaredSum = 0.;
    for (int i = 0; i < rawYields.size(); i++) {
        if (sysLabels.size() != rawYields.size()) {
            std::cerr << "Mismatch between systematics labels and raw yields size!" << std::endl;
        }
        hRawYieldSys->GetXaxis()->SetBinLabel(i+1, sysLabels[i].c_str());
        if (rawYieldsErr[i] == 0) continue;
        hRawYieldSys->SetBinContent(i+1, rawYields[i]);
        hRawYieldSys->SetBinError(i+1, rawYieldsErr[i]);
        GetStatistics(hRawYieldSys, average, averageStat, averageSyst);
        double deviation = ratioYields[i] - 1;
        squaredSum += deviation*deviation;
    }
    std::cout << "Average: " << average << ", average stat: " << averageStat << ", average syst: " << averageSyst << endl;

    hRawYieldSys->SetMarkerColor(kBlue);
    hRawYieldSys->SetMarkerStyle(21);
    hRawYieldSys->SetLineColor(kBlue);
    hRawYieldSys->GetYaxis()->SetRangeUser(average/2, average*1.5);

    TLine *line = new TLine(0, average, rawYields.size(), average);

    TLegend *legend = new TLegend(0.15, 0.7, 0.3, 0.8);
    legend->AddEntry(hRawYieldSys, "trial/main result", "PE");
    TCanvas *c = new TCanvas("c", "Raw Yield Systematics", 800, 600);
    gStyle->SetOptStat(0);
    hRawYieldSys->Draw("PE");
    line->Draw("same");
    legend->Draw();
    std::string outPng = outname + ".png";
    std::string outRoot = outname + ".root";
    c->SaveAs(outPng.c_str());
    c->SaveAs(outRoot.c_str());
    c->Close();
    fmain->Close();
    delete hRawYieldSys;
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