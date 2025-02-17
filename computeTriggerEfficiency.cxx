/// Quick macro to calculate the luminosity of the triggered dataset

#include <iostream>
#include <functional>

#include "TFile.h"
#include "TH1.h"
#include "TKey.h"


using std::cout;
using std::endl;

bool verbose = true;
bool debug = false;

struct HistogramValues {
    double analysedTriggersValue;
    double inspectedTVXValue;
    double scalerValue;
};

void computeTriggerEfficiency() {

    TFile *fileMassMB = TFile::Open("~/MyMacros/Correlations_github/trigger_efficiency/rawYields_MB.root", "read");
    TFile *fileMassTrigger = TFile::Open("~/MyMacros/Correlations_github/trigger_efficiency/rawYields_singleCharm.root", "read");
    TFile *fileAnaResMB = TFile::Open("~/MyMacros/Correlations_github/trigger_efficiency/AnalysisResults_2024_sampled_triggered.root", "read");
    TFile *fileAnaResTrigger = TFile::Open("~/MyMacros/Correlations_github/trigger_efficiency/AnalysisResults_singleCharm_triggered.root", "read");

    TFile *fout = TFile::Open("triggerEfficiency_singleCharm.root", "RECREATE");

    if (!fileMassMB || !fileMassTrigger) {
        cout << ">> ERROR Files with inv. mass fits not well readout" << endl;
        return;
    }
    if (!fileAnaResMB || !fileAnaResTrigger) {
        cout << ">> ERROR AnalysisResults files not well readout" << endl;
        return;
    }

    TDirectory *dir = dynamic_cast<TDirectory*>(fileAnaResTrigger->Get("hf-candidate-creator-2prong/Zorro"));
    // Map to store histogram values for each directory
    std::map<std::string, HistogramValues> directoryValues;

    // Recursive function to traverse directories and process histograms
    std::function<void(TDirectory*, const std::string&)> processDir =
        [&](TDirectory* dir, const std::string& path) {
            TIter nextKey(dir->GetListOfKeys());
            TKey* key;

            HistogramValues values = {-1, -1, -1}; // Default invalid values

            while ((key = (TKey*)nextKey())) {
                TObject* obj = key->ReadObj();

                if (obj->InheritsFrom("TDirectory")) {
                    // If it's a directory, recurse into it
                    TDirectory* subDir = (TDirectory*)obj;
                    std::string newPath = path + "/" + subDir->GetName();
                    processDir(subDir, newPath);
                } else if (obj->InheritsFrom("TH1")) {
                    // If it's a histogram, process it
                    TH1* hist = (TH1*)obj;

                    if (std::string(hist->GetName()) == "AnalysedTriggers") {
                        int binAnalysedTriggers = hist->GetXaxis()->FindBin("fHfSingleCharm2P");
                        values.analysedTriggersValue = hist->GetBinContent(binAnalysedTriggers);
                    } else if (std::string(hist->GetName()) == "InspectedTVX") {
                        values.inspectedTVXValue = hist->Integral();
                    } else if (std::string(hist->GetName()) == "Scalers") {
                        //int bin = hist->GetNbinsX(); // the last bin contains the triggered events
                        int bin = hist->GetXaxis()->FindBin("fHfSingleCharm2P"); // the last bin contains the triggered events
                        values.scalerValue = hist->GetBinContent(bin);
                    }
                }
            }

            // If any values were found, store them
            if (values.analysedTriggersValue != -1 || values.inspectedTVXValue != -1 || values.scalerValue != -1) {
                directoryValues[path] = values;
            }
        };

    // Start processing from the root directory
    processDir(dir, "");

    // Get normalisation of triggered data
    double totAnalysedTriggers = 0.;
    double totScalers = 0.;
    double totInspectedTVX = 0.;
    double totNevTrigger = 0.;
    for (const auto& [dirPath, values] : directoryValues) {

        double nEvents = (values.inspectedTVXValue * values.analysedTriggersValue) / values.scalerValue;
        //totNev += nEvents;

        totAnalysedTriggers += values.analysedTriggersValue;
        totScalers += values.scalerValue;
        totInspectedTVX += values.inspectedTVXValue;

        if (debug) {
            std::cout << "Directory: " << dirPath << std::endl;
            std::cout << "  AnalysedTriggersOfInterest: " << values.analysedTriggersValue << std::endl;
            std::cout << "  InspectedTVX: " << values.inspectedTVXValue << std::endl;
            std::cout << "  Scalers: " << values.scalerValue << std::endl;
            std::cout << "  Number of events: " << nEvents << std::endl;
        }
    }

    totNevTrigger = totInspectedTVX * totAnalysedTriggers / totScalers;

    if (verbose) {
        std::cout << "  NORMALISATION OF TRIGGERED DATA" << std::endl;
        std::cout << "  Total AnalysedTriggersOfInterest: " << totAnalysedTriggers << std::endl;
        std::cout << "  Total InspectedTVX: " << totInspectedTVX << std::endl;
        std::cout << "  Total Scalers: " << totScalers << std::endl;
        std::cout << "  Number of events: " << totNevTrigger << std::endl;
    }

    // Normalisation of MB data
    TDirectory *dirMB = dynamic_cast<TDirectory*>(fileAnaResMB->Get("bc-selection-task"));
    TH1F *hCounterTVX = (TH1F *)dirMB->Get("hCounterTVX");
    double countsTVX = hCounterTVX->GetEntries();

    if (verbose) {
        std::cout << "\n  NORMALISATION OF MB DATA" << std::endl;
        std::cout << "  Total TVX entries: " << countsTVX << std::endl;
    }
    double ratioTriggerMB = totNevTrigger / countsTVX;
    std::cout << "\n trigger / MB = " << totNevTrigger / countsTVX << endl; // Fabrizio dijo que el singleCharm tenía un downscale muy grande

    // Divide normalised raw yields
    TH1F *hRawYieldMB = (TH1F *)fileMassMB->Get("hRawYields");
    TH1F *hRawYieldTrigger = (TH1F *)fileMassTrigger->Get("hRawYields");

    if (!hRawYieldMB || !hRawYieldTrigger) {
        cerr << "ERROR: raw yields not found" << endl;
        return;
    }

    TH1F *hTriggerEfficiency = (TH1F *)hRawYieldMB->Clone();

    for (int i = 1; i <= hRawYieldMB->GetNbinsX(); i++) {
        double triggerEff = (hRawYieldTrigger->GetBinContent(i) / hRawYieldMB->GetBinContent(i)) / ratioTriggerMB;
        double triggerEffError = triggerEff * sqrt(pow(hRawYieldTrigger->GetBinError(i) / hRawYieldTrigger->GetBinContent(i), 2) +
                                      pow(hRawYieldMB->GetBinError(i)/ hRawYieldMB->GetBinContent(i), 2));
        cout << triggerEff<< endl;
        hTriggerEfficiency->SetBinContent(i, triggerEff);
        hTriggerEfficiency->SetBinError(i, triggerEffError);
    }

    fout->cd();
    TCanvas *c = new TCanvas("c", "c", 1000, 700);
    c->SetGrid();  // Add grid
    hTriggerEfficiency->SetTitle("Trigger efficiency");

    hTriggerEfficiency->SetLineColor(kBlue+3);
    hTriggerEfficiency->SetLineWidth(2);
    hTriggerEfficiency->SetMarkerStyle(20);
    hTriggerEfficiency->SetMarkerSize(1);
    hTriggerEfficiency->SetMarkerColor(kBlue+3);

    hTriggerEfficiency->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hTriggerEfficiency->GetYaxis()->SetTitle("Trigger Efficiency");

    hTriggerEfficiency->GetYaxis()->SetRangeUser(0.2, 1.4); // Adjust Y range
    hTriggerEfficiency->GetXaxis()->SetRangeUser(0, 16.0); // Adjust Y range

    // Draw a red line at y=1
    TLine *line = new TLine(hTriggerEfficiency->GetXaxis()->GetXmin(), 1, hTriggerEfficiency->GetXaxis()->GetXmax(), 1);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    //line->SetLineStyle(2);  // Dashed line

    hTriggerEfficiency->Draw();
    line->Draw("same");
    c->SaveAs("trigger_eff.png");
    hTriggerEfficiency->Write();
}
