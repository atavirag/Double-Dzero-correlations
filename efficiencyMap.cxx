#include "TStyle.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TMath.h"
#include "TSystem.h"
#include "nlohmann/json.hpp"

using namespace std;
using json = nlohmann::json;

bool verbose, debug;
void efficiencyMap() {
    // Load JSON config
    ifstream jsonFile("config-efficiencyMap.json");
    if (!jsonFile.is_open()) {
        cerr << "Error: Could not open JSON file." << endl;
        return;
    }

    // Read the JSON file into a string.
    std::string jsonString;
    jsonFile.seekg(0, std::ios::end);
    jsonString.reserve(jsonFile.tellg());
    jsonFile.seekg(0, std::ios::beg);
    jsonString.assign((std::istreambuf_iterator<char>(jsonFile)), std::istreambuf_iterator<char>());

    // Parse the JSON data.
    json jsonData = json::parse(jsonString);
    verbose = jsonData["verbose"];
    debug = jsonData["debug"];

    TString const fileNameRecoPrompt = jsonData["files"]["filenameRecoPrompt"];
    TString const fileNameRecoNonPrompt = jsonData["files"]["filenameRecoNonPrompt"];
    TString const fileNameGenPrompt = jsonData["files"]["filenameGenPrompt"];
    TString const fileNameGenNonPrompt = jsonData["files"]["filenameGenNonPrompt"];
    std::string const outputDirName = jsonData["outputDir"];
    TString const hPromptName = jsonData["histos"]["hPromptName"];
    TString const hNonPromptName = jsonData["histos"]["hNonPromptName"];
    TString const hAllName = jsonData["histos"]["hAllName"];

    // Load data
    TFile *fRecoPrompt = TFile::Open(fileNameRecoPrompt);
    TFile *fRecoNonPrompt = TFile::Open(fileNameRecoNonPrompt);
    TFile *fGenPrompt = TFile::Open(fileNameGenPrompt);
    TFile *fGenNonPrompt = TFile::Open(fileNameGenNonPrompt);
    if (!fRecoPrompt || !fRecoNonPrompt) {
        cerr << "ERROR: files reco not found\n";
        return;
    }
    if (!fGenPrompt || !fGenNonPrompt) {
        cerr << "ERROR: files gen not found\n";
        return;
    }
     auto outputName = outputDirName + "Eff_times_Acc_Map.root";
    TFile *outputFile = new TFile(outputName.c_str(), "RECREATE");
    outputFile->cd();

    TString dirnameGen = "hf-task-d0";
    TDirectory *dirGenPrompt = (TDirectory *)fGenPrompt->Get(dirnameGen);
    TDirectory *dirGenNonPrompt = (TDirectory *)fGenNonPrompt->Get(dirnameGen);
    TDirectory *dirRecoPrompt = (TDirectory *)fRecoPrompt->Get(dirnameGen);
    TDirectory *dirRecoNonPrompt = (TDirectory *)fRecoNonPrompt->Get(dirnameGen);

    // Load MC histos
    TH2F *hPtVsYGenPrompt = (TH2F *)dirGenPrompt->Get("hPtVsYGenPrompt");
    TH2F *hPtVsYGenNonPrompt = (TH2F *)dirGenNonPrompt->Get("hPtVsYGenNonPrompt");
    TH2F *hPtVsYGenAll = (TH2F *)dirGenNonPrompt->Get("hPtVsYGen");

    TH2F *hPtVsYRecoPrompt = (TH2F *)dirRecoPrompt->Get(hPromptName);
    TH2F *hPtVsYRecoNonPrompt = (TH2F *)dirRecoNonPrompt->Get(hNonPromptName);
    TH2F *hPtVsYRecoAll = (TH2F *)dirRecoNonPrompt->Get(hAllName);
    if (verbose) {
        cout << "            Data loaded \n";
    }

    // Get number of bins
    int nbinsPt = hPtVsYGenPrompt->GetXaxis()->GetNbins();
    int const nbinsY = hPtVsYGenPrompt->GetYaxis()->GetNbins();

    TH2F *hEfficiencyMapPrompt = (TH2F *)hPtVsYRecoPrompt->Clone();
    hEfficiencyMapPrompt->Reset();
    hEfficiencyMapPrompt->SetNameTitle("hEfficiencyMapPrompt", "Prompt efficiency map");

    TH2F *hEfficiencyMapNonPrompt = (TH2F *)hPtVsYRecoNonPrompt->Clone();
    hEfficiencyMapNonPrompt->Reset();
    hEfficiencyMapNonPrompt->SetNameTitle("hEfficiencyMapNonPrompt", "Non-prompt efficiency map");

    TH2F *hEfficiencyMapAll = (TH2F *)hPtVsYRecoAll->Clone();
    hEfficiencyMapAll->Reset();
    hEfficiencyMapAll->SetNameTitle("hEfficiencyMapAll", "efficiency map");

    if (debug) {
        std::cout << "Number of pT bins: " << nbinsPt << std::endl;
        std::cout << "Min pT bin: " << hPtVsYGenPrompt->GetXaxis()->GetBinLowEdge(0)
        << "     ; Max pT bin: " << hPtVsYGenPrompt->GetXaxis()->GetBinUpEdge(nbinsPt) << std::endl;
        std::cout << "Number of bins Y: " << nbinsY << std::endl;
    }

    // Divide histos point by point
    for (int ipt = 1; ipt <= nbinsPt; ipt++) {
        for (int iy = 1; iy <= nbinsY; iy++) {
            auto effPrompt = 0., effNonPrompt = 0., effAll = 0.;
            auto effPromptError = 0., effNonPromptError = 0., effAllError = 0.;

            if (hPtVsYRecoPrompt->GetBinContent(ipt, iy) != 0 && hPtVsYGenPrompt->GetBinContent(ipt, iy) != 0) {
                effPrompt = hPtVsYRecoPrompt->GetBinContent(ipt, iy) / hPtVsYGenPrompt->GetBinContent(ipt, iy);
                effPromptError = effPrompt * sqrt(pow(hPtVsYRecoPrompt->GetBinError(ipt, iy)/hPtVsYRecoPrompt->GetBinContent(ipt, iy), 2)
                          + pow(hPtVsYGenPrompt->GetBinError(ipt, iy)/hPtVsYGenPrompt->GetBinContent(ipt, iy), 2));
                hEfficiencyMapPrompt->SetBinContent(ipt, iy, effPrompt);
                hEfficiencyMapPrompt->SetBinError(ipt, iy, effPromptError);
            }
            if (hPtVsYRecoNonPrompt->GetBinContent(ipt, iy) != 0 && hPtVsYGenNonPrompt->GetBinContent(ipt, iy) != 0) {
                effNonPrompt = hPtVsYRecoNonPrompt->GetBinContent(ipt, iy) / hPtVsYGenNonPrompt->GetBinContent(ipt, iy);
                effNonPromptError = effNonPrompt * sqrt(pow(hPtVsYRecoNonPrompt->GetBinError(ipt, iy)/hPtVsYRecoNonPrompt->GetBinContent(ipt, iy), 2)
                          + pow(hPtVsYGenNonPrompt->GetBinError(ipt, iy)/hPtVsYGenNonPrompt->GetBinContent(ipt, iy), 2));
                hEfficiencyMapNonPrompt->SetBinContent(ipt, iy, effNonPrompt);
                hEfficiencyMapNonPrompt->SetBinError(ipt, iy, effNonPromptError);
            }
            if (hPtVsYRecoAll->GetBinContent(ipt, iy) != 0 && hPtVsYGenAll->GetBinContent(ipt, iy) != 0) {
                effAll = hPtVsYRecoAll->GetBinContent(ipt, iy) / hPtVsYGenAll->GetBinContent(ipt, iy);
                effAllError = effAll * sqrt(pow(hPtVsYRecoAll->GetBinError(ipt, iy)/hPtVsYRecoAll->GetBinContent(ipt, iy), 2)
                          + pow(hPtVsYGenAll->GetBinError(ipt, iy)/hPtVsYGenAll->GetBinContent(ipt, iy), 2));
                hEfficiencyMapAll->SetBinContent(ipt, iy, effAll);
                hEfficiencyMapAll->SetBinError(ipt, iy, effAllError);
            }
        }
    }

    // Plot efficienct map
    /* TCanvas *cPrompt = new TCanvas("cPrompt", "prompt canvas", 1000, 700);
    hEfficiencyMapPrompt->Draw("colz");
    cPrompt->SetLogz();

    TCanvas *cNonPrompt = new TCanvas("cNonPrompt", "non-prompt canvas", 1000, 700);
    hEfficiencyMapNonPrompt->Draw("colz");
    cNonPrompt->SetLogz(); */

    outputFile->cd();
    hEfficiencyMapPrompt->Write();
    hEfficiencyMapNonPrompt->Write();
    hEfficiencyMapAll->Write();
    outputFile->Close();

    cout << "Plots saved. Exit!" << endl;
}