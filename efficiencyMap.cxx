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
    TString const fileNameGenAll = jsonData["files"]["filenameGenAll"];
    TString const fileNameRecoAll = jsonData["files"]["filenameRecoAll"];

    std::string const outputDirName = jsonData["outputDir"];

    TString const hPromptName = jsonData["histos"]["hPromptName"];
    TString const hNonPromptName = jsonData["histos"]["hNonPromptName"];
    TString const hAllName = jsonData["histos"]["hAllName"];

    // Load data
    TFile *fRecoPrompt = TFile::Open(fileNameRecoPrompt);
    TFile *fRecoNonPrompt = TFile::Open(fileNameRecoNonPrompt);
    TFile *fGenPrompt = TFile::Open(fileNameGenPrompt);
    TFile *fGenNonPrompt = TFile::Open(fileNameGenNonPrompt);
    TFile *fGenAll = TFile::Open(fileNameGenAll);
    TFile *fRecoAll = TFile::Open(fileNameRecoAll);

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
    TString dirname = "hf-correlator-d-meson-pairs";
    TDirectory *dirGenPrompt = (TDirectory *)fGenPrompt->Get(dirnameGen);
    TDirectory *dirGenNonPrompt = (TDirectory *)fGenNonPrompt->Get(dirnameGen);
    TDirectory *dirRecoPrompt = (TDirectory *)fRecoPrompt->Get(dirnameGen);
    TDirectory *dirRecoNonPrompt = (TDirectory *)fRecoNonPrompt->Get(dirnameGen);
    TDirectory *dirRecoAll = (TDirectory *)fRecoAll->Get(dirname);
    TDirectory *dirGenAll = (TDirectory *)fGenAll->Get(dirname);

    // Load MC histos
    TH2F *hPtVsYGenPrompt = (TH2F *)dirGenPrompt->Get("hPtVsYGenPrompt");
    TH2F *hPtVsYGenNonPrompt = (TH2F *)dirGenNonPrompt->Get("hPtVsYGenNonPrompt");
    TH2F *hPtVsYGenAll = (TH2F *)dirGenAll->Get("hPtVsYMcGen");
    //TH2F *hPtVsYGenAll = (TH2F *)dirGenAll->Get("hPtVsYGen");

    TH2F *hPtVsYRecoPrompt = (TH2F *)dirRecoPrompt->Get(hPromptName);
    TH2F *hPtVsYRecoNonPrompt = (TH2F *)dirRecoNonPrompt->Get(hNonPromptName);
    TH2F *hPtVsYRecoAll = (TH2F *)dirRecoAll->Get(hAllName);
    if (verbose) {
        cout << "            Data loaded \n";
    }

    // Get number of bins
    int nbinsPt = hPtVsYGenPrompt->GetXaxis()->GetNbins();
    int const nbinsY = hPtVsYGenPrompt->GetYaxis()->GetNbins();
    int y_min = hPtVsYGenAll->GetYaxis()->FindBin(-0.5);
    int y_max = hPtVsYGenAll->GetYaxis()->FindBin(0.5);

    int nbinsPtAll = hPtVsYGenAll->GetXaxis()->GetNbins();
    int const nbinsYAll = hPtVsYGenAll->GetYaxis()->GetNbins();

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
            auto effPrompt = 0., effNonPrompt = 0.;
            auto effPromptError = 0., effNonPromptError = 0.;

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
        }
    }

    for (int ipt = 1; ipt <= nbinsPtAll; ipt++) {
        for (int iy = y_min; iy < y_max; iy++) {
            auto effAll = 0.;
            auto effAllError = 0.;
    if (hPtVsYRecoAll->GetBinContent(ipt, iy) != 0 && hPtVsYGenAll->GetBinContent(ipt, iy) != 0) {
                effAll = hPtVsYRecoAll->GetBinContent(ipt, iy) / hPtVsYGenAll->GetBinContent(ipt, iy);
                effAllError = effAll * sqrt(pow(hPtVsYRecoAll->GetBinError(ipt, iy)/hPtVsYRecoAll->GetBinContent(ipt, iy), 2)
                          + pow(hPtVsYGenAll->GetBinError(ipt, iy)/hPtVsYGenAll->GetBinContent(ipt, iy), 2));
                hEfficiencyMapAll->SetBinContent(ipt, iy, effAll);
                hEfficiencyMapAll->SetBinError(ipt, iy, effAllError);
            }
        }
    }

    // Define new binning
int newNbinsX = 72; // For example, 0 to 36 GeV with 0.2 GeV width = 180 bins
int newNbinsY = hEfficiencyMapAll->GetNbinsY(); // Keep the same Y bins or adjust as needed
double newXMin = hEfficiencyMapAll->GetXaxis()->GetXmin();
double newXMax = hEfficiencyMapAll->GetXaxis()->GetXmax();

// Create a new histogram with fewer bins
TH2F* reducedEfficiencyMap = new TH2F("reducedEfficiencyMap", "Reduced Efficiency Map",
                                       newNbinsX, newXMin, newXMax,
                                       newNbinsY, -1., 1.);

// Loop over original histogram bins and fill new histogram
for (int i = 1; i <= nbinsPtAll; ++i) {
    for (int j = 1; j <= nbinsYAll; ++j) {
        // Get the contents of the original bin
        double content = hEfficiencyMapAll->GetBinContent(i, j);
        double error = hEfficiencyMapAll->GetBinError(i, j);
        
        // Determine the new bin index based on your new binning strategy
        int newBinX = (i - 1) / 5 + 1; // Merging every 2 original bins
        
        // Only fill the new histogram if the new bin index is valid
        if (newBinX <= newNbinsX) {
            // Add the content to the new histogram
            reducedEfficiencyMap->SetBinContent(newBinX, j, content);  // Accumulate content
            reducedEfficiencyMap->SetBinError(newBinX, j, 
                TMath::Sqrt(TMath::Power(reducedEfficiencyMap->GetBinError(newBinX, j), 2) + TMath::Power(error, 2))); // Accumulate error
        }
    }
}
    outputFile->cd();
    hEfficiencyMapPrompt->Write();
    hEfficiencyMapNonPrompt->Write();
    hEfficiencyMapAll->Write();
    reducedEfficiencyMap->Write();
    outputFile->Close();

    cout << "Plots saved. Exit!" << endl;

}