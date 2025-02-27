#include "TStyle.h"
#include "TString.h"
#include "TH2F.h"
#include "TH3F.h"
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
void efficiencyMap_PVContrib()
{
    // Load JSON config
    ifstream jsonFile("config-efficiencyMap_PVContrib.json");
    if (!jsonFile.is_open())
    {
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

    TString const fileNameGen = jsonData["files"]["filenameGen"];
    TString const fileNameReco = jsonData["files"]["filenameReco"];
    TString const fileNameData = jsonData["files"]["filenameData"];

    std::string const outputDirName = jsonData["outputDir"];

    TString const hNameReco = jsonData["histos"]["hNameReco"];
    TString const hNameGen = jsonData["histos"]["hNameGen"];
    TString const hNameData = jsonData["histos"]["hNameData"];

    // Load data
    TFile *fGen = TFile::Open(fileNameGen);
    TFile *fReco = TFile::Open(fileNameReco);
    TFile *fData = TFile::Open(fileNameData);

    if (!fReco || !fGen || !fData)
    {
        cerr << "ERROR: files not found\n";
        return;
    }

    auto outputName = outputDirName + "Eff_times_Acc_Map_weighted.root";
    TFile *outputFile = new TFile(outputName.c_str(), "RECREATE");
    outputFile->cd();

    //TString dirnameGen = "hf-task-d0";
    TString dirname = "hf-correlator-d-meson-pairs";

    TDirectory *dirReco = (TDirectory *)fReco->Get(dirname);
    TDirectory *dirGen = (TDirectory *)fGen->Get(dirname);
    TDirectory *dirData = (TDirectory *)fData->Get(dirname);

    // Load MC histos
    TH3F *hPtVsYReco = (TH3F *)dirReco->Get(hNameReco);
    TH3F *hPtVsYGen = (TH3F *)dirGen->Get(hNameGen);
    // Load data histo
    TH1F *hNContrib = (TH1F *)dirData->Get(hNameData);


    double dataIntegral = hNContrib->Integral();
    double recoIntegral = hPtVsYReco->ProjectionZ()->Integral();
    double genIntegral = hPtVsYGen->ProjectionZ()->Integral();

    TH1F *hNContribReco = (TH1F *)hPtVsYReco->ProjectionZ();
    TH1F *hNContribGen = (TH1F *)hPtVsYGen->ProjectionZ();

    if (verbose)
    {
        cout << "            Data loaded \n";
    }

    // Get number of bins
    int const nbinsPt = hPtVsYGen->GetXaxis()->GetNbins();
    int const nbinsY = hPtVsYGen->GetYaxis()->GetNbins();
    int const nbinsNContrib = hPtVsYGen->GetZaxis()->GetNbins();

    int y_min = hPtVsYGen->GetYaxis()->FindBin(-0.5);
    int y_max = hPtVsYGen->GetYaxis()->FindBin(0.5);

    TH2F *hEfficiencyMap = (TH2F *)hPtVsYReco->Project3D("yx");
    hEfficiencyMap->Reset();
    hEfficiencyMap->SetNameTitle("hEfficiencyMap", "efficiency map");

    TH2F *hPtVsYRecoWeighted = (TH2F *)hPtVsYReco->Project3D("yx");
    hPtVsYRecoWeighted->Reset();
    hPtVsYRecoWeighted->SetNameTitle("hPtVsYRecoWeighted", "hPtVsYRecoWeighted");

    TH2F *hPtVsYGenWeighted = (TH2F *)hPtVsYGen->Project3D("yx");
    hPtVsYGenWeighted->Reset();
    hPtVsYGenWeighted->SetNameTitle("hPtVsYGenWeighted", "hPtVsYGenWeighted");

    if (debug)
    {
        std::cout << "Number of pT bins: " << nbinsPt << std::endl;
        std::cout << "Min pT bin: " << hPtVsYGen->GetXaxis()->GetBinLowEdge(0)
                  << "     ; Max pT bin: " << hPtVsYGen->GetXaxis()->GetBinUpEdge(nbinsPt) << std::endl;
        std::cout << "Number of bins Y: " << nbinsY << std::endl;
        std::cout << "Number of PV contrib bins: " << nbinsNContrib << std::endl;
    }

    // Divide histos point by point
    for (int ipt = 1; ipt <= nbinsPt; ipt++)
    {
        for (int iy = y_min; iy < y_max; iy++)
        {
            // Weight the distributions using nContrib info
            double weightSumReco = 0.0, weightSumGen = 0.0;
            for (int in = 1; in <= nbinsNContrib; in++)
            {
                double binContentReco = hPtVsYReco->GetBinContent(ipt, iy, in);
                double binContentGen = hPtVsYGen->GetBinContent(ipt, iy, in);

                double nContribValueData = hNContrib->GetBinContent(in);
                double nContribValueMCReco = hNContribReco->GetBinContent(in); // MC shape in NContrib
                double nContribValueMCGen = hNContribGen->GetBinContent(in); // MC shape in NContrib

                double dataRecoWeight = 1.0;
                double dataGenWeight = 1.0;
                if (nContribValueData > 0) {
                    dataRecoWeight = (nContribValueMCReco / recoIntegral) / (nContribValueData/dataIntegral);
                    dataGenWeight = (nContribValueMCGen / genIntegral) / (nContribValueData/dataIntegral);
                }

                weightSumReco += (dataRecoWeight * binContentReco);
                weightSumGen += ( dataGenWeight * binContentGen);
            }
            // Compute the weighted average (or leave as summed weights)
            if (valueSumReco > 0)
            {
                hPtVsYRecoWeighted->SetBinContent(ipt, iy, weightSumReco);
                cout << "weightSumReco " << weightSumReco << " valueSumReco " << valueSumReco << endl;
            }
            if (valueSumGen > 0)
            {
                hPtVsYGenWeighted->SetBinContent(ipt, iy, weightSumGen);
            }
        }
    }


    // Divide histos point by point
    for (int ipt = 1; ipt <= nbinsPt; ipt++)
    {
        for (int iy = y_min; iy < y_max; iy++)
        {
            auto eff = 0.;
            auto effError = 0.;

            if (hPtVsYRecoWeighted->GetBinContent(ipt, iy) != 0 && hPtVsYGenWeighted->GetBinContent(ipt, iy) != 0)
            {
                eff = hPtVsYRecoWeighted->GetBinContent(ipt, iy) / hPtVsYGenWeighted->GetBinContent(ipt, iy);
                effError = eff * sqrt(pow(hPtVsYRecoWeighted->GetBinError(ipt, iy) / hPtVsYRecoWeighted->GetBinContent(ipt, iy), 2) + pow(hPtVsYGenWeighted->GetBinError(ipt, iy) / hPtVsYGenWeighted->GetBinContent(ipt, iy), 2));
                hEfficiencyMap->SetBinContent(ipt, iy, eff);
                hEfficiencyMap->SetBinError(ipt, iy, effError);
            }
        }
    }

    // Define new binning
    int newNbinsX = 72;                          // For example, 0 to 36 GeV with 0.2 GeV width = 180 bins
    int newNbinsY = hEfficiencyMap->GetNbinsY(); // Keep the same Y bins or adjust as needed
    double newXMin = hEfficiencyMap->GetXaxis()->GetXmin();
    double newXMax = hEfficiencyMap->GetXaxis()->GetXmax();

    // Create a new histogram with fewer bins
    TH2F *reducedEfficiencyMap = new TH2F("reducedEfficiencyMap", "Reduced Efficiency Map",
                                          newNbinsX, newXMin, newXMax,
                                          newNbinsY, -1., 1.);

    // Loop over original histogram bins and fill new histogram
    for (int i = 1; i <= nbinsPt; ++i)
    {
        for (int j = 1; j <= nbinsY; ++j)
        {
            // Get the contents of the original bin
            double content = hEfficiencyMap->GetBinContent(i, j);
            double error = hEfficiencyMap->GetBinError(i, j);

            // Determine the new bin index based on your new binning strategy
            int newBinX = (i - 1) / 5 + 1; // Merging every 2 original bins

            // Only fill the new histogram if the new bin index is valid
            if (newBinX <= newNbinsX)
            {
                // Add the content to the new histogram
                reducedEfficiencyMap->SetBinContent(newBinX, j, content); // Accumulate content
                reducedEfficiencyMap->SetBinError(newBinX, j,
                                                  TMath::Sqrt(TMath::Power(reducedEfficiencyMap->GetBinError(newBinX, j), 2) + TMath::Power(error, 2))); // Accumulate error
            }
        }
    }
    outputFile->cd();
    hEfficiencyMap->Write();
    hPtVsYRecoWeighted->Write();
    hPtVsYGenWeighted->Write();
    reducedEfficiencyMap->Write();
    outputFile->Close();

    cout << "Plots saved. Exit!" << endl;
}