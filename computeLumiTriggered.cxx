/// Quick macro to calculate the luminosity of the triggered dataset

#include <iostream>
#include <functional>

#include "TFile.h"
#include "TH1.h"
#include "TKey.h"


using std::cout;
using std::endl;

bool verbose = false;

struct HistogramValues {
    double analysedTriggersValue;
    double inspectedTVXValue;
    double scalerValue;
};

void computeLumiTriggered() {

    TFile *file = TFile::Open("~/MyMacros/Correlations_v2/triggered_data/AnalysisResults_Data_Full2024.root", "read");
    TFile *fout = TFile::Open("lumiTriggered_LHC24.root", "RECREATE");

    if (!file) {
        cout << ">> ERROR File not well readout" << endl;
        return;
    }

    TDirectory *dir = dynamic_cast<TDirectory*>(file->Get("hf-candidate-creator-2prong/Zorro"));
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

                    if (std::string(hist->GetName()) == "AnalysedTriggersOfInterest") {
                        values.analysedTriggersValue = hist->GetBinContent(1);
                    } else if (std::string(hist->GetName()) == "InspectedTVX") {
                        values.inspectedTVXValue = hist->GetBinContent(1);
                    } else if (std::string(hist->GetName()) == "Scalers") {
                        const char* triggerLabel = "fHfDoubleCharm2P";
                        int bin = hist->GetXaxis()->FindBin(triggerLabel);
                        if (bin > 0 && bin <= hist->GetNbinsX()) {
                            values.scalerValue = hist->GetBinContent(bin);
                        } else {
                            std::cout << ">> ERROR: Label '" << triggerLabel << "' not found in Scalers histogram" << std::endl;
                        }
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

    // Example: Print out stored values
    double totNev = 0.;
    for (const auto& [dirPath, values] : directoryValues) {

        double nEvents = values.inspectedTVXValue * (values.analysedTriggersValue / values.scalerValue);
        totNev += nEvents;

        if (verbose) {
            std::cout << "Directory: " << dirPath << std::endl;
            std::cout << "  AnalysedTriggersOfInterest: " << values.analysedTriggersValue << std::endl;
            std::cout << "  InspectedTVX: " << values.inspectedTVXValue << std::endl;
            std::cout << "  Scalers[fHfDoubleCharm2P]: " << values.scalerValue << std::endl;
            std::cout << "  Number of events: " << nEvents << std::endl;
        }
    }

    // TVX cross-section
    float xSecTVX = 59.4; // mb
    // Integrated luminosity TVX = Nevents / xSecTVX
    double lumi = totNev / (xSecTVX * pow(10.,9)); // pb⁻1
    if (verbose) {
        std::cout << "  Integrated Luminosity TVX BEFORE BC CORRECTION: " << lumi << std::endl;
    }

    // Apply correction for the BC cuts
    TDirectory *dirBC = dynamic_cast<TDirectory*>(file->Get("bc-selection-task"));
    TH1F *hCounterTVX = (TH1F *)(dirBC->Get("hCounterTVX"));
    TH1F *hCounterTVXafterBCcuts = (TH1F *)(dirBC->Get("hCounterTVXafterBCcuts"));

    // Get the number of entries in the counter histos
    double countsTVX = hCounterTVX->GetEntries();
    double countsTVXafterBCcuts = hCounterTVXafterBCcuts->GetEntries();

    // Corrected integrated lumi: (Nevents / xSecTVX) * (countsTVX / countsTVXadterBCcuts)
    double corrLumi = lumi * (countsTVX / countsTVXafterBCcuts); // pb⁻1
    std::cout << "  Integrated Luminosity: " << corrLumi << std::endl;
}
