#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TDirectoryFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TString.h"
#include "TSystem.h"

#include "nlohmann/json.hpp"
#include "InvMassFitter2D.h"

using namespace RooFit;
using json = nlohmann::json;
using std::cout;
using std::endl;

bool verbose;
void runInvMassFitter2D() {

    // Load JSON config
    ifstream jsonFile("config_invMassFitter2D.json");
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

    // Use .get<std::string>() before assigning to TString
    TString const AO2D_data_name  = TString(jsonData["files"]["AO2D_data"].get<std::string>());

    TFile *AO2D_data = TFile::Open(AO2D_data_name, "read");
    if (!AO2D_data) {
        cout << ">> ERROR: AO2D_data not well readout" << endl;
        return;
    }
    TKey *key = (TKey*)AO2D_data->GetListOfKeys()->At(0);
    TDirectoryFile *dir = (TDirectoryFile*)AO2D_data->Get(key->GetName());

    TTree *tree = (TTree *)dir->Get("O2d0pair");
    if (!tree) {
        cout << ">> ERROR Tree not well readout" << endl;
        return;
    }

    TString const efficiencyMap_name  = TString(jsonData["files"]["efficiency_map"].get<std::string>());
    TFile *fEfficiencies = TFile::Open(efficiencyMap_name, "read");
    if (!fEfficiencies) {
        cout << ">> ERROR: efficiency map file not well readout" << endl;
        return;
    }

    TString const anaRes_data_name  = TString(jsonData["files"]["AnalysisResults_data"].get<std::string>());
    TFile *fAnaRes = TFile::Open(anaRes_data_name, "read");
    if (!fAnaRes) {
        cout << ">> ERROR: analysis results data not well readout" << endl;
        return;
    }

    TString const singleD0_fitResults_name  = TString(jsonData["files"]["singleD0_fitResults"].get<std::string>());
    TFile *file1DFit = TFile::Open(singleD0_fitResults_name, "read");
    if (!file1DFit) {
        cout << ">> ERROR: single D fit not well readout" << endl;
        return;
    }

    TString const intEff_name  = TString(jsonData["files"]["integrated_efficiency"].get<std::string>());
    TFile *fileEffInt = TFile::Open(intEff_name, "read");
    if (!fileEffInt) {
        cout << ">> ERROR: integrated efficiency file not well readout" << endl;
        return;
    }

    if (verbose) {
        cout << "Loaded files." << endl;
    }

    std::vector<std::pair<double, double>> ptSingleRanges;
    std::vector<std::pair<double, double>> ptPairRange;
    std::vector<std::pair<double, double>> massRanges;

    auto loadRanges = [](const nlohmann::json& j, const std::string& key) {
        std::vector<std::pair<double, double>> ranges;
        for (const auto& range : j["ranges"][key]) {
            if (range.is_array() && range.size() == 2) {
                ranges.emplace_back(range[0].get<double>(), range[1].get<double>());
            }
        }
        return ranges;
    };

    ptSingleRanges = loadRanges(jsonData, "ptSingle");
    ptPairRange = loadRanges(jsonData, "ptPair");
    massRanges = loadRanges(jsonData, "mass");

    std::string ptPairLimStrFirst = Form("%.2f", ptPairRange[0].first);
    std::string ptPairLimStrSecond = Form("%.2f", ptPairRange[0].second);

    std::string bkgFunc = jsonData["functions"]["bkgFunc"].get<std::string>();
    std::string sgnFunc = jsonData["functions"]["sgnFunc"].get<std::string>();
    std::string reflFunc = jsonData["functions"]["reflFunc"].get<std::string>();

    if (verbose) {
        cout << "  Selected functions:" << endl;
        cout << "   - Background function: " << bkgFunc << endl;
        cout << "   - Signal function: " << sgnFunc << endl;
        cout << "   - Reflected function: " << reflFunc << endl;
        cout << endl;
    }

        TH1F *hReflOverSgn = (TH1F *)file1DFit->Get("hReflectionOverSignal");
    double reflOverSgn = hReflOverSgn->GetBinContent(1);
    //double reflOverSgn = 0.0;
    if (verbose) {
        cout << "Reflection over signal: " << reflOverSgn << endl;
    }

    TH1F *hEffInt = (TH1F *)fileEffInt->Get("hAccEffPreselD0ToKPi_bothBDTs_no_ambiguous_pt_integratedAll");
    double integratedEff = hEffInt->GetBinContent(1);

    //TH2F *hEffMap = dynamic_cast<TH2F *>(fEfficiencies->Get("hEfficiencyMapAll"));
    TH2D *hEffMap = dynamic_cast<TH2D *>(fEfficiencies->Get("hEfficiencyMap"));
    if (!hEffMap) {
        std::cerr << "Error: Histogram 'efficiency map' not found or not a TH2F." << std::endl;
        return;
    }

    std::vector<bool> removeAmbiguous = jsonData["options"]["removeAmbiguous"];
    std::vector<bool> fitReflections = jsonData["options"]["fitReflections"];
    bool doFit = jsonData["options"]["doFit"];
    bool analyseKinematics = jsonData["options"]["analyseKinematics"];
    bool isMc = jsonData["options"]["isMC"];

    for (size_t i = 0; i < ptSingleRanges.size(); ++i) {
        for (const auto& mass : massRanges) {
            const auto& ptSingle = ptSingleRanges[i];

            std::string ptSingleLimStrFirst = Form("%.f", ptSingle.first);
            std::string ptSingleLimStrSecond = Form("%.f", ptSingle.second);

            std::string massLimStrFirst = Form("%.2f", mass.first);
            std::string massLimStrSecond = Form("%.2f", mass.second);

            if (verbose) {
                std::cout << "Running fit for:" << std::endl;
                std::cout << "  ptSingle: [" << ptSingle.first << ", " << ptSingle.second << "]" << std::endl;
                std::cout << "  ptPair:   [" << ptPairRange[0].first << ", " << ptPairRange[0].second << "]" << std::endl;
                std::cout << "  mass:     [" << mass.first << ", " << mass.second << "]" << std::endl;
            }

            std::string filenameLS = "correlations_LS_" + bkgFunc + "_" +
                                    massLimStrFirst + "_" + massLimStrSecond + "_pt_" + ptSingleLimStrFirst + "_" + ptSingleLimStrSecond + ".root";
            std::string filenameOS = "correlations_OS_" + bkgFunc + "_" +
                                    massLimStrFirst + "_" + massLimStrSecond + "_pt_" + ptSingleLimStrFirst + "_" + ptSingleLimStrSecond + ".root";

            TFile *foutLS = TFile::Open(filenameLS.c_str(), "RECREATE");
            TFile *foutOS = TFile::Open(filenameOS.c_str(), "RECREATE");

            if (verbose) {
                cout << " !! Saving LS results in " << filenameLS << endl;
                cout << " !! Saving OS results in " << filenameOS << endl;
            }

            InvMassFitter2D fitterLS(tree, "LS");
            InvMassFitter2D fitterOS(tree, "OS");
            cout << "fitter objects created and tree data loaded" << endl;

    
            fitterLS.removeAmbiguous(removeAmbiguous[i]);
            fitterOS.removeAmbiguous(removeAmbiguous[i]);

            fitterLS.setPtLims(ptSingle.first, ptSingle.second);
            fitterOS.setPtLims(ptSingle.first, ptSingle.second);

            fitterLS.setPtPairLims(ptPairRange[0].first, ptPairRange[0].second);
            fitterOS.setPtPairLims(ptPairRange[0].first, ptPairRange[0].second);

            fitterLS.setMassLims(mass.first, mass.second);
            fitterOS.setMassLims(mass.first, mass.second);

            // Set signal function for fit: gaus, CB
            // Both work fine for data, but for MC CB works significantly better
            fitterLS.setSgnFunc(sgnFunc);
            fitterOS.setSgnFunc(sgnFunc);
            // Set background function for fit: expo, poly0, poly1, poly2, expPoly1, expPoly2, expPoly1, exp2
            fitterLS.setBkgFunc(bkgFunc);
            fitterOS.setBkgFunc(bkgFunc);
            // Set reflection function for fit: gaus, doubleGaus
            fitterLS.setReflFunc(reflFunc);
            fitterOS.setReflFunc(reflFunc);


            fitterLS.set1DParameters(reflOverSgn, integratedEff);
            fitterOS.set1DParameters(reflOverSgn, integratedEff);

            fitterLS.setEfficiencyMap(hEffMap);
            fitterOS.setEfficiencyMap(hEffMap);

            // do2DFit(Bool_t draw, Bool_t doReflections, Bool_t isMc, TFile *fout);

            if (doFit) {
                fitterLS.do2DFit(true, fitReflections[i], isMc, foutLS);
                fitterOS.do2DFit(true, fitReflections[i], isMc, foutOS);
            }

            if (analyseKinematics) {
                fitterLS.analyseKinematicDistributions(foutLS, false, "beforeEffs_LS");
                fitterLS.analyseKinematicDistributions(foutLS, true, "afterEffs_LS");

                fitterOS.analyseKinematicDistributions(foutOS, false, "beforeEffs_OS");
                fitterOS.analyseKinematicDistributions(foutOS, true, "afterEffs_OS");
            }
        }
    }

    cout << "Programa terminado" << endl;

}
