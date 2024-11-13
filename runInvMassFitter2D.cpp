#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooCurve.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TH2.h"

#include "InvMassFitter2D.h"

using namespace RooFit;
using std::cout;
using std::endl;

void runInvMassFitter2D() {
    //const char *fname = "~/MyMacros/Correlations/AO2D/AO2D_highIR.root";
    //const char *fname = "~/MyMacros/Correlations_v2/triggered_data/AO2D_LHC24_pass1.root";
    const char *fname = "~/MyMacros/Correlations_v2/triggered_data/AO2D_LHC23_pass4.root";

    TFile *file = TFile::Open(fname, "read");
    TFile *foutLS = TFile::Open("test_LS_triggered_LHC23_pass4.root", "RECREATE");
    TFile *foutOS = TFile::Open("test_OS_triggered_LHC23_pass4.root", "RECREATE");
    TFile *fEfficiencies = TFile::Open("Eff_times_Acc_Map.root", "read");

    if (!file) {
        cout << ">> ERROR File not well readout" << endl;
        return;
    }

    TFile *fAnaRes = TFile::Open("~/MyMacros/Correlations_v2/triggered_data/AnalysisResults_LHC23_pass4.root");

    TString dirnameBcCount = "bc-selection-task";
    TDirectory *dirBcData = (TDirectory *)fAnaRes->Get(dirnameBcCount);
    TH1F *hCounterTVX = (TH1F *)dirBcData->Get("hCounterTVX");
    cout << "Number of TVX entries: " << hCounterTVX->GetEntries() << endl;
    float const LumiTVX = hCounterTVX->GetEntries()/(59.4*1000);    //  micro b^-1

    //TDirectoryFile *dir = (TDirectoryFile *)file->Get("DF_2261906152687232"); // AO2D_highIR.root
    //TDirectoryFile *dir = (TDirectoryFile *)file->Get("DF_2363808317917888"); // AO2D_LHC24.root
    TDirectoryFile *dir = (TDirectoryFile *)file->Get("DF_2298103932356800"); // AO2D_LHC24.root

    TTree *tree = (TTree *)dir->Get("O2d0pair");
    if (!tree) {
        cout << ">> ERROR Tree not well readout" << endl;
        return;
    }

    // load 1D fit results
    //TFile *file1DFit = TFile::Open("/home/andrea/MyMacros/Correlations_v2/rawYields_D0_correlations_fromTree_cut2GeV.root", "read");
    //TFile *file1DFit = TFile::Open("/home/andrea/MyMacros/Correlations_v2/rawYields_D0_correlations_fromTree_LHC24.root", "read");
    TFile *file1DFit = TFile::Open("/home/andrea/MyMacros/Correlations_v2/rawYields_D0_correlations_fromTree_LHC23_pass4_1GeV.root", "read");
    TH1F *hReflOverSgn = (TH1F *)file1DFit->Get("hReflectionOverSignal");
    double reflOverSgn = hReflOverSgn->GetBinContent(1);
    cout << "Reflection over signal: " << reflOverSgn << endl;

    TFile *fileEffInt = TFile::Open("/home/andrea/MyMacros/Correlations_v2/AccEffPreselD0ToKPi_correlations_integrated.root", "read");
    TH1F *hEffInt = (TH1F *)fileEffInt->Get("hAccEffPreselD0ToKPi_pos0All");
    double integratedEff = hEffInt->GetBinContent(1);

    //TFile *fileWorkspace = TFile::Open("/home/andrea/MyMacros/Correlations_v2/workspace_massFitter_cut2GeV.root", "read");
    //TFile *fileWorkspace = TFile::Open("/home/andrea/MyMacros/Correlations_v2/workspace_massFitter_LHC24.root", "read");
    TFile *fileWorkspace = TFile::Open("/home/andrea/MyMacros/Correlations_v2/workspace_massFitter_cut1GeV.root", "read");
    if (!fileWorkspace || fileWorkspace->IsZombie()) {
        std::cerr << "Error opening workspace file!" << std::endl;
        return;
    }
    RooWorkspace* w = (RooWorkspace*)fileWorkspace->Get("wOut");
    if (!w) {
        std::cerr << "Workspace not found in file!" << std::endl;
        fileWorkspace->Close();
        return;
    }
    RooDataSet* dataSaved = (RooDataSet*)w->data("dataToSave");

    // Get all variables in the workspace
    const RooArgSet* varsSaved = dataSaved->get();

    InvMassFitter2D fitterLS(tree, "LS");
    InvMassFitter2D fitterOS(tree, "OS");
    cout << "fitter objects created and tree data loaded" << endl;

    //TH2F *hEffMap = dynamic_cast<TH2F *>(fEfficiencies->Get("hEfficiencyMapAll"));
    TH2F *hEffMap = dynamic_cast<TH2F *>(fEfficiencies->Get("reducedEfficiencyMap"));
    if (!hEffMap) {
        std::cerr << "Error: Histogram 'reducedEfficiencyMap' not found or not a TH2F." << std::endl;
        return;
    }

    fitterLS.setPtLims(1., 24.);
    fitterOS.setPtLims(1., 24.);

    fitterLS.setLumi(LumiTVX);
    fitterOS.setLumi(LumiTVX);

    fitterLS.setMassLims(1.7, 2.05);
    fitterOS.setMassLims(1.7, 2.05);

    fitterLS.set1DParameters(varsSaved, reflOverSgn, integratedEff);
    fitterOS.set1DParameters(varsSaved, reflOverSgn, integratedEff);

    fitterLS.setEfficiencyMap(hEffMap);
    fitterOS.setEfficiencyMap(hEffMap);

    fitterLS.do2DFit(true, true, foutLS);
    fitterOS.do2DFit(true, true, foutOS);

    cout << "Programa terminado" << endl;

}
