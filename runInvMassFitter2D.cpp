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
#include "TH2D.h"
#include "RooBernstein.h"
#include "TKey.h"


#include "InvMassFitter2D.h"

using namespace RooFit;
using std::cout;
using std::endl;

void runInvMassFitter2D() {
    //const char *fname =  "~/MyMacros/Correlations_v2/triggered_data/AO2D_Data_Full2024.root";
    const char *fname =  "~/MyMacros/corelations_with_phi/AO2D_fullData_oldBDTs.root";
    //const char *fname =  "~/MyMacros/Double-Dzero-correlations/Datasets/AO2D_LHC24k3_full.root";

    TFile *file = TFile::Open(fname, "read");
    TFile *foutLS = TFile::Open("~/MyMacros/corelations_with_phi/correlations_LS_oldBDTs_main.root", "RECREATE");
    TFile *foutOS = TFile::Open("~/MyMacros/corelations_with_phi/correlations_OS_oldBDTs_main.root", "RECREATE");
    TFile *fEfficiencies = TFile::Open("~/MyMacros/Double-Dzero-correlations/Datasets/Eff_times_Acc_Map_weighted_no_ambiguous.root", "read");

    if (!file) {
        cout << ">> ERROR File not well readout" << endl;
        return;
    }

    //TFile *fAnaRes = TFile::Open("~/MyMacros/Correlations_github/AnalysisResults_LHC23_pass4.root");
    TFile *fAnaRes = TFile::Open("~/MyMacros/Correlations_v2/triggered_data/AnalysisResults_Data_Full2024.root");

    TKey *key = (TKey*)file->GetListOfKeys()->At(0);
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(key->GetName());

    TTree *tree = (TTree *)dir->Get("O2d0pair");
    if (!tree) {
        cout << ">> ERROR Tree not well readout" << endl;
        return;
    }

    // load 1D fit results
    TFile *file1DFit = TFile::Open("/home/andrea/MyMacros/Correlations_github/rawYields_LHC23_pass4.root", "read");
    //TFile *file1DFit = TFile::Open("/home/andrea/MyMacros/Correlations_v2/rawYields_LHC24_full_pt23.root", "read");
    TH1F *hReflOverSgn = (TH1F *)file1DFit->Get("hReflectionOverSignal");
    //double reflOverSgn = hReflOverSgn->GetBinContent(1);
    double reflOverSgn = 0.0;
    cout << "Reflection over signal: " << reflOverSgn << endl;

    TFile *fileEffInt = TFile::Open("/home/andrea/MyMacros/Common/AccEffPreselD0ToKPi_k3_SecondBDT_ptIntegrated.root", "read");
    TH1F *hEffInt = (TH1F *)fileEffInt->Get("hAccEffPreselD0ToKPi_k3_SecondBDT_ptIntegratedAll");
    double integratedEff = hEffInt->GetBinContent(1);

    //TFile *fileWorkspace = TFile::Open("/home/andrea/MyMacros/Correlations_github/workspace_massFitter_LHC23_pass4.root", "read");
    TFile *fileWorkspace = TFile::Open("/home/andrea/MyMacros/Correlations_github/workspace_massFitter_full2024.root", "read");
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

    if (dataSaved && varsSaved)
    {
        // Loop over all variables in the dataset
        RooFIter iter = varsSaved->fwdIterator();
        RooAbsArg *var = nullptr;

        while ((var = iter.next()))
        {
            // Cast to RooRealVar to access value, if applicable
            RooRealVar *realVar = dynamic_cast<RooRealVar *>(var);
            if (realVar)
            {
                std::cout << "Variable: " << realVar->GetName()
                          << ", Value: " << realVar->getVal() << std::endl;
            }
            else
            {
                std::cout << "Variable: " << var->GetName()
                          << " is not a RooRealVar." << std::endl;
            }
        }
    }
    else
    {
        std::cerr << "Data or variables set is null!" << std::endl;
    }

    InvMassFitter2D fitterLS(tree, "LS");
    InvMassFitter2D fitterOS(tree, "OS");
    cout << "fitter objects created and tree data loaded" << endl;

    //TH2F *hEffMap = dynamic_cast<TH2F *>(fEfficiencies->Get("hEfficiencyMapAll"));
    TH2F *hEffMap = dynamic_cast<TH2F *>(fEfficiencies->Get("reducedEfficiencyMap"));
    if (!hEffMap) {
        std::cerr << "Error: Histogram 'reducedEfficiencyMap' not found or not a TH2F." << std::endl;
        return;
    }

    fitterLS.removeAmbiguous(true);
    fitterOS.removeAmbiguous(true);

    fitterLS.setPtLims(1., 24.);
    fitterOS.setPtLims(1., 24.);

    fitterLS.setPtPairLims(-100., 100.);
    fitterOS.setPtPairLims(-100.0, 100.);

    //fitterLS.setMassLims(1.8, 1.95);
    //fitterOS.setMassLims(1.8, 1.95);
    fitterLS.setMassLims(1.74, 2.04);
    fitterOS.setMassLims(1.74, 2.04);

    // Set signal function for fit: gaus, CB
    // Both work fine for data, but for MC CB works significantly better
    fitterLS.setSgnFunc("gaus");
    fitterOS.setSgnFunc("gaus");
    // Set background function for fit: expo, poly0, poly1, poly2, cheby,expPoly1, expPoly2, expPoly1, bern, exp2
    fitterLS.setBkgFunc("expPoly2");
    fitterOS.setBkgFunc("expPoly2");
    // Set reflection function for fit: gaus, doubleGaus
    fitterLS.setReflFunc("gaus");
    fitterOS.setReflFunc("gaus");


    fitterLS.set1DParameters(varsSaved, reflOverSgn, integratedEff);
    fitterOS.set1DParameters(varsSaved, reflOverSgn, integratedEff);

    fitterLS.setEfficiencyMap(hEffMap);
    fitterOS.setEfficiencyMap(hEffMap);

    // do2DFit(Bool_t draw, Bool_t doReflections, Bool_t isMc, TFile *fout);
    fitterOS.do2DFit(true, false, false, foutOS);
    fitterLS.do2DFit(true, false, false, foutLS);

    cout << "Programa terminado" << endl;

}
