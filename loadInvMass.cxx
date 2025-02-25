#include "TFile.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TObject.h"
#include "Riostream.h"

#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooCategory.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>

#include "TCanvas.h"
#include "RooGaussian.h"
#include "RooWorkspace.h"
#include "RooPlot.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace RooFit;

using std::cout;
using std::endl;

enum PairTypeSel
{
    DD = 0,   // D0-D0
    DbarDbar, // D0bar-D0bar
    DDbar,
    DbarD
};

enum CandidateType
{
    SelectedD = 0, // This particle is selected as a D
    SelectedDbar,  // This particle is selected as a Dbar
    TrueD,         // This particle is a true D
    TrueDbar       // This particle is a true Dbar
};

bool removeAmbiguous = false;
float mlCut = 0.5;

void loadInvMass()
{
    const char *fname = "~/MyMacros/Correlations_v2/triggered_data/AO2D_Data_Full2024.root";
    TFile *file = TFile::Open(fname, "read");

    const char *fname_MC = "~/MyMacros/Correlations_github/AO2D_MC_noAmbiguous.root";
    TFile *file_MC = TFile::Open(fname_MC, "read");

    TFile *fout = TFile::Open("fitFromTree.root", "RECREATE");

    // Create a new tree to store the merged data
    TTree *mergedTree = nullptr;
    TTree *mergedTreeMcInfo = nullptr;
    TTree *mergedTreeMl = nullptr;

    // Get the list of keys (directories) in the input file
    TList *listOfKeys = file_MC->GetListOfKeys();
    TIter next(listOfKeys);
    TObject *key;

    // Vector to store the trees from each directory
    vector<TTree *> trees;
    vector<TTree *> treesMcInfo;
    vector<TTree *> treesMl;

    // Loop over each directory
    while ((key = next())) {
        TDirectory *dir = dynamic_cast<TDirectory *>(file_MC->Get(key->GetName()));
        if (dir) {
            // Get the tree in the current directory
            TTree *tree = dynamic_cast<TTree *>(dir->Get("O2d0pair"));
            TTree *treeMcInfo = dynamic_cast<TTree *>(dir->Get("O2d0pairmcinfo"));
            TTree *treeMl = dynamic_cast<TTree *>(dir->Get("O2d0pairml"));
            if (tree && treeMcInfo && treeMl) {
                trees.push_back(tree);
                cout << "number of entries in directory" << key->GetName() << "  " << tree->GetEntries() << endl;
                treesMcInfo.push_back(treeMcInfo);
                treesMl.push_back(treeMl);
                cout << "Tree " << "O2D0pair" << " found in directory " << key->GetName() << endl;
                cout << "Tree " << "O2D0pairmcinfo" << " found in directory " << key->GetName() << endl;
                cout << "Tree " << "O2D0pairml" << " found in directory " << key->GetName() << endl;
            }
            else {
                cerr << "Error: Tree " << "O2D0pair, O2D0pairmcinfo or O2D0pairml" << " not found in directory " << key->GetName() << endl;
            }
        }
    }

    // Check if we found any trees
    if (trees.empty() || treesMcInfo.empty() || treesMl.empty()) {
        cerr << "Error: No trees found to merge" << endl;
        return;
    }

    // Merge the trees
    for (auto &tree : trees) {
        if (!mergedTree) {
            // Clone the structure of the first tree
            mergedTree = tree->CloneTree(0);
        }
        mergedTree->CopyEntries(tree);
    }

    // Merge the treesMcInfo
    for (auto &tree : treesMcInfo) {
        if (!mergedTreeMcInfo) {
            // Clone the structure of the first tree
            mergedTreeMcInfo = tree->CloneTree(0);
        }
        mergedTreeMcInfo->CopyEntries(tree);
    }

    // Merge the treesMcInfo
    for (auto &tree : treesMl) {
        if (!mergedTreeMl) {
            // Clone the structure of the first tree
            mergedTreeMl = tree->CloneTree(0);
        }
        mergedTreeMl->CopyEntries(tree);
    }

    TTree *mergedTreeTotal = nullptr;
    if (!mergedTreeTotal)
    {
        mergedTreeTotal = new TTree("mergedTreeTotal", "Merged Tree Total");

        // Variables from O2d0pair
        float ptCand1 = 0., ptCand2 = 0., yCand1 = 0., yCand2 = 0., phiCand1 = 0., phiCand2 = 0.;
        float mDCand1 = 0., mDCand2 = 0., mDbarCand1 = 0., mDbarCand2 = 0.;
        float bdtScoreD0Cand1 = 0., bdtScoreD0Cand2 = 0., bdtScoreD0barCand1 = 0., bdtScoreD0barCand2 = 0.;
        uint8_t typeCand1 = 0, typeCand2 = 0, typePair = 0;

        // Variables from O2d0pairmcinfo
        uint8_t origin1 = 0, origin2 = 0;

        // Set branch addresses for input trees
        mergedTree->SetBranchAddress("fPtCand1", &ptCand1);
        mergedTree->SetBranchAddress("fPtCand2", &ptCand2);
        mergedTree->SetBranchAddress("fYCand1", &yCand1);
        mergedTree->SetBranchAddress("fYCand2", &yCand2);
        mergedTree->SetBranchAddress("fPhiCand1", &phiCand1);
        mergedTree->SetBranchAddress("fPhiCand2", &phiCand2);
        mergedTree->SetBranchAddress("fMDCand1", &mDCand1);
        mergedTree->SetBranchAddress("fMDCand2", &mDCand2);
        mergedTree->SetBranchAddress("fMDbarCand1", &mDbarCand1);
        mergedTree->SetBranchAddress("fMDbarCand2", &mDbarCand2);
        mergedTree->SetBranchAddress("fCandidateType1", &typeCand1);
        mergedTree->SetBranchAddress("fCandidateType2", &typeCand2);
        mergedTree->SetBranchAddress("fPairType", &typePair);
        // MC info tree
        mergedTreeMcInfo->SetBranchAddress("fOrigin1", &origin1);
        mergedTreeMcInfo->SetBranchAddress("fOrigin2", &origin2);
        // ML tree
        mergedTreeMl->SetBranchAddress("fMlProbD0Cand1", &bdtScoreD0Cand1);
        mergedTreeMl->SetBranchAddress("fMlProbD0Cand2", &bdtScoreD0Cand2);
        mergedTreeMl->SetBranchAddress("fMlProbD0barCand1", &bdtScoreD0barCand1);
        mergedTreeMl->SetBranchAddress("fMlProbD0barCand2", &bdtScoreD0barCand2);

        // Create branches in the merged tree
        mergedTreeTotal->Branch("fPtCand1", &ptCand1);
        mergedTreeTotal->Branch("fPtCand2", &ptCand2);
        mergedTreeTotal->Branch("fYCand1", &yCand1);
        mergedTreeTotal->Branch("fYCand2", &yCand2);
        mergedTreeTotal->Branch("fPhiCand1", &phiCand1);
        mergedTreeTotal->Branch("fPhiCand2", &phiCand2);
        mergedTreeTotal->Branch("fMDCand1", &mDCand1);
        mergedTreeTotal->Branch("fMDCand2", &mDCand2);
        mergedTreeTotal->Branch("fMDbarCand1", &mDbarCand1);
        mergedTreeTotal->Branch("fMDbarCand2", &mDbarCand2);
        mergedTreeTotal->Branch("fCandidateType1", &typeCand1);
        mergedTreeTotal->Branch("fCandidateType2", &typeCand2);
        mergedTreeTotal->Branch("fPairType", &typePair);
        mergedTreeTotal->Branch("fOrigin1", &origin1);
        mergedTreeTotal->Branch("fOrigin2", &origin2);
        mergedTreeTotal->SetBranchAddress("fMlProbD0Cand1", &bdtScoreD0Cand1);
        mergedTreeTotal->SetBranchAddress("fMlProbD0Cand2", &bdtScoreD0Cand2);
        mergedTreeTotal->SetBranchAddress("fMlProbD0barCand1", &bdtScoreD0barCand1);
        mergedTreeTotal->SetBranchAddress("fMlProbD0barCand2", &bdtScoreD0barCand2);
    }

    // Loop over entries and fill the merged tree
    Long64_t nEntries = mergedTree->GetEntries();
    cout << "nentries: " << nEntries << endl;
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        mergedTree->GetEntry(i);
        mergedTreeMcInfo->GetEntry(i);
        mergedTreeTotal->Fill();
    }

    float ptCand1 = 0., ptCand2 = 0., yCand1 = 0., yCand2 = 0., phiCand1 = 0., phiCand2 = 0.;
    float mDCand1 = 0., mDCand2 = 0., mDbarCand1 = 0., mDbarCand2 = 0.;
    float bdtScoreD0Cand1 = 0., bdtScoreD0Cand2 = 0., bdtScoreD0barCand1 = 0., bdtScoreD0barCand2 = 0.;
    uint8_t typeCand1 = 0, typeCand2 = 0, typePair = 0;
    uint8_t origin1 = 0, origin2 = 0;

    mergedTreeTotal->SetBranchAddress("fPtCand1", &ptCand1);
    mergedTreeTotal->SetBranchAddress("fPtCand2", &ptCand2);
    mergedTreeTotal->SetBranchAddress("fYCand1", &yCand1);
    mergedTreeTotal->SetBranchAddress("fYCand2", &yCand2);
    mergedTreeTotal->SetBranchAddress("fPhiCand1", &phiCand1);
    mergedTreeTotal->SetBranchAddress("fPhiCand2", &phiCand2);
    mergedTreeTotal->SetBranchAddress("fMDCand1", &mDCand1);
    mergedTreeTotal->SetBranchAddress("fMDCand2", &mDCand2);
    mergedTreeTotal->SetBranchAddress("fMDbarCand1", &mDbarCand1);
    mergedTreeTotal->SetBranchAddress("fMDbarCand2", &mDbarCand2);
    mergedTreeTotal->SetBranchAddress("fCandidateType1", &typeCand1);
    mergedTreeTotal->SetBranchAddress("fCandidateType2", &typeCand2);
    mergedTreeTotal->SetBranchAddress("fPairType", &typePair);
    mergedTreeTotal->SetBranchAddress("fOrigin1", &origin1);
    mergedTreeTotal->SetBranchAddress("fOrigin2", &origin2);
    mergedTreeTotal->SetBranchAddress("fMlProbD0Cand1", &bdtScoreD0Cand1);
    mergedTreeTotal->SetBranchAddress("fMlProbD0Cand2", &bdtScoreD0Cand2);
    mergedTreeTotal->SetBranchAddress("fMlProbD0barCand1", &bdtScoreD0barCand1);
    mergedTreeTotal->SetBranchAddress("fMlProbD0barCand2", &bdtScoreD0barCand2);

    // Counter
    int nentries = mergedTreeTotal->GetEntries();
    int nDD = 0, nDbarDbar = 0, nDDbar = 0, nDbarD = 0, nDDbarAll = 0;

    // Create a histogram
    gROOT->cd();                       // just a precaution
    delete gROOT->FindObject("hLego"); // prevent "memory leak"
    TH2F *hLego = new TH2F("hLego", "Invariant Mass Distributions", 100, 1.7, 2.05, 100, 1.7, 2.05);
    TH1F *hPrompt = new TH1F("hPrompt", "Mc Prompt Invariant Mass Distributions", 100, 1.7, 2.05);
    TH1F *hNonPrompt = new TH1F("hNonPrompt", "Mc Non-prompt Invariant Mass Distributions", 100, 1.7, 2.05);
    TH1F *hReflection = new TH1F("hReflection", "Reflection Invariant Mass Distributions", 100, 1.7, 2.05);
    mergedTreeTotal->Draw("(fMDbarCand2) : (fMDbarCand1) >> hLego", "", "colz");

    float counterSelD = 0.;
    float counterTrueD = 0.;
    float counterSelTrueD = 0.;
    float counterSelDbar = 0.;
    float counterTrueDbar = 0.;
    float counterSelTrueDbar = 0.;
    for (int i = 0; i < nentries; i++)
    {
        mergedTreeTotal->GetEntry(i);

        // Select pT range
        if ((ptCand1 < 1.0 || ptCand2 < 1.0) || (ptCand1 > 24.0 || ptCand2 > 24.0)) {
            continue;
        }

        // Cut to remove ambiguous candidates
        if (removeAmbiguous && ((TESTBIT(typeCand1, SelectedD) && TESTBIT(typeCand1, SelectedDbar)))) {
            continue;
        }


        // Select signal events
        if ((TESTBIT(typeCand1, SelectedD) && TESTBIT(typeCand1, TrueD))) {
            // Cut on the ML score
            if (bdtScoreD0Cand1 < mlCut) { // Dummy cut
                continue;
            }
            if (origin1 == 1) {
                hPrompt->Fill(mDCand1);
            }
            if (origin1 == 2) {
                hNonPrompt->Fill(mDCand1);
            }
        }
        if (TESTBIT(typeCand1, SelectedDbar) && TESTBIT(typeCand1, TrueDbar))
        {
            // Cut on the ML score
            if (bdtScoreD0barCand1 < mlCut) { // Dummy cut -- same as the previous if
                continue;
            }
            if (origin1 == 1) {
                hPrompt->Fill(mDbarCand1);
            }
            if (origin1 == 2) {
                hNonPrompt->Fill(mDbarCand1);
            }
        }
        // Select reflected events
        if (!removeAmbiguous) {
            if (TESTBIT(typeCand1, SelectedDbar) && TESTBIT(typeCand1, TrueD)) {
                // Cut on the ML score
                if (bdtScoreD0barCand1 < mlCut) { // Dummy cut -- same as the previous if
                    continue;
                }
                if ((origin1 == 1) || (origin1 == 2))
                    hReflection->Fill(mDbarCand1);
            }
            if (TESTBIT(typeCand1, SelectedD) && TESTBIT(typeCand1, TrueDbar)) {
                // Cut on the ML score
                if (bdtScoreD0Cand1 < mlCut) { // Dummy cut -- same as the previous if
                    continue;
                }
                if ((origin1 == 1) || (origin1 == 2))
                    hReflection->Fill(mDCand1);
            }
        }

        // Counters
        if ((TESTBIT(typeCand1, SelectedDbar) && TESTBIT(typeCand1, TrueDbar))){
            counterSelTrueDbar++;
        }
        if ((TESTBIT(typeCand1, SelectedDbar))){
            counterSelDbar++;
        }
        if ((TESTBIT(typeCand1, TrueDbar))){
            counterTrueDbar++;
        }

        if ((TESTBIT(typeCand1, SelectedD) && TESTBIT(typeCand1, TrueD))){
            counterSelTrueD++;
        }
        if ((TESTBIT(typeCand1, SelectedD))){
            counterSelD++;
        }
        if ((TESTBIT(typeCand1, TrueD))){
            counterTrueD++;
        }
    }
    cout << "CounterSel D: " << counterSelD << "; CounterTrue D: " << counterTrueD << "; Counter SelTrue D: " << counterSelTrueD << std::endl;
    cout << "CounterSel Dbar: " << counterSelDbar << "; CounterTrue Dbar: " << counterTrueDbar << "; Counter SelTrue Dbar: " << counterSelTrueDbar << std::endl;

    TKey *keyData = (TKey*)file->GetListOfKeys()->At(0);
    TDirectoryFile *dirData = (TDirectoryFile*)file->Get(keyData->GetName());

    TTree *mergedTreeTotalData = new TTree("mergedTreeTotalData", "Merged Tree Total data");

    if (!dirData)
    {
        cout << ">> ERROR Directory not well readout" << endl;
        return;
    }
    TTree *treeData = (TTree *)dirData->Get("O2d0pair");
    if (!treeData)
    {
        cout << ">> ERROR Tree not well readout" << endl;
        return;
    }
    TTree *treeDataMl = (TTree *)dirData->Get("O2d0pairml");
    if (!treeDataMl)
    {
        cout << ">> ERROR data ML Tree not well readout" << endl;
        return;
    }
    float ptCand1Data = 0., ptCand2Data = 0., yCand1Data = 0., yCand2Data = 0., phiCand1Data = 0., phiCand2Data = 0.;
    float mDCand1Data = 0., mDCand2Data = 0., mDbarCand1Data = 0., mDbarCand2Data = 0.;
    float bdtScoreD0Cand1Data = 0., bdtScoreD0Cand2Data = 0., bdtScoreD0barCand1Data = 0., bdtScoreD0barCand2Data = 0.;
    uint8_t typeCand1Data = 0, typeCand2Data = 0, typePairData = 0;

    treeData->SetBranchAddress("fPtCand1", &ptCand1Data);
    treeData->SetBranchAddress("fPtCand2", &ptCand2Data);
    treeData->SetBranchAddress("fYCand1", &yCand1Data);
    treeData->SetBranchAddress("fYCand2", &yCand2Data);
    treeData->SetBranchAddress("fPhiCand1", &phiCand1Data);
    treeData->SetBranchAddress("fPhiCand2", &phiCand2Data);
    treeData->SetBranchAddress("fMDCand1", &mDCand1Data);
    treeData->SetBranchAddress("fMDCand2", &mDCand2Data);
    treeData->SetBranchAddress("fMDbarCand1", &mDbarCand1Data);
    treeData->SetBranchAddress("fMDbarCand2", &mDbarCand2Data);
    treeData->SetBranchAddress("fCandidateType1", &typeCand1Data);
    treeData->SetBranchAddress("fCandidateType2", &typeCand2Data);
    treeData->SetBranchAddress("fPairType", &typePairData);

    treeDataMl->SetBranchAddress("fMlProbD0Cand1", &bdtScoreD0Cand1);
    treeDataMl->SetBranchAddress("fMlProbD0Cand2", &bdtScoreD0Cand2);
    treeDataMl->SetBranchAddress("fMlProbD0barCand1", &bdtScoreD0barCand1);
    treeDataMl->SetBranchAddress("fMlProbD0barCand2", &bdtScoreD0barCand2);

    mergedTreeTotalData->SetBranchAddress("fPtCand1", &ptCand1Data);
    mergedTreeTotalData->SetBranchAddress("fPtCand2", &ptCand2Data);
    mergedTreeTotalData->SetBranchAddress("fYCand1", &yCand1Data);
    mergedTreeTotalData->SetBranchAddress("fYCand2", &yCand2Data);
    mergedTreeTotalData->SetBranchAddress("fPhiCand1", &phiCand1Data);
    mergedTreeTotalData->SetBranchAddress("fPhiCand2", &phiCand2Data);
    mergedTreeTotalData->SetBranchAddress("fMDCand1", &mDCand1Data);
    mergedTreeTotalData->SetBranchAddress("fMDCand2", &mDCand2Data);
    mergedTreeTotalData->SetBranchAddress("fMDbarCand1", &mDbarCand1Data);
    mergedTreeTotalData->SetBranchAddress("fMDbarCand2", &mDbarCand2Data);
    mergedTreeTotalData->SetBranchAddress("fCandidateType1", &typeCand1Data);
    mergedTreeTotalData->SetBranchAddress("fCandidateType2", &typeCand2Data);
    mergedTreeTotalData->SetBranchAddress("fPairType", &typePairData);
    mergedTreeTotalData->SetBranchAddress("fMlProbD0Cand1", &bdtScoreD0Cand1Data);
    mergedTreeTotalData->SetBranchAddress("fMlProbD0Cand2", &bdtScoreD0Cand2Data);
    mergedTreeTotalData->SetBranchAddress("fMlProbD0barCand1", &bdtScoreD0barCand1Data);
    mergedTreeTotalData->SetBranchAddress("fMlProbD0barCand2", &bdtScoreD0barCand2Data);

    // Loop over entries and fill the merged tree
    Long64_t nEntriesData = mergedTreeTotalData->GetEntries();
    for (Long64_t i = 0; i < nEntriesData; ++i)
    {
        treeData->GetEntry(i);
        treeDataMl->GetEntry(i);
        mergedTreeTotalData->Fill();
    }

    cout << "Number of entries in data: " << nEntriesData << endl;
    TH1F *hData = new TH1F("hData", "Data Invariant Mass Distributions", 100, 1.7, 2.05);

    for (int i = 0; i < nEntriesData; i++)
    {
        mergedTreeTotalData->GetEntry(i);

        // Select pT range
        if ((ptCand1Data < 1.0 || ptCand2Data < 1.0) || (ptCand1Data > 50.0 || ptCand2Data > 50.0)) {
            continue;
        }
        if (removeAmbiguous && (TESTBIT(typeCand1Data, SelectedD) && TESTBIT(typeCand1Data, SelectedDbar))) {
            continue;
        }
        if (TESTBIT(typeCand1Data, SelectedD)) {
            // Cut on the ML score
            if (bdtScoreD0Cand1Data < mlCut) { // Dummy cut -- same as the previous if
                continue;
            }
            hData->Fill(mDCand1Data);
        }
        if (TESTBIT(typeCand1Data, SelectedDbar)) {
            // Cut on the ML score
            if (bdtScoreD0barCand1Data < mlCut) { // Dummy cut -- same as the previous if
                continue;
            }
            hData->Fill(mDbarCand1Data);
        }
    }

    hData->Draw();

    fout->cd();
    hData->Write();
    hPrompt->Write();
    hNonPrompt->Write();
    hReflection->Write();
    hLego->Write();
    fout->Close();
}

