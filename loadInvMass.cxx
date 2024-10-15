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

void loadInvMass()
{
    //const char *fname = "~/MyMacros/Correlations/AO2D/AO2D_highIR.root";
    //const char *fname = "~/MyMacros/Correlations_v2/triggered_data/AO2D_LHC24_pass1.root";
    const char *fname = "~/MyMacros/Correlations_v2/triggered_data/AO2D_LHC23_pass4.root";
    TFile *file = TFile::Open(fname, "read");
    //const char *fname_MC = "~/MyMacros/Correlations/AO2D/AO2D_23b_fixed.root";
    const char *fname_MC = "~/Downloads/AO2D.root";
    TFile *file_MC = TFile::Open(fname_MC, "read");

    //TFile *fout = TFile::Open("fitFromTree_classDevelopment_cut2GeV.root", "RECREATE");
    //TFile *fout = TFile::Open("fitFromTree_triggeredData_LHC24.root", "RECREATE");
    TFile *fout = TFile::Open("fitFromTree_triggeredData_LHC23_pass4.root", "RECREATE");

    // Create a new tree to store the merged data
    TTree *mergedTree = nullptr;
    TTree *mergedTreeMcInfo = nullptr;

    // Get the list of keys (directories) in the input file
    TList *listOfKeys = file_MC->GetListOfKeys();
    TIter next(listOfKeys);
    TObject *key;

    // Vector to store the trees from each directory
    vector<TTree *> trees;
    vector<TTree *> treesMcInfo;

    // Loop over each directory
    while ((key = next()))
    {
        TDirectory *dir = dynamic_cast<TDirectory *>(file_MC->Get(key->GetName()));
        if (dir)
        {
            // Get the tree in the current directory
            TTree *tree = dynamic_cast<TTree *>(dir->Get("O2d0pair"));
            TTree *treeMcInfo = dynamic_cast<TTree *>(dir->Get("O2d0pairmcinfo"));
            if (tree && treeMcInfo)
            {
                trees.push_back(tree);
                cout << "number of entries in directory" << key->GetName() << "  " << tree->GetEntries() << endl;
                treesMcInfo.push_back(treeMcInfo);
                cout << "Tree " << "O2D0pair" << " found in directory " << key->GetName() << endl;
                cout << "Tree " << "O2D0pairmcinfo" << " found in directory " << key->GetName() << endl;
            }
            else
            {
                cerr << "Error: Tree " << "O2D0pair or O2D0pairMcInfo" << " not found in directory " << key->GetName() << endl;
            }
        }
    }

    // Check if we found any trees
    if (trees.empty() || treesMcInfo.empty())
    {
        cerr << "Error: No trees found to merge" << endl;
        return;
    }

    // Merge the trees
    for (auto &tree : trees)
    {
        if (!mergedTree)
        {
            // Clone the structure of the first tree
            mergedTree = tree->CloneTree(0);
        }
        mergedTree->CopyEntries(tree);
    }

    // Merge the treesMcInfo
    for (auto &tree : treesMcInfo)
    {
        if (!mergedTreeMcInfo)
        {
            // Clone the structure of the first tree
            mergedTreeMcInfo = tree->CloneTree(0);
        }
        mergedTreeMcInfo->CopyEntries(tree);
    }

    TTree *mergedTreeTotal = nullptr;
    if (!mergedTreeTotal)
    {
        mergedTreeTotal = new TTree("mergedTreeTotal", "Merged Tree Total");

        // Variables from O2d0pair
        float ptCand1 = 0., ptCand2 = 0., yCand1 = 0., yCand2 = 0., mDCand1 = 0., mDCand2 = 0., mDbarCand1 = 0., mDbarCand2 = 0.;
        uint8_t typeCand1 = 0, typeCand2 = 0, typePair = 0;

        // Variables from O2d0pairmcinfo
        uint8_t origin1 = 0, origin2 = 0;

        // Set branch addresses for input trees
        mergedTree->SetBranchAddress("fPtCand1", &ptCand1);
        mergedTree->SetBranchAddress("fPtCand2", &ptCand2);
        mergedTree->SetBranchAddress("fYCand1", &yCand1);
        mergedTree->SetBranchAddress("fYCand2", &yCand2);
        mergedTree->SetBranchAddress("fMDCand1", &mDCand1);
        mergedTree->SetBranchAddress("fMDCand2", &mDCand2);
        mergedTree->SetBranchAddress("fMDbarCand1", &mDbarCand1);
        mergedTree->SetBranchAddress("fMDbarCand2", &mDbarCand2);
        mergedTree->SetBranchAddress("fCandidateType1", &typeCand1);
        mergedTree->SetBranchAddress("fCandidateType2", &typeCand2);
        mergedTree->SetBranchAddress("fPairType", &typePair);

        mergedTreeMcInfo->SetBranchAddress("fOrigin1", &origin1);
        mergedTreeMcInfo->SetBranchAddress("fOrigin2", &origin2);

        // Create branches in the merged tree
        mergedTreeTotal->Branch("fPtCand1", &ptCand1);
        mergedTreeTotal->Branch("fPtCand2", &ptCand2);
        mergedTreeTotal->Branch("fYCand1", &yCand1);
        mergedTreeTotal->Branch("fYCand2", &yCand2);
        mergedTreeTotal->Branch("fMDCand1", &mDCand1);
        mergedTreeTotal->Branch("fMDCand2", &mDCand2);
        mergedTreeTotal->Branch("fMDbarCand1", &mDbarCand1);
        mergedTreeTotal->Branch("fMDbarCand2", &mDbarCand2);
        mergedTreeTotal->Branch("fCandidateType1", &typeCand1);
        mergedTreeTotal->Branch("fCandidateType2", &typeCand2);
        mergedTreeTotal->Branch("fPairType", &typePair);
        mergedTreeTotal->Branch("fOrigin1", &origin1);
        mergedTreeTotal->Branch("fOrigin2", &origin2);
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

    float ptCand1 = 0., ptCand2 = 0., yCand1 = 0., yCand2 = 0., mDCand1 = 0., mDCand2 = 0., mDbarCand1 = 0., mDbarCand2 = 0.;
    uint8_t typeCand1 = 0, typeCand2 = 0, typePair = 0, origin1 = 0, origin2 = 0;
    mergedTreeTotal->SetBranchAddress("fPtCand1", &ptCand1);
    mergedTreeTotal->SetBranchAddress("fPtCand2", &ptCand2);
    mergedTreeTotal->SetBranchAddress("fYCand1", &yCand1);
    mergedTreeTotal->SetBranchAddress("fYCand2", &yCand2);
    mergedTreeTotal->SetBranchAddress("fMDCand1", &mDCand1);
    mergedTreeTotal->SetBranchAddress("fMDCand2", &mDCand2);
    mergedTreeTotal->SetBranchAddress("fMDbarCand1", &mDbarCand1);
    mergedTreeTotal->SetBranchAddress("fMDbarCand2", &mDbarCand2);
    mergedTreeTotal->SetBranchAddress("fCandidateType1", &typeCand1);
    mergedTreeTotal->SetBranchAddress("fCandidateType2", &typeCand2);
    mergedTreeTotal->SetBranchAddress("fPairType", &typePair);
    mergedTreeTotal->SetBranchAddress("fOrigin1", &origin1);
    mergedTreeTotal->SetBranchAddress("fOrigin2", &origin2);

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
        ptCand1 = 0.;
        ptCand2 = 0.;
        yCand1 = 0.;
        mDCand1 = 0.;
        mDCand2 = 0.;
        mDbarCand1 = 0.;
        mDbarCand2 = 0.;
        typeCand1 = 0;
        typeCand2 = 0;
        typePair = 0;
        origin1 = 0;
        origin2 = 0;
        mergedTreeTotal->GetEntry(i);

        // Select pT range
        if ((ptCand1 < 2.0 || ptCand2 < 2.0) || (ptCand1 > 24.0 || ptCand2 > 24.0))
        {
            continue;
        }

        // Select signal events
        if ((TESTBIT(typeCand2, SelectedD) && TESTBIT(typeCand2, TrueD)))
        {
            if (origin2 == 1)
            {
                hPrompt->Fill(mDCand2);
            }
            if (origin2 == 2)
            {
                hNonPrompt->Fill(mDCand2);
            }
        }
        if (TESTBIT(typeCand2, SelectedDbar) && TESTBIT(typeCand2, TrueDbar))
        {
            if (origin2 == 1)
            {
                hPrompt->Fill(mDbarCand2);
            }
            if (origin2 == 2)
            {
                hNonPrompt->Fill(mDbarCand2);
            }
        }
        // Select reflected events
        if ((TESTBIT(typeCand2, SelectedDbar) && TESTBIT(typeCand2, TrueDbar)))
        {
            counterSelTrueDbar++;
        }
        if ((TESTBIT(typeCand2, SelectedDbar)))
        {
            counterSelDbar++;
        }
        if ((TESTBIT(typeCand2, TrueDbar)))
        {
            counterTrueDbar++;
        }

        if ((TESTBIT(typeCand2, SelectedD) && TESTBIT(typeCand2, TrueD)))
        {
            counterSelTrueD++;
        }
        if ((TESTBIT(typeCand2, SelectedD)))
        {
            counterSelD++;
        }
        if ((TESTBIT(typeCand2, TrueD)))
        {
            counterTrueD++;
        }
        if (TESTBIT(typeCand2, SelectedDbar) && TESTBIT(typeCand2, TrueD))
        {
            if ((origin2 == 1) || (origin2 == 2))
                hReflection->Fill(mDbarCand2);
        }
        if (TESTBIT(typeCand2, SelectedD) && TESTBIT(typeCand2, TrueDbar))
        {
            if ((origin2 == 1) || (origin2 == 2))
                hReflection->Fill(mDCand2);
        }
    }
    cout << "CounterSel D: " << counterSelD << "; CounterTrue D: " << counterTrueD << "; Counter SelTrue D: " << counterSelTrueD << std::endl;
    cout << "CounterSel Dbar: " << counterSelDbar << "; CounterTrue Dbar: " << counterTrueDbar << "; Counter SelTrue Dbar: " << counterSelTrueDbar << std::endl;

    //TDirectoryFile *dirData = (TDirectoryFile *)file->Get("DF_2261906152687232"); // highIR
    //TDirectoryFile *dirData = (TDirectoryFile *)file->Get("DF_2363808317917888"); // LHC24_pass1
    TDirectoryFile *dirData = (TDirectoryFile *)file->Get("DF_2298103932356800"); // LHC23_pass4

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
    float ptCand1Data = 0., ptCand2Data = 0., yCand1Data = 0., yCand2Data = 0., mDCand1Data = 0., mDCand2Data = 0., mDbarCand1Data = 0., mDbarCand2Data = 0.;
    uint8_t typeCand1Data = 0, typeCand2Data = 0, typePairData = 0;
    treeData->SetBranchAddress("fPtCand1", &ptCand1Data);
    treeData->SetBranchAddress("fPtCand2", &ptCand2Data);
    treeData->SetBranchAddress("fYCand1", &yCand1Data);
    treeData->SetBranchAddress("fYCand2", &yCand2Data);
    treeData->SetBranchAddress("fMDCand1", &mDCand1Data);
    treeData->SetBranchAddress("fMDCand2", &mDCand2Data);
    treeData->SetBranchAddress("fMDbarCand1", &mDbarCand1Data);
    treeData->SetBranchAddress("fMDbarCand2", &mDbarCand2Data);
    treeData->SetBranchAddress("fCandidateType1", &typeCand1Data);
    treeData->SetBranchAddress("fCandidateType2", &typeCand2Data);
    treeData->SetBranchAddress("fPairType", &typePairData);

    Long64_t nEntriesData = treeData->GetEntries();
    cout << "Number of entries in data: " << nEntriesData << endl;
    TH1F *hData = new TH1F("hData", "Data Invariant Mass Distributions", 100, 1.7, 2.05);

    for (int i = 0; i < nEntriesData; i++)
    {
        ptCand1Data = 0.;
        ptCand2Data = 0.;
        yCand1Data = 0.;
        mDCand1Data = 0.;
        mDCand2Data = 0.;
        mDbarCand1Data = 0.;
        mDbarCand2Data = 0.;
        typeCand1Data = 0;
        typeCand2Data = 0;
        typePairData = 0;
        treeData->GetEntry(i);

        // Select pT range
        if ((ptCand1Data < 1.0 || ptCand2Data < 1.0) || (ptCand1Data > 50.0 || ptCand2Data > 50.0))
        {
            continue;
        }
        if (TESTBIT(typeCand2Data, SelectedD)) {
            hData->Fill(mDCand2Data);
        }
        if (TESTBIT(typeCand2Data, SelectedDbar)) {
          hData->Fill(mDbarCand2Data);
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