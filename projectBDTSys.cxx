#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "Riostream.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TSystem.h"
#include "nlohmann/json.hpp"

using namespace std;
using json = nlohmann::json;

bool verbose;
void DrawMassHist(TH1F *h, TString title, TString xAxisLabel, TFile *f, TString dirName);
void Draw2DHist(TH2F *h, TString title, TString xAxisLabel, TString yAxisLabel, TFile *f, TString dirName);
void fillPtVsYHist(TH2F *hAdded, TH2F *h);
void fillAddedMassHist(TH1F *hAdded, TH1F *h);

// SPARSE AXIS:
// 0:Score Bkg
// 1:Score Signal
// 2:Inv. Mass
// 3:Pt
// 4:Y
// 5:PV contrib
// 6:Origin
// 7:CandType

enum candType {
    SelectedD =0,
    SelectedDbar,
    TrueD,
    TrueDbar
};

enum origin {
  Default = 0,
  Prompt,
  NonPrompt
};

TH3F *get3DHisto(THnSparseF* sparse, bool isMc, bool isSignal, std::string histoname, int const& origin);

void projectBDTSys()
{
  // Load JSON file
  ifstream jsonFile("config_bdtSys.json");
  if (!jsonFile.is_open())
  {
    cerr << "Error: Could not open JSON file." << endl;
    return;
  }

  std::string jsonString;
  jsonFile.seekg(0, std::ios::end);
  jsonString.reserve(jsonFile.tellg());
  jsonFile.seekg(0, std::ios::beg);
  jsonString.assign((std::istreambuf_iterator<char>(jsonFile)), std::istreambuf_iterator<char>());

  // Parse the JSON data.
  json jsonData = json::parse(jsonString);
  verbose = jsonData["verbose"];
  TString outputName = jsonData["outputfile"];
  TFile *outputFile = new TFile(outputName, "RECREATE");

  // Load AnalysisResults files
  TString filenameData = jsonData["filenameData"];
  TString filenamePrompt = jsonData["filenameMcPrompt"];
  TString filenameNonPrompt = jsonData["filenameMcNonPrompt"];

  TFile *fData = TFile::Open(filenameData);
  TFile *fPrompt = TFile::Open(filenamePrompt);
  TFile *fNonPrompt = TFile::Open(filenameNonPrompt);

  if (!fData) {
    cerr << "ERROR: file not found" << endl;
    return;
  }

  TString dirName = "hf-correlator-d-meson-pairs";
  TString sparseName = "hnDMesonMl";

  TDirectory *dirData = (TDirectory *)fData->Get(dirName);
  TDirectory *dirPrompt = (TDirectory *)fPrompt->Get(dirName);
  TDirectory *dirNonPrompt = (TDirectory *)fNonPrompt->Get(dirName);
  // Load sparses
  THnSparseF *sparseData = (THnSparseF *)dirData->Get(sparseName);
  THnSparseF *sparseMcPrompt = (THnSparseF *)dirPrompt->Get(sparseName);
  THnSparseF *sparseMcNonPrompt = (THnSparseF *)dirNonPrompt->Get(sparseName);

  if (!sparseData)
  {
    cerr << "ERROR: sparse data not found \n";
    return;
  }

  // GetOrigin
  sparseMcNonPrompt->GetAxis(6)->SetRangeUser(1.5, 2.5); // Non-Prompt
  cout << "Before setting range: " << sparseMcPrompt->GetEntries() << endl;
  sparseMcPrompt->GetAxis(6)->SetRangeUser(0.5, 1.5);
  // Project onto a 1D histogram (e.g., along the axis of interest)
  TH1D *hProjection = (TH1D *)sparseMcPrompt->Projection(6);
  cout << "Entries in selected range: " << hProjection->Integral() << endl;

  if (verbose)
  {
    cout << "Prompt and non-prompt axis selected" << endl;
  }

  // Get pT bins
  std::vector<int> nbinptmin = {}, nbinptmax = {};
  std::vector<float> nptmax = jsonData["nPtMax"];
  std::vector<float> nptmin = jsonData["nPtMin"];
  for (float i : nptmin) {
    nbinptmin.push_back(sparseMcPrompt->GetAxis(3)->FindBin(i));
    if (verbose) {
      cout << "nBinsPtMin: " << nbinptmin.back() << endl;
    }
  }

  for (float i : nptmax) {
    nbinptmax.push_back(sparseMcPrompt->GetAxis(3)->FindBin(i) - 1);
    if (verbose) {
      cout << "nBinsPtMax: " << nbinptmax.back() << endl;
    }
  }
  int const nPtBins = nbinptmin.size();
  if (verbose) {
    std::cout << "Bin number: " << nPtBins << std::endl;
  }

  TH3F *hData = get3DHisto(sparseData, false, false, "h3Ddata", Default);
  TH3F *hPromptSignal = get3DHisto(sparseMcPrompt, true, true, "h3DPromptSignal", Prompt);
  TH3F *hNonPromptSignal = get3DHisto(sparseMcNonPrompt, true, true, "h3DNonPromptSignal", NonPrompt);
  //TH3F *hPromptReflections = get3DHisto(sparseMcPrompt, true, false, "h3DPromptReflections");
  //TH3F *hNonPromptReflections = get3DHisto(sparseMcNonPrompt, true, false, "h3DNonPromptReflections");

  // Make cuts by projecting on the BDT and pT axes
  vector<vector<double>> vBkgScoreCuts = jsonData["bkgScoreCuts"];
  int const nBkgCuts = vBkgScoreCuts.size();
  int nPositions = vBkgScoreCuts[0].size();
  bool saveScorePlots = jsonData["saveScorePlots"];

  if (nBkgCuts != nPtBins) {
    cerr << "ERROR: must have a score cut per pT bin (nPtBins != nCuts) \n";
    return;
  }

  // Define histos to be filled
  // Data histos
  TH1F *hMassData[nPtBins][nPositions];
  TH1F *hBkgScore[nPtBins][nPositions];
  TH1F *hMassDataAdded[nPositions];
  // Mc prompt histos
  TH1F *hMassMcPrompt[nPtBins][nPositions];
  TH1F *hMassMcPromptAdded[nPositions];
  TH1F *hBkgScoreMcPrompt[nPtBins][nPositions];
  // Mc non-prompt histos
  TH1F *hMassMcNonPrompt[nPtBins][nPositions];
  TH1F *hMassMcNonPromptAdded[nPositions];
  TH1F *hBkgScoreMcNonPrompt[nPtBins][nPositions];

  for (int iPos = 0; iPos < nPositions; iPos++)
  {
    cout << "----------- POSITION " << iPos << " ----------- \n";
    outputFile->cd();

    auto posDir = "pos" + std::to_string(iPos) + "/";
    auto massDataDirName = posDir + "InvMassData/";
    auto massMcPromptDirName = posDir + "InvMassMcPrompt/";
    auto massMcNonPromptDirName = posDir + "InvMassMcNonPrompt/";

    auto bkgScoreDataDirName = posDir + "BkgScoreData/";

    gDirectory->mkdir(massDataDirName.c_str());
    gDirectory->mkdir(massMcPromptDirName.c_str());
    gDirectory->mkdir(massMcNonPromptDirName.c_str());

    if (saveScorePlots) {
      gDirectory->mkdir(bkgScoreDataDirName.c_str());
    }

    hMassDataAdded[iPos] = (TH1F *)hData->ProjectionX();
    hMassDataAdded[iPos]->SetNameTitle(Form("hMassDataAdded_pos%d", iPos), Form("hMassDataAdded_pos%d", iPos));
    hMassDataAdded[iPos]->Reset();
    hMassDataAdded[iPos]->SetNameTitle(Form("hMassDataAdded_pos%d", iPos), Form("hMassDataAdded_pos%d", iPos));

    hMassMcPromptAdded[iPos] = (TH1F *)hPromptSignal->ProjectionX();
    hMassMcPromptAdded[iPos]->SetNameTitle(Form("hMassMcPromptAdded_pos%d", iPos), Form("hMassMcPromptAdded_pos%d", iPos));
    hMassMcPromptAdded[iPos]->Reset();
    hMassMcPromptAdded[iPos]->SetNameTitle(Form("hMassMcPromptAdded_pos%d", iPos), Form("hMassMcPromptAdded_pos%d", iPos));

    hMassMcNonPromptAdded[iPos] = (TH1F *)hNonPromptSignal->ProjectionX();
    hMassMcNonPromptAdded[iPos]->SetNameTitle(Form("hMassMcNonPromptAdded_pos%d", iPos), Form("hMassMcNonPromptAdded_pos%d", iPos));
    hMassMcNonPromptAdded[iPos]->Reset();
    hMassMcNonPromptAdded[iPos]->SetNameTitle(Form("hMassMcNonPromptAdded_pos%d", iPos), Form("hMassMcNonPromptAdded_pos%d", iPos));

    int nBinsPtAdded = hMassDataAdded[iPos]->GetNbinsX();
    for (int ipt = 0; ipt < nPtBins; ipt++)
    {
      cout << "                PT RANGE: (" << nptmin[ipt] << ", " << nptmax[ipt] << ") \n";
      if (verbose) {
        cout << "Bin range: (" << nbinptmin[ipt] << ", " << nbinptmax[ipt] << ")" << endl;
      }

      // Set pT range
      hData->GetYaxis()->SetRange(nbinptmin[ipt], nbinptmax[ipt]);
      hPromptSignal->GetYaxis()->SetRange(nbinptmin[ipt], nbinptmax[ipt]);
      hNonPromptSignal->GetYaxis()->SetRange(nbinptmin[ipt], nbinptmax[ipt]);

      // Cut on BDT scores
      auto bkgCut = vBkgScoreCuts[ipt][iPos];

      hPromptSignal->GetZaxis()->SetRangeUser(bkgCut, 1.0);
      hNonPromptSignal->GetZaxis()->SetRangeUser(bkgCut, 1.0);
      hData->GetZaxis()->SetRange(hData->GetZaxis()->FindBin(bkgCut), hData->GetZaxis()->FindBin(1.0));
      cout << "bkg cut " << bkgCut << endl;

      // Get invariant-mass plots - Data
      auto massDataName = Form("InvMassData_%0.1f_%0.1f", nptmin[ipt], nptmax[ipt]);
      auto massDataTitle = "Invariant-mass Data, position " + std::to_string(iPos) + ", pt (" + std::to_string(nptmin[ipt]) + ", " + std::to_string(nptmax[ipt]) + ")";
      hMassData[ipt][iPos] = (TH1F *)hData->ProjectionX();
      hMassData[ipt][iPos]->Reset();
      hMassData[ipt][iPos]->SetName(massDataName);

      for (int xBin = 1; xBin <= hData->GetNbinsX(); xBin++) {
        for (int yBin = 1; yBin <= hData->GetNbinsY(); yBin++) {
            for (int zBin = 1; zBin <= hData->GetNbinsZ(); zBin++) {
    
                // Retrieve pT, mass, and BDT score values
                double massValue = hData->GetXaxis()->GetBinCenter(xBin);
                double pTValue = hData->GetYaxis()->GetBinCenter(yBin);
                double bdtValue = hData->GetZaxis()->GetBinCenter(zBin);
    
                // Apply BDT and pT cuts
                if (bdtValue < bkgCut || bdtValue > 1.0) {
                    continue;
                }
                if (pTValue < nptmin[ipt] || pTValue > nptmax[ipt]) {
                    continue;
                }
    
                // Get bin content from hData
                double binContent = hData->GetBinContent(xBin, yBin, zBin);
                if (binContent > 0) {
                    // Find the corresponding bin in the 1D mass histogram
                    int massBin = hMassData[ipt][iPos]->FindBin(massValue);
                    hMassData[ipt][iPos]->AddBinContent(massBin, binContent);
                }
            }
        }
    }

    hMassData[ipt][iPos]->SetEntries(hMassData[ipt][iPos]->Integral());

      DrawMassHist(hMassData[ipt][iPos], massDataTitle, "#it{M}_{inv}", outputFile, massDataDirName);
      fillAddedMassHist(hMassDataAdded[iPos], hMassData[ipt][iPos]);

      if (verbose) {
        cout << "Invariant-mass Data plots drawn and saved" << endl;
      }

      // Get invariant-mass plots - MC Prompt
      auto massMcPromptName = Form("InvMassMcPrompt_%0.1f_%0.1f", nptmin[ipt], nptmax[ipt]);
      auto massMcPromptTitle = "Invariant-mass MC Prompt, position " + std::to_string(iPos) + ", pt (" + std::to_string(nptmin[ipt]) + ", " + std::to_string(nptmax[ipt]) + ")";
      hMassMcPrompt[ipt][iPos] = (TH1F *)hPromptSignal->ProjectionX();
      hMassMcPrompt[ipt][iPos]->SetName(massMcPromptName);
      DrawMassHist(hMassMcPrompt[ipt][iPos], massMcPromptTitle, "#it{M}_{inv}", outputFile, massMcPromptDirName);
      fillAddedMassHist(hMassMcPromptAdded[iPos], hMassMcPrompt[ipt][iPos]);

      if (verbose) {
        cout << "Invariant-mass MC prompt plots drawn and saved" << endl;
      }

      // Get invariant-mass plots - MC Non Prompt
      auto massMcNonPromptName = Form("InvMassMcNonPrompt_%0.1f_%0.1f", nptmin[ipt], nptmax[ipt]);
      auto massMcNonPromptTitle = Form("Invariant-mass MC Non-Prompt, position %d, pT (%.0f, %.0f)", iPos, nptmin[ipt], nptmax[ipt]);
      hMassMcNonPrompt[ipt][iPos] = (TH1F *)hNonPromptSignal->ProjectionX();
      hMassMcNonPrompt[ipt][iPos]->SetName(massMcNonPromptName);
      DrawMassHist(hMassMcNonPrompt[ipt][iPos], massMcNonPromptTitle, "#it{M}_{inv}", outputFile, massMcNonPromptDirName);
      fillAddedMassHist(hMassMcNonPromptAdded[iPos], hMassMcNonPrompt[ipt][iPos]);

      if (verbose) {
        cout << "Invariant-mass MC non-prompt plots drawn and saved" << endl;
      }

      // Plot and save bkg and non-prompt BDT score plots
      if (saveScorePlots)
      {
        // Get bkg score plots - Data
        auto bkgScoreName = Form("BkgScoreData_%0.1f_%0.1f", nptmin[ipt], nptmax[ipt]);
        auto bkgScoreTitle = Form("Bkg Score Data, position %d, pT (%.0f, %.0f)", iPos, nptmin[ipt], nptmax[ipt]);
        hBkgScore[ipt][iPos] = (TH1F *)hData->ProjectionZ();
        hBkgScore[ipt][iPos]->SetName(bkgScoreName);
        DrawMassHist(hBkgScore[ipt][iPos], bkgScoreTitle, "Bkg Score", outputFile, bkgScoreDataDirName);

        if (verbose) {
          cout << "Bkg score Data plots drawn and saved" << endl;
        }
      }

      if (verbose) {
        cout << "Resetting pT range \n";
      }
      // Reset pT range
      hData->GetYaxis()->SetRange(1, sparseData->GetAxis(3)->GetNbins());
      hPromptSignal->GetYaxis()->SetRange(1, sparseData->GetAxis(3)->GetNbins());
      hNonPromptSignal->GetYaxis()->SetRange(1, sparseData->GetAxis(3)->GetNbins());

      if (verbose) {
        cout << "Resetting score cuts \n";
      }
      // Reset score cuts
      hData->GetZaxis()->SetRangeUser(0.0, 1.0);
      hPromptSignal->GetZaxis()->SetRangeUser(0.0, 1.0);
      hNonPromptSignal->GetZaxis()->SetRangeUser(0.0, 1.0);
    }

    auto massMcPromptAddedTitle = Form("Inv. mass MC Reco Prompt, position %.0d", iPos);
    auto massMcNonPromptAddedTitle = Form("Inv. mass MC Reco Non-prompt, position %.0d", iPos);
    auto massDataAddedTitle = Form("Inv. mass Data, position %.0d", iPos);

    DrawMassHist(hMassMcPromptAdded[iPos], massMcPromptAddedTitle, "#it{M}_{inv}", outputFile, posDir);
    DrawMassHist(hMassMcNonPromptAdded[iPos], massMcNonPromptAddedTitle, "#it{M}_{inv}", outputFile, posDir);
    DrawMassHist(hMassDataAdded[iPos], massDataAddedTitle, "#it{M}_{inv}", outputFile, posDir);
  }
  outputFile->Close();
}


TH3F *get3DHisto(THnSparseF* sparse, bool isMc, bool isSignal, std::string histoname, int const& origin) {

    // Define histogram: TH3F(name, title, bins in mass, pT, BDT)
    int nbinsMass = sparse->GetAxis(2)->GetNbins();
    double massMin = sparse->GetAxis(2)->GetXmin();
    double massMax = sparse->GetAxis(2)->GetXmax();

    int nbinsPt = sparse->GetAxis(3)->GetNbins();
    double ptMin = sparse->GetAxis(3)->GetXmin();
    double ptMax = sparse->GetAxis(3)->GetXmax();

    int nbinsBDT = sparse->GetAxis(1)->GetNbins();
    double bdtMin = sparse->GetAxis(1)->GetXmin();
    double bdtMax = sparse->GetAxis(1)->GetXmax();

    TH3F* h3D = new TH3F(histoname.c_str(), "Invariant Mass vs pT vs BDT Score",
                      nbinsMass, massMin, massMax,
                      nbinsPt, ptMin, ptMax,
                      nbinsBDT, bdtMin, bdtMax);

    // Get relevant axes
    auto axisBitmap = sparse->GetAxis(7);
    auto axisOrigin = sparse->GetAxis(6);
    auto axisMass = sparse->GetAxis(2);
    auto axisPt = sparse->GetAxis(3);
    auto axisBDT = sparse->GetAxis(1);

    // Create an iterator over the THnSparse
    // Cut for Data
    Long64_t nEntries = sparse->GetNbins();
    std::vector<int> coord(sparse->GetNdimensions(), 0);

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        // Get the global bin index from the entry
        sparse->GetBinContent(entry, coord.data());

        // Retrieve the actual bitmap value for this entry
        int bitmapValue = axisBitmap->GetBinCenter(coord[7]);
        int originValue = axisOrigin->GetBinCenter(coord[6]);
        
        // cut on the origin
        if (origin == Prompt && originValue != Prompt) {
          continue;
        }

        if (origin == NonPrompt && originValue != NonPrompt) {
          continue;
        }

        if (origin == Default && originValue != Default) {
          continue;
        }

        // Cut to remove ambiguous
        if ((TESTBIT(bitmapValue, SelectedD) && TESTBIT(bitmapValue, SelectedDbar))) {
            continue; // Skip this entry if it doesn't meet selection
        }
        // Cuts on MC
        if (isMc) {
            if (isSignal) {
                // Cut to remove non-signal candidates (MC only)
                if (!((TESTBIT(bitmapValue, SelectedD) && TESTBIT(bitmapValue, TrueD)) ||
                    (TESTBIT(bitmapValue, SelectedDbar) && TESTBIT(bitmapValue, TrueDbar)))) {
                    continue; // Skip this entry if it doesn't meet selection
                }
            } else {
                // Cut to remove non-reflected candidates (MC only)
                if (!((TESTBIT(bitmapValue, SelectedD) && TESTBIT(bitmapValue, TrueDbar)) ||
                    (TESTBIT(bitmapValue, SelectedDbar) && TESTBIT(bitmapValue, TrueD)))) {
                    continue; // Skip this entry if it doesn't meet selection
                }
            }
        }

        // Retrieve values for mass, pT, and BDT score
        double massValue = axisMass->GetBinCenter(coord[2]);
        double ptValue = axisPt->GetBinCenter(coord[3]);
        double bdtScore = axisBDT->GetBinCenter(coord[1]);

        // Get the bin content for this entry
        double binContent = sparse->GetBinContent(entry);

        // Fill TH3F histogram
        if (binContent > 0) {
            h3D->Fill(massValue, ptValue, bdtScore, binContent);
        }
    }

    return h3D;

}

// Función para pintar un único histograma y salvarlo
void DrawMassHist(TH1F *h, TString title, TString xAxisLabel, TFile *f, TString dirName)
{
  h->GetXaxis()->SetLabelFont(63); // font in pixels
  h->GetXaxis()->SetLabelSize(16); // in pixels
  h->GetYaxis()->SetLabelFont(63); // font in pixels
  h->GetYaxis()->SetLabelSize(16); // in pixels

  h->SetTitle(title);
  h->GetXaxis()->SetTitle(xAxisLabel);
  h->SetMarkerColor(kBlue - 3);
  h->SetLineColor(kBlue - 3);
  h->SetMarkerStyle(20);
  h->Draw("EP");

  f->cd(dirName);
  h->Write();
}

void fillAddedMassHist(TH1F *hAdded, TH1F *h)
{
  for (int binx = 1; binx <= hAdded->GetNbinsX(); binx++)
  {
    double mass = h->GetBinContent(binx);
    double massError = h->GetBinError(binx);

    double initMass = hAdded->GetBinContent(binx);
    double initMassError = hAdded->GetBinError(binx);

    hAdded->SetBinContent(binx, mass + initMass);
    hAdded->SetBinError(binx, massError + initMassError);
  }
}