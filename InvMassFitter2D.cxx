// Implementation of InvMassFitter2D class
// Author: Andrea Tavira García, IJCLab (Orsay, France)

#include "InvMassFitter2D.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "RooPlot.h"
#include "TPaveText.h"
#include "RooFitResult.h"
#include "TStyle.h"
#include "TBox.h"
#include "RooCurve.h"
#include "TLegend.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include <cstring> // for strcmp
#include "Math/Vector4D.h"
#include "BifurcatedCB.h"
#include "TH2D.h"
#include "TLine.h"
#include "RooHist.h"

using namespace RooFit;
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

enum fitOrder {
  Prefit = 0,
  Rawfit,
  Corrfit,
  PrefitCorr
};

// Constructors
InvMassFitter2D::InvMassFitter2D() : _tree(nullptr), _pairType("OS"), nentries(0),
                                     _massMin(1.70), _massMax(2.05), _ptMin(0.), _ptMax(50.),
                                     rooPtCand1("fPtCand1", "pt of candidate 1", -36., 36.),
                                     rooPtCand2("fPtCand2", "pt of candidate 2", -36., 36.),
                                     rooMCand1("fMCand1", "invariant-mass of the first candidate", 1.70, 2.05),
                                     rooMCand2("fMCand2", "invariant-mass of the second candidate", 1.70, 2.05),
                                     rooYCand1("fYCand1", "y of candidate 1", -1., 1.),
                                     rooYCand2("fYCand2", "y of candidate 2", -1., 1.),
                                     rooPhiCand1("fPhiCand1", "phi of candidate 1", -36., 36.),
                                     rooPhiCand2("fPhiCand2", "phi of candidate 2", -36., 36.),
                                     rooPtPair("fPtPair", "pt of the pair", 0., 50.0),
                                     _efficiencyMap(0x0),
                                     _removeAmbiguous(false),
                                     _mean(0x0), _sigma(0x0), _meanRefl(0x0), _sigmaRefl(0x0), _tau(0x0), _fracRefl(0x0),
                                     _meanReflDoubleGaus(0x0), _sigmaReflDoubleGaus(0x0), _rawYield(0x0), _reflOverSgn(0),
                                     _workspace(0x0), _bkgPdfCand1(0x0), _sgnPdfCand1(0x0), _reflPdfCand1(0x0), _totPdfCand1(0x0),
                                     _bkgPdfCand2(0x0), _sgnPdfCand2(0x0), _reflPdfCand2(0x0), _totPdfCand2(0x0) {}
InvMassFitter2D::InvMassFitter2D(TTree *tree, const char *pairType) : _tree(tree), _pairType(pairType), nentries(0),
                                                                      _massMin(1.70), _massMax(2.05), _ptMin(0.), _ptMax(50.),
                                                                      rooPtCand1("fPtCand1", "pt of candidate 1", -36., 36.),
                                                                      rooPtCand2("fPtCand2", "pt of candidate 2", -36., 36.),
                                                                      rooMCand1("fMCand1", "invariant-mass of the first candidate", 1.70, 2.05),
                                                                      rooMCand2("fMCand2", "invariant-mass of the second candidate", 1.70, 2.05),
                                                                      rooYCand1("fYCand1", "y of candidate 1", -1., 1.),
                                                                      rooYCand2("fYCand2", "y of candidate 2", -1., 1.),
                                                                      rooPhiCand1("fPhiCand1", "phi of candidate 1", -36., 36.),
                                                                      rooPhiCand2("fPhiCand2", "phi of candidate 2", -36., 36.),
                                                                      rooPtPair("fPtPair", "pt of the pair", 0., 50.0),
                                                                      _efficiencyMap(0x0),
                                                                      _removeAmbiguous(false),
                                                                      _mean(0x0), _sigma(0x0), _meanRefl(0x0), _sigmaRefl(0x0), _tau(0x0), _fracRefl(0x0),
                                                                      _meanReflDoubleGaus(0x0), _sigmaReflDoubleGaus(0x0), _rawYield(0x0), _reflOverSgn(0),
                                                                      _workspace(0x0), _bkgPdfCand1(0x0), _sgnPdfCand1(0x0), _reflPdfCand1(0x0), _totPdfCand1(0x0),
                                                                      _bkgPdfCand2(0x0), _sgnPdfCand2(0x0), _reflPdfCand2(0x0), _totPdfCand2(0x0)
{
  if (_tree)
  {
    _loadTreeInfo();
  }
}

void InvMassFitter2D::_loadTreeInfo()
{
  if (!_tree)
  {
    cerr << "ERROR: tree not found!" << endl;
    return;
  }
  ////// PT Y PHI ESTÁN INVERTIDOS EN LOS TREES QUE USARON ML EN EL CORRELATOR EN DATOS!!!!! -- Sólo en el código del Preliminary
  _tree->SetBranchAddress("fPtCand1", &ptCand1);
  _tree->SetBranchAddress("fPtCand2", &ptCand2);
  _tree->SetBranchAddress("fYCand1", &yCand1);
  _tree->SetBranchAddress("fYCand2", &yCand2);
  _tree->SetBranchAddress("fPhiCand1", &phiCand1);
  _tree->SetBranchAddress("fPhiCand2", &phiCand2);
  _tree->SetBranchAddress("fMDCand1", &mDCand1);
  _tree->SetBranchAddress("fMDCand2", &mDCand2);
  _tree->SetBranchAddress("fMDbarCand1", &mDbarCand1);
  _tree->SetBranchAddress("fMDbarCand2", &mDbarCand2);
  _tree->SetBranchAddress("fCandidateType1", &typeCand1);
  _tree->SetBranchAddress("fCandidateType2", &typeCand2);
  _tree->SetBranchAddress("fPairType", &typePair);

  nentries = _tree->GetEntries();

  int nDD = 0, nDbarDbar = 0, nDDbar = 0, nDbarD = 0, nDDbarAll = 0;

  if (nentries == 0)
  {
    cerr << "ERROR: the tree contains no entries!" << endl;
    return;
  }

  for (int i = 0; i < nentries; i++)
  {
    _tree->GetEntry(i);

    if (TESTBIT(typePair, DD))
      nDD++;
    if (TESTBIT(typePair, DbarDbar))
      nDbarDbar++;
    if (TESTBIT(typePair, DDbar))
      nDDbar++;
    if (TESTBIT(typePair, DbarD))
      nDbarD++;
    if (TESTBIT(typePair, DDbar) | TESTBIT(typePair, DbarD))
      nDDbarAll++;
  }
  cout << " A total of " << nentries << " pairs analyzed " << endl;
  cout << "   nDD: " << nDD << ", nDbarDbar: " << nDbarDbar << "\n   nDDbar: " << nDDbar << ", nDbarD: " << nDbarD << ", nDDbarAll: " << nDDbarAll << endl;
}

void InvMassFitter2D::removeAmbiguous(bool remove=false) {
  _removeAmbiguous = remove;
}

void InvMassFitter2D::setPtLims(double const &ptMin, double const &ptMax)
{
  _ptMin = ptMin;
  _ptMax = ptMax;
}

void InvMassFitter2D::setPtPairLims(double const &ptMinPair, double const &ptMaxPair)
{
  _ptMinPair = ptMinPair;
  _ptMaxPair = ptMaxPair;
}

void InvMassFitter2D::setMassLims(double const &massMin, double const &massMax)
{
  _massMin = massMin;
  _massMax = massMax;
}

void InvMassFitter2D::setEfficiencyMap(TH2D *h)
{
  _efficiencyMap = (TH2D *)h->Clone();
}

void InvMassFitter2D::setSgnFunc(TString const &sgnFunc)
{
  _sgnFuncOption = sgnFunc;
}

void InvMassFitter2D::setBkgFunc(TString const &bkgFunc)
{
  _bkgFuncOption = bkgFunc;
}

void InvMassFitter2D::setReflFunc(TString const &reflFunc)
{
  _reflFuncOption = reflFunc;
}

void InvMassFitter2D::createDataset()
{
  RooCategory rooPairType("rooPairType", "pair type enum in RooFit");
  RooCategory rooCandType1("rooCandType1", "candidate 1 type enum in RooFit");
  RooCategory rooCandType2("rooCandType2", "candidate 2 type enum in RooFit");

  rooPairType.defineType("DD", DD);
  rooPairType.defineType("DbarDbar", DbarDbar);
  rooPairType.defineType("DDbar", DDbar);
  rooPairType.defineType("DbarD", DbarD);

  rooCandType1.defineType("SelectedD", SelectedD);
  rooCandType1.defineType("SelectedDbar", SelectedDbar);
  rooCandType1.defineType("TrueD", TrueD);
  rooCandType1.defineType("TrueDbar", TrueDbar);

  rooCandType2.defineType("SelectedD", SelectedD);
  rooCandType2.defineType("SelectedDbar", SelectedDbar);
  rooCandType2.defineType("TrueD", TrueD);
  rooCandType2.defineType("TrueDbar", TrueDbar);
  // Create a rooArgSet containing our variables of interest
  RooArgSet vars(rooPtCand1, rooPtCand2, rooMCand1, rooMCand2, rooYCand1, rooYCand2, rooPhiCand1, rooPhiCand2, rooPtPair); // Cand 1: D, Cand 2: Dbar (OS)

  RooRealVar weightCand1("weightCand1", "weights of cand 1", 1., 0., 100.);
  RooRealVar weightCand2("weightCand2", "weights of cand 2", 1., 0., 100.);
  RooArgSet weightedVars(rooPtCand1, rooPtCand2, rooMCand1, rooMCand2, rooYCand1, rooYCand2, rooPhiCand1, rooPhiCand2, rooPtPair, weightCand1, weightCand2);
  RooFormulaVar combinedWeight("combinedWeight", "combined weight", "weightCand1 * weightCand2", RooArgList(weightCand1, weightCand2));

  // Create an empty dataset with the variables and category
  RooDataSet data("data", "data", vars);
  RooDataSet weightedData("weightedData", "weightedData", weightedVars, RooFit::WeightVar(combinedWeight.GetName()));

  fillDataset(data, vars);
  fillDataset(weightedData, weightedVars);
}

void InvMassFitter2D::fillDataset(RooDataSet &data, RooArgSet &vars)
{

  RooRealVar *weightsCand1 = nullptr;
  RooRealVar *weightsCand2 = nullptr;

  double totalWeight = 0.; // Variable to store the total weight for normalization
  if (_efficiencyMap && vars.getSize() != 9)
  {
    weightsCand1 = dynamic_cast<RooRealVar *>(vars.find("weightCand1"));
    weightsCand2 = dynamic_cast<RooRealVar *>(vars.find("weightCand2"));
  }

  // Fill the dataset with info from the tree
  cout << "Number of tree entries " << _tree->GetEntries() << endl;
  ////// PT Y PHI ESTÁN INVERTIDOS EN LOS TREES QUE USARON ML EN EL CORRELATOR EN DATOS!!!!! -- Sólo en el código del Preliminary
  _tree->SetBranchAddress("fPtCand1", &ptCand1);
  _tree->SetBranchAddress("fPtCand2", &ptCand2);
  _tree->SetBranchAddress("fYCand1", &yCand1);
  _tree->SetBranchAddress("fYCand2", &yCand2);
  _tree->SetBranchAddress("fPhiCand1", &phiCand1);
  _tree->SetBranchAddress("fPhiCand2", &phiCand2);
  _tree->SetBranchAddress("fMDCand1", &mDCand1);
  _tree->SetBranchAddress("fMDCand2", &mDCand2);
  _tree->SetBranchAddress("fMDbarCand1", &mDbarCand1);
  _tree->SetBranchAddress("fMDbarCand2", &mDbarCand2);
  _tree->SetBranchAddress("fCandidateType1", &typeCand1);
  _tree->SetBranchAddress("fCandidateType2", &typeCand2);
  _tree->SetBranchAddress("fPairType", &typePair);

  int counterBeforeAmbiguous = 0;
  int counterAfterAmbiguous = 0;

  for (int i = 0; i < nentries; i++)
  {
    ptCand1 = 0.;
    ptCand2 = 0.;
    yCand1 = 0.;
    yCand2 = 0.;
    phiCand1 = 0.;
    phiCand2 = 0.;
    mDCand1 = 0.;
    mDCand2 = 0.;
    mDbarCand1 = 0.;
    mDbarCand2 = 0.;
    typeCand1 = 0;
    typeCand2 = 0;
    typePair = 0;
    _tree->GetEntry(i);

    // Select pT range using the pT of the pair
    if ((ptCand1 < _ptMin || ptCand2 < _ptMin) || (ptCand1 > _ptMax || ptCand2 > _ptMax)) {
      continue;
    }

    if ((yCand1 > 0.5 || yCand1 < -0.5) || (yCand2 > 0.5 || yCand2 < -0.5)) {
      continue;
    }

    counterBeforeAmbiguous++;
    // Add cuts to avoid ambiguous candidates
    if (_removeAmbiguous && (TESTBIT(typeCand1, SelectedD) && TESTBIT(typeCand1, SelectedDbar))) {
      continue;
    }
    if (_removeAmbiguous && (TESTBIT(typeCand2, SelectedD) && TESTBIT(typeCand2, SelectedDbar))) {
      continue;
    }
    counterAfterAmbiguous++;

    std::vector<std::string> matchedTypes;

    // Check all matching pair types
    if (TESTBIT(typePair, DD)) {
        if (!((mDCand1 < _massMin || mDCand1 > _massMax) || (mDCand2 < _massMin || mDCand2 > _massMax))) {
            matchedTypes.push_back("DD");
        }
    }
    if (TESTBIT(typePair, DbarDbar)) {
        if (!((mDbarCand1 < _massMin || mDbarCand1 > _massMax) || (mDbarCand2 < _massMin || mDbarCand2 > _massMax))) {
            matchedTypes.push_back("DbarDbar");
        }
    }
    if (TESTBIT(typePair, DDbar)) {
        if (!((mDCand1 < _massMin || mDCand1 > _massMax) || (mDbarCand2 < _massMin || mDbarCand2 > _massMax))) {
            matchedTypes.push_back("DDbar");
        }
    }
    if (TESTBIT(typePair, DbarD)) {
        if (!((mDbarCand1 < _massMin || mDbarCand1 > _massMax) || (mDCand2 < _massMin || mDCand2 > _massMax))) {
            matchedTypes.push_back("DbarD");
        }
    }

    // Process all matching pair types
    for (const auto &pairType : matchedTypes) {
      ROOT::Math::PxPyPzMVector vLorentzCand1; ROOT::Math::PxPyPzMVector vLorentzCand2;
      // Create Lorentz vectors
      if (pairType == "DD") {
        vLorentzCand1 = createLorentzVector(phiCand1, yCand1, ptCand1, mDCand1);
        vLorentzCand2 = createLorentzVector(phiCand2, yCand2, ptCand2, mDCand2);
      } else if (pairType == "DbarDbar") {
        vLorentzCand1 = createLorentzVector(phiCand1, yCand1, ptCand1, mDCand1);
        vLorentzCand2 = createLorentzVector(phiCand2, yCand2, ptCand2, mDbarCand2);
      } else if (pairType == "DDbar") {
        vLorentzCand1 = createLorentzVector(phiCand1, yCand1, ptCand1, mDbarCand1);
        vLorentzCand2 = createLorentzVector(phiCand2, yCand2, ptCand2, mDCand2);
      } else if (pairType == "DbarD") {
        vLorentzCand1 = createLorentzVector(phiCand1, yCand1, ptCand1, mDbarCand1);
        vLorentzCand2 = createLorentzVector(phiCand2, yCand2, ptCand2, mDbarCand2);
      }
      ROOT::Math::PxPyPzMVector vLorentzPair = vLorentzCand1 + vLorentzCand2;

      if (!(vLorentzPair.Pt() >= _ptMinPair && vLorentzPair.Pt() <= _ptMaxPair)) {
        continue;
      }

      if (pairType == "DD") {
        rooMCand1.setVal(mDCand1);
        rooMCand2.setVal(mDCand2);
      } else if (pairType == "DbarDbar") {
        rooMCand1.setVal(mDbarCand1);
        rooMCand2.setVal(mDbarCand2);
      } else if (pairType == "DDbar") {
        rooMCand1.setVal(mDCand1);
        rooMCand2.setVal(mDbarCand2);
      } else if (pairType == "DbarD") {
        rooMCand1.setVal(mDCand2);
        rooMCand2.setVal(mDbarCand1);
      }

      // Set kinematic variables
      rooPtCand1.setVal(ptCand1);
      rooPtCand2.setVal(ptCand2);
      rooYCand1.setVal(yCand1);
      rooYCand2.setVal(yCand2);
      rooPhiCand1.setVal(phiCand1);
      rooPhiCand2.setVal(phiCand2);
      rooPtPair.setVal(vLorentzPair.Pt());

      double weightCand1 = 1., weightCand2 = 1.;
      double combinedWeight = 1.;

      totalWeight += combinedWeight;

      if (_efficiencyMap) {
        weightCand1 = calculateWeights(yCand1, ptCand1);
        weightCand2 = calculateWeights(yCand2, ptCand2);
        combinedWeight = weightCand1 * weightCand2;
        if (vars.getSize() != 9) {
          weightsCand1->setVal(weightCand1);
          weightsCand2->setVal(weightCand2);
        }
      }

      if (strcmp(_pairType, "OS") == 0) {
        if (pairType == "DDbar" || pairType == "DbarD") {
          if (vars.getSize() != 9)
            data.add(vars, combinedWeight);
          else
            data.add(vars);
        }
      } else if (strcmp(_pairType, "LS") == 0) {
        if (pairType == "DD" || pairType == "DbarDbar") {
          if (vars.getSize() != 9)
            data.add(vars, combinedWeight);
          else
            data.add(vars);
        }
      } else {
        cerr << "ERROR: wrong pairType assigned. Please choose OS or LS" << endl;
        return;
      }
    }
  }

  // Import dataset after processing all pair types
  _workspace.import(data);
  cout << "data loaded" << endl;

  std::cout << "Counter before removing ambiguous: " << counterBeforeAmbiguous << endl;
  std::cout << "Counter after removing ambiguous: " << counterAfterAmbiguous << endl;
}

void InvMassFitter2D::set1DParameters(double const &reflOverSgn, double const &integratedEfficiency)
{
  _reflOverSgn = reflOverSgn;
  _integratedEfficiency = integratedEfficiency;
}

/// @brief Fill workspace with relevant functions for signal, bkg and reflections
/// Functions are triplicated to make the prefit, raw fit and efficiecy-corrected fit
void InvMassFitter2D::fillWorkspace(RooDataSet *dataset)
{
  const RooArgSet *vars = dataset->get();
  ////// PT Y PHI ESTÁN INVERTIDOS EN LOS TREES QUE USARON ML EN EL CORRELATOR EN DATOS!!!!! --- Sólo en el código del Preliminary
  RooRealVar *massCand1 = dynamic_cast<RooRealVar *>(vars->find("fMCand1"));
  RooRealVar *massCand2 = dynamic_cast<RooRealVar *>(vars->find("fMCand2"));
  RooRealVar *ptCand1 = dynamic_cast<RooRealVar *>(vars->find("fPtCand1"));
  RooRealVar *ptCand2 = dynamic_cast<RooRealVar *>(vars->find("fPtCand2"));
  RooRealVar *yCand1 = dynamic_cast<RooRealVar *>(vars->find("fYCand1"));
  RooRealVar *yCand2 = dynamic_cast<RooRealVar *>(vars->find("fYCand2"));

  massCand1->setRange(_massMin, _massMax);
  massCand2->setRange(_massMin, _massMax);

  /// | ------------------------------------------------------------------ |
  /// | ------------------- BACKGROUND FUNCTIONS ------------------------- |
  /// | ------------------------------------------------------------------ |
  cout << "loading functions" << endl;
  // bkg expo
  RooRealVar tauCand1Prefit("tauCand1Prefit", "tauCand1Prefit", -1., -3., 1.);
  RooRealVar tauCand2Prefit("tauCand2Prefit", "tauCand2Prefit", -1, -3., 1.);

  RooRealVar tauCand1("tauCand1", "tauCand1", -1., -3., 1.);
  RooRealVar tauCand2("tauCand2", "tauCand2", -1, -3., 1.);

  RooRealVar tauCand1Corr("tauCand1Corr", "tauCand1Corr", -1., -3., 1.);
  RooRealVar tauCand2Corr("tauCand2Corr", "tauCand2Corr", -1, -3., 1.);

  RooRealVar tauCand1CorrPrefit("tauCand1CorrPrefit", "tauCand1CorrPrefit", -1., -3., 1.);
  RooRealVar tauCand2CorrPrefit("tauCand2CorrPrefit", "tauCand2CorrPrefit", -1, -3., 1.);

  RooAbsPdf *bkgFuncExpoCand1Prefit = new RooExponential("bkgFuncExpoCand1Prefit", "background prefit function of candidate 1", *massCand1, tauCand1Prefit);
  RooAbsPdf *bkgFuncExpoCand2Prefit = new RooExponential("bkgFuncExpoCand2Prefit", "background prefit function of candidate 2", *massCand2, tauCand2Prefit);

  RooAbsPdf *bkgFuncExpoCand1 = new RooExponential("bkgFuncExpoCand1", "background fit function of candidate 1", *massCand1, tauCand1);
  RooAbsPdf *bkgFuncExpoCand2 = new RooExponential("bkgFuncExpoCand2", "background fit function of candidate 2", *massCand2, tauCand2);

  RooAbsPdf *bkgFuncExpoCand1Corr = new RooExponential("bkgFuncExpoCand1Corr", "Corrected background fit function of candidate 1", *massCand1, tauCand1Corr);
  RooAbsPdf *bkgFuncExpoCand2Corr = new RooExponential("bkgFuncExpoCand2Corr", "Corrected background fit function of candidate 2", *massCand2, tauCand2Corr);

  RooAbsPdf *bkgFuncExpoCand1CorrPrefit = new RooExponential("bkgFuncExpoCand1CorrPrefit", "Prefit Corrected background fit function of candidate 1", *massCand1, tauCand1CorrPrefit);
  RooAbsPdf *bkgFuncExpoCand2CorrPrefit = new RooExponential("bkgFuncExpoCand2CorrPrefit", "Prefit Corrected background fit function of candidate 2", *massCand2, tauCand2CorrPrefit);

  _workspace.import(*bkgFuncExpoCand1Prefit);
  _workspace.import(*bkgFuncExpoCand2Prefit);
  _workspace.import(*bkgFuncExpoCand1);
  _workspace.import(*bkgFuncExpoCand2);
  _workspace.import(*bkgFuncExpoCand1Corr);
  _workspace.import(*bkgFuncExpoCand2Corr);
  _workspace.import(*bkgFuncExpoCand1CorrPrefit);
  _workspace.import(*bkgFuncExpoCand2CorrPrefit);

  // bkg poly0: c1 (useful for MC)
  RooRealVar c1Prefit("c1Prefit", "Prefit Linear coefficient", 10, -10., 1e10);      // c1: constant coefficient
  RooRealVar c1Cand2Prefit("c1Cand2Prefit", "Prefit Linear coefficient", 10, -10., 1e10);      // c1: constant coefficient

  RooRealVar c1("c1", "Linear coefficient", -1, -10., 1e10);      // c1: constant coefficient
  RooRealVar c1Cand2("c1Cand2", "Linear coefficient", -1, -10., 1e10);      // c1: constant coefficient

  RooRealVar c1Corr("c1Corr", "Corrected Linear coefficient", 1, -10., 1e10);      // c1: constant coefficient
  RooRealVar c1Cand2Corr("c1Cand2Corr", "Corrected Linear coefficient", 1, -10., 1e10);      // c1: constant coefficient

  RooRealVar c1CorrPrefit("c1CorrPrefit", "Corrected Prefit Linear coefficient", 1, -10., 1e10);      // c1: constant coefficient
  RooRealVar c1Cand2CorrPrefit("c1Cand2CorrPrefit", "Corrected Prefit Linear coefficient", 1, -10., 1e10);      // c1: constant coefficient

  RooAbsPdf *bkgFuncPoly0Cand1Prefit = new RooPolynomial("bkgFuncPoly0Cand1Prefit", "background prefit function of candidate 1", *massCand1, RooArgList(c1Prefit));
  RooAbsPdf *bkgFuncPoly0Cand2Prefit = new RooPolynomial("bkgFuncPoly0Cand2Prefit", "background prefit function of candidate 2", *massCand2, RooArgList(c1Cand2Prefit));

  RooAbsPdf *bkgFuncPoly0Cand1 = new RooPolynomial("bkgFuncPoly0Cand1", "background fit function of candidate 1", *massCand1, RooArgList(c1));
  RooAbsPdf *bkgFuncPoly0Cand2 = new RooPolynomial("bkgFuncPoly0Cand2", "background fit function of candidate 2", *massCand2, RooArgList(c1Cand2));

  RooAbsPdf *bkgFuncPoly0Cand1Corr = new RooPolynomial("bkgFuncPoly0Cand1Corr", "Corrected background fit function of candidate 1", *massCand1, RooArgList(c1Corr));
  RooAbsPdf *bkgFuncPoly0Cand2Corr = new RooPolynomial("bkgFuncPoly0Cand2Corr", "Corrected background fit function of candidate 2", *massCand2, RooArgList(c1Cand2Corr));

  RooAbsPdf *bkgFuncPoly0Cand1CorrPrefit = new RooPolynomial("bkgFuncPoly0Cand1CorrPrefit", "Prefit Corrected background fit function of candidate 1", *massCand1, RooArgList(c1CorrPrefit));
  RooAbsPdf *bkgFuncPoly0Cand2CorrPrefit = new RooPolynomial("bkgFuncPoly0Cand2CorrPrefit", "Prefit Corrected background fit function of candidate 2", *massCand2, RooArgList(c1Cand2CorrPrefit));

  _workspace.import(*bkgFuncPoly0Cand1Prefit);
  _workspace.import(*bkgFuncPoly0Cand2Prefit);
  _workspace.import(*bkgFuncPoly0Cand1);
  _workspace.import(*bkgFuncPoly0Cand2);
  _workspace.import(*bkgFuncPoly0Cand1Corr);
  _workspace.import(*bkgFuncPoly0Cand2Corr);
  _workspace.import(*bkgFuncPoly0Cand1CorrPrefit);
  _workspace.import(*bkgFuncPoly0Cand2CorrPrefit);

  // bkg poly1: c1 + c2*x
  RooRealVar c2Prefit("c2Prefit", "Prefit Linear coefficient", -3.15, -10, 1.); // c2: linear coefficient
  RooRealVar c2Cand2Prefit("c2Cand2Prefit", "Prefit Linear coefficient", -3.15, -10, 1.); // c2: linear coefficient

  RooRealVar c2("c2", "Linear coefficient", -3.15, -10, 1.); // c2: linear coefficient
  RooRealVar c2Cand2("c2Cand2", "Linear coefficient", -3.15, -10, 1.); // c2: linear coefficient

  RooRealVar c2Corr("c2Corr", "Corrected Linear coefficient", -3.15, -10, 1.); // c2: linear coefficient
  RooRealVar c2Cand2Corr("c2Cand2Corr", "Corrected Linear coefficient", -3.15, -10, 1.); // c2: linear coefficient

  RooRealVar c2CorrPrefit("c2CorrPrefit", "Corrected Prefit Linear coefficient", -3.15, -10, 1.); // c2: linear coefficient
  RooRealVar c2Cand2CorrPrefit("c2Cand2CorrPrefit", "Corrected Prefit Linear coefficient", -3.15, -10, 1.); // c2: linear coefficient

  RooAbsPdf *bkgFuncPoly1Cand1Prefit = new RooPolynomial("bkgFuncPoly1Cand1Prefit", "background prefit function of candidate 1", *massCand1, RooArgList(c1Prefit, c2Prefit));
  RooAbsPdf *bkgFuncPoly1Cand2Prefit = new RooPolynomial("bkgFuncPoly1Cand2Prefit", "background prefit function of candidate 2", *massCand2, RooArgList(c1Cand2Prefit, c2Cand2Prefit));

  RooAbsPdf *bkgFuncPoly1Cand1 = new RooPolynomial("bkgFuncPoly1Cand1", "background fit function of candidate 1", *massCand1, RooArgList(c1, c2));
  RooAbsPdf *bkgFuncPoly1Cand2 = new RooPolynomial("bkgFuncPoly1Cand2", "background fit function of candidate 2", *massCand2, RooArgList(c1Cand2, c2Cand2));

  RooAbsPdf *bkgFuncPoly1Cand1Corr = new RooPolynomial("bkgFuncPoly1Cand1Corr", "Corrected background fit function of candidate 1", *massCand1, RooArgList(c1Corr, c2Corr));
  RooAbsPdf *bkgFuncPoly1Cand2Corr = new RooPolynomial("bkgFuncPoly1Cand2Corr", "Corrected background fit function of candidate 2", *massCand2, RooArgList(c1Cand2Corr, c2Cand2Corr));

  RooAbsPdf *bkgFuncPoly1Cand1CorrPrefit = new RooPolynomial("bkgFuncPoly1Cand1CorrPrefit", "Prefit Corrected background fit function of candidate 1", *massCand1, RooArgList(c1CorrPrefit, c2CorrPrefit));
  RooAbsPdf *bkgFuncPoly1Cand2CorrPrefit = new RooPolynomial("bkgFuncPoly1Cand2CorrPrefit", "Prefit Corrected background fit function of candidate 2", *massCand2, RooArgList(c1Cand2CorrPrefit, c2Cand2CorrPrefit));

  _workspace.import(*bkgFuncPoly1Cand1Prefit);
  _workspace.import(*bkgFuncPoly1Cand2Prefit);
  _workspace.import(*bkgFuncPoly1Cand1);
  _workspace.import(*bkgFuncPoly1Cand2);
  _workspace.import(*bkgFuncPoly1Cand1Corr);
  _workspace.import(*bkgFuncPoly1Cand2Corr);
  _workspace.import(*bkgFuncPoly1Cand1CorrPrefit);
  _workspace.import(*bkgFuncPoly1Cand2CorrPrefit);

  // bkg poly2: c1 + c2*x + c3*x^2
  RooRealVar c3Prefit("c3Prefit", "Prefit Quadratic coefficient", 0.3, 0., 1.);           // c3: quadratic coefficient
  RooRealVar c3Cand2Prefit("c3Cand2Prefit", "Prefit Quadratic coefficient", 0.3, 0., 1.); // c3: quadratic coefficient

  RooRealVar c3("c3", "Quadratic coefficient", 0.3, 0., 1.);           // c3: quadratic coefficient
  RooRealVar c3Cand2("c3Cand2", "Quadratic coefficient", 0.3, 0., 1.); // c3: quadratic coefficient

  RooRealVar c3Corr("c3Corr", "Corrected Quadratic coefficient", 0.3, -1., 1.);           // c3: quadratic coefficient
  RooRealVar c3Cand2Corr("c3Cand2Corr", "Corrected Quadratic coefficient", 0.3, -1., 1.); // c3: quadratic coefficient

  RooRealVar c3CorrPrefit("c3CorrPrefit", "Corrected Prefit Quadratic coefficient", 0.3, -1., 1.);           // c3: quadratic coefficient
  RooRealVar c3Cand2CorrPrefit("c3Cand2CorrPrefit", "Corrected Prefit Quadratic coefficient", 0.3, -1., 1.); // c3: quadratic coefficient

  RooAbsPdf *bkgFuncPoly2Cand1Prefit = new RooPolynomial("bkgFuncPoly2Cand1Prefit", "background prefit function of candidate 1", *massCand1, RooArgList(c1Prefit, c2Prefit, c3Prefit));
  RooAbsPdf *bkgFuncPoly2Cand2Prefit = new RooPolynomial("bkgFuncPoly2Cand2Prefit", "background prefit function of candidate 2", *massCand2, RooArgList(c1Cand2Prefit, c2Cand2Prefit, c3Cand2Prefit));

  RooAbsPdf *bkgFuncPoly2Cand1 = new RooPolynomial("bkgFuncPoly2Cand1", "background fit function of candidate 1", *massCand1, RooArgList(c1, c2, c3));
  RooAbsPdf *bkgFuncPoly2Cand2 = new RooPolynomial("bkgFuncPoly2Cand2", "background fit function of candidate 2", *massCand2, RooArgList(c1Cand2, c2Cand2, c3Cand2));

  RooAbsPdf *bkgFuncPoly2Cand1Corr = new RooPolynomial("bkgFuncPoly2Cand1Corr", "Corrected background fit function of candidate 1", *massCand1, RooArgList(c1Corr, c2Corr, c3Corr));
  RooAbsPdf *bkgFuncPoly2Cand2Corr = new RooPolynomial("bkgFuncPoly2Cand2Corr", "Correctedbackground fit function of candidate 2", *massCand2, RooArgList(c1Corr, c2Cand2Corr, c3Cand2Corr));

  RooAbsPdf *bkgFuncPoly2Cand1CorrPrefit = new RooPolynomial("bkgFuncPoly2Cand1CorrPrefit", "Corrected prefit background fit function of candidate 1", *massCand1, RooArgList(c1CorrPrefit, c2CorrPrefit, c3CorrPrefit));
  RooAbsPdf *bkgFuncPoly2Cand2CorrPrefit = new RooPolynomial("bkgFuncPoly2Cand2CorrPrefit", "Corrected prefit background fit function of candidate 2", *massCand2, RooArgList(c1CorrPrefit, c2Cand2CorrPrefit, c3Cand2CorrPrefit));

  _workspace.import(*bkgFuncPoly2Cand1Prefit);
  _workspace.import(*bkgFuncPoly2Cand2Prefit);
  _workspace.import(*bkgFuncPoly2Cand1);
  _workspace.import(*bkgFuncPoly2Cand2);
  _workspace.import(*bkgFuncPoly2Cand1Corr);
  _workspace.import(*bkgFuncPoly2Cand2Corr);
  _workspace.import(*bkgFuncPoly2Cand1CorrPrefit);
  _workspace.import(*bkgFuncPoly2Cand2CorrPrefit);

  // bkg poly2 * expo: (c1 + c2*x + c3*x^2) * exp(tau*x)
  RooAbsPdf *bkgFuncExpPoly2Cand1Prefit = new RooProdPdf("bkgFuncExpPoly2Cand1Prefit", "Prefit Exponential * Polynomial", RooArgList(*bkgFuncExpoCand1Prefit, *bkgFuncPoly2Cand1Prefit));
  RooAbsPdf *bkgFuncExpPoly2Cand2Prefit = new RooProdPdf("bkgFuncExpPoly2Cand2Prefit", "Prefit Exponential * Polynomial", RooArgList(*bkgFuncExpoCand2Prefit, *bkgFuncPoly2Cand2Prefit));

  RooAbsPdf *bkgFuncExpPoly2Cand1 = new RooProdPdf("bkgFuncExpPoly2Cand1", "Exponential * Polynomial", RooArgList(*bkgFuncExpoCand1, *bkgFuncPoly2Cand1));
  RooAbsPdf *bkgFuncExpPoly2Cand2 = new RooProdPdf("bkgFuncExpPoly2Cand2", "Exponential * Polynomial", RooArgList(*bkgFuncExpoCand2, *bkgFuncPoly2Cand2));

  RooAbsPdf *bkgFuncExpPoly2Cand1Corr = new RooProdPdf("bkgFuncExpPoly2Cand1Corr", "Corrected Exponential * Polynomial", RooArgList(*bkgFuncExpoCand1Corr, *bkgFuncPoly2Cand1Corr));
  RooAbsPdf *bkgFuncExpPoly2Cand2Corr = new RooProdPdf("bkgFuncExpPoly2Cand2Corr", "Corrected Exponential * Polynomial", RooArgList(*bkgFuncExpoCand2Corr, *bkgFuncPoly2Cand2Corr));

  RooAbsPdf *bkgFuncExpPoly2Cand1CorrPrefit = new RooProdPdf("bkgFuncExpPoly2Cand1CorrPrefit", "Prefit Corrected Exponential * Polynomial", RooArgList(*bkgFuncExpoCand1CorrPrefit, *bkgFuncPoly2Cand1CorrPrefit));
  RooAbsPdf *bkgFuncExpPoly2Cand2CorrPrefit = new RooProdPdf("bkgFuncExpPoly2Cand2CorrPrefit", "Prefit Corrected Exponential * Polynomial", RooArgList(*bkgFuncExpoCand2CorrPrefit, *bkgFuncPoly2Cand2CorrPrefit));

  _workspace.import(*bkgFuncExpPoly2Cand1Prefit, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly2Cand2Prefit, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly2Cand1, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly2Cand2, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly2Cand1Corr, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly2Cand2Corr, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly2Cand1CorrPrefit, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly2Cand2CorrPrefit, RooFit::RecycleConflictNodes());

  // bkg poly1 * expo: (c1 + c2*x) * exp(tau*x)
  RooAbsPdf *bkgFuncExpPoly1Cand1Prefit = new RooProdPdf("bkgFuncExpPoly1Cand1Prefit", "Prefit Exponential * Polynomial", RooArgList(*bkgFuncExpoCand1Prefit, *bkgFuncPoly1Cand1Prefit));
  RooAbsPdf *bkgFuncExpPoly1Cand2Prefit = new RooProdPdf("bkgFuncExpPoly1Cand2Prefit", "Prefit Exponential * Polynomial", RooArgList(*bkgFuncExpoCand2Prefit, *bkgFuncPoly1Cand2Prefit));

  RooAbsPdf *bkgFuncExpPoly1Cand1 = new RooProdPdf("bkgFuncExpPoly1Cand1", "Exponential * Polynomial", RooArgList(*bkgFuncExpoCand1, *bkgFuncPoly1Cand1));
  RooAbsPdf *bkgFuncExpPoly1Cand2 = new RooProdPdf("bkgFuncExpPoly1Cand2", "Exponential * Polynomial", RooArgList(*bkgFuncExpoCand2, *bkgFuncPoly1Cand2));

  RooAbsPdf *bkgFuncExpPoly1Cand1Corr = new RooProdPdf("bkgFuncExpPoly1Cand1Corr", "Corr Exponential * Polynomial", RooArgList(*bkgFuncExpoCand1Corr, *bkgFuncPoly1Cand1Corr));
  RooAbsPdf *bkgFuncExpPoly1Cand2Corr = new RooProdPdf("bkgFuncExpPoly1Cand2Corr", "Corr Exponential * Polynomial", RooArgList(*bkgFuncExpoCand2Corr, *bkgFuncPoly1Cand2Corr));

  RooAbsPdf *bkgFuncExpPoly1Cand1CorrPrefit = new RooProdPdf("bkgFuncExpPoly1Cand1CorrPrefit", "Prefit Corr Exponential * Polynomial", RooArgList(*bkgFuncExpoCand1CorrPrefit, *bkgFuncPoly1Cand1CorrPrefit));
  RooAbsPdf *bkgFuncExpPoly1Cand2CorrPrefit = new RooProdPdf("bkgFuncExpPoly1Cand2CorrPrefit", "Prefit Corr Exponential * Polynomial", RooArgList(*bkgFuncExpoCand2CorrPrefit, *bkgFuncPoly1Cand2CorrPrefit));

  _workspace.import(*bkgFuncExpPoly1Cand1Prefit, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly1Cand2Prefit, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly1Cand1, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly1Cand2, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly1Cand1Corr, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly1Cand2Corr, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly1Cand1CorrPrefit, RooFit::RecycleConflictNodes());
  _workspace.import(*bkgFuncExpPoly1Cand2CorrPrefit, RooFit::RecycleConflictNodes());


  // bkg a *exp(b*x + c*x^2)
  RooRealVar aExp2Prefit("aExp2Prefit", "aExp2Prefit", 5000, -1e3, 1e5); // aExp2: normalisation of the function
  RooRealVar bExp2Prefit("bExp2Prefit", "bExp2Prefit", 0.04, -0.5, 0.5); // bExp2: linear coefficient in the exponential
  RooRealVar cExp2Prefit("cExp2Prefit", "cExp2Prefit", -0.6, -3.5, 1.); // cExp2: quadratic coefficient in the exponential

  RooRealVar aExp2Cand2Prefit("aExp2Cand2Prefit", "aExp2Cand2Prefit", 5000, -1e3, 1e5); // aExp2: normalisation of the function
  RooRealVar bExp2Cand2Prefit("bExp2Cand2Prefit", "bExp2Cand2Prefit", 0.04, -0.5, 0.5); // bExp2: linear coefficient in the exponential
  RooRealVar cExp2Cand2Prefit("cExp2Cand2Prefit", "cExp2Cand2Prefit", -0.6, -3.5, 1.); // cExp2: quadratic coefficient in the exponential

  RooRealVar aExp2("aExp2", "aExp2", 5000, -1e3, 1e5); // aExp2: normalisation of the function
  RooRealVar bExp2("bExp2", "bExp2", 0.04, -0.5, 0.5); // bExp2: linear coefficient in the exponential
  RooRealVar cExp2("cExp2", "cExp2", -0.6, -3.5, 1.); // cExp2: quadratic coefficient in the exponential

  RooRealVar aExp2Cand2("aExp2Cand2", "aExp2Cand2", 5000, -1e3, 1e5); // aExp2: normalisation of the function
  RooRealVar bExp2Cand2("bExp2Cand2", "bExp2Cand2", 0.04, -0.5, 0.5); // bExp2: linear coefficient in the exponential
  RooRealVar cExp2Cand2("cExp2Cand2", "cExp2Cand2", -0.6, -3.5, 1.); // cExp2: quadratic coefficient in the exponential

  RooRealVar aExp2Corr("aExp2Corr", "aExp2Corr", 5000, -1e3, 1e5); // aExp2: normalisation of the function
  RooRealVar bExp2Corr("bExp2Corr", "bExp2Corr", 0.04, -0.5, 0.5); // bExp2: linear coefficient in the exponential
  RooRealVar cExp2Corr("cExp2Corr", "cExp2Corr", -0.6, -3.5, 1.); // cExp2: quadratic coefficient in the exponential

  RooRealVar aExp2Cand2Corr("aExp2Cand2Corr", "aExp2Cand2Corr", 5000, -1e3, 1e5); // aExp2: normalisation of the function
  RooRealVar bExp2Cand2Corr("bExp2Cand2Corr", "bExp2Cand2Corr", 0.04, -0.5, 0.5); // bExp2: linear coefficient in the exponential
  RooRealVar cExp2Cand2Corr("cExp2Cand2Corr", "cExp2Cand2Corr", -0.6, -3.5, 1.); // cExp2: quadratic coefficient in the exponential

  RooRealVar aExp2CorrPrefit("aExp2CorrPrefit", "aExp2CorrPrefit", 5000, -1e3, 1e5); // aExp2: normalisation of the function
  RooRealVar bExp2CorrPrefit("bExp2CorrPrefit", "bExp2CorrPrefit", 0.04, -0.5, 0.5); // bExp2: linear coefficient in the exponential
  RooRealVar cExp2CorrPrefit("cExp2CorrPrefit", "cExp2CorrPrefit", -0.6, -3.5, 1.); // cExp2: quadratic coefficient in the exponential

  RooRealVar aExp2Cand2CorrPrefit("aExp2Cand2CorrPrefit", "aExp2Cand2CorrPrefit", 5000, -1e3, 1e5); // aExp2: normalisation of the function
  RooRealVar bExp2Cand2CorrPrefit("bExp2Cand2CorrPrefit", "bExp2Cand2CorrPrefit", 0.04, -0.5, 0.5); // bExp2: linear coefficient in the exponential
  RooRealVar cExp2Cand2CorrPrefit("cExp2Cand2CorrPrefit", "cExp2Cand2CorrPrefit", -0.6, -3.5, 1.); // cExp2: quadratic coefficient in the exponential

  RooAbsPdf *bkgFuncExp2Cand1Prefit = new RooGenericPdf("bkgFuncExp2Cand1Prefit", "bkgFuncExp2Cand1Prefit", "@1*exp(@2*@0+@3*@0*@0)", RooArgList(*massCand1, aExp2Prefit, bExp2Prefit, cExp2Prefit));
  RooAbsPdf *bkgFuncExp2Cand2Prefit = new RooGenericPdf("bkgFuncExp2Cand2Prefit", "bkgFuncExp2Cand2Prefit", "@1*exp(@2*@0+@3*@0*@0)", RooArgList(*massCand2, aExp2Prefit, bExp2Cand2Prefit, cExp2Cand2Prefit));

  RooAbsPdf *bkgFuncExp2Cand1 = new RooGenericPdf("bkgFuncExp2Cand1", "bkgFuncExp2Cand1", "@1*exp(@2*@0+@3*@0*@0)", RooArgList(*massCand1, aExp2, bExp2, cExp2));
  RooAbsPdf *bkgFuncExp2Cand2 = new RooGenericPdf("bkgFuncExp2Cand2", "bkgFuncExp2Cand2", "@1*exp(@2*@0+@3*@0*@0)", RooArgList(*massCand2, aExp2, bExp2Cand2, cExp2Cand2));

  RooAbsPdf *bkgFuncExp2Cand1Corr = new RooGenericPdf("bkgFuncExp2Cand1Corr", "bkgFuncExp2Cand1Corr", "@1*exp(@2*@0+@3*@0*@0)", RooArgList(*massCand1, aExp2Corr, bExp2Corr, cExp2Corr));
  RooAbsPdf *bkgFuncExp2Cand2Corr = new RooGenericPdf("bkgFuncExp2Cand2Corr", "bkgFuncExp2Cand2Corr", "@1*exp(@2*@0+@3*@0*@0)", RooArgList(*massCand2, aExp2Corr, bExp2Cand2Corr, cExp2Cand2Corr));

  RooAbsPdf *bkgFuncExp2Cand1CorrPrefit = new RooGenericPdf("bkgFuncExp2Cand1CorrPrefit", "bkgFuncExp2Cand1CorrPrefit", "@1*exp(@2*@0+@3*@0*@0)", RooArgList(*massCand1, aExp2CorrPrefit, bExp2CorrPrefit, cExp2CorrPrefit));
  RooAbsPdf *bkgFuncExp2Cand2CorrPrefit = new RooGenericPdf("bkgFuncExp2Cand2CorrPrefit", "bkgFuncExp2Cand2CorrPrefit", "@1*exp(@2*@0+@3*@0*@0)", RooArgList(*massCand2, aExp2CorrPrefit, bExp2Cand2CorrPrefit, cExp2Cand2CorrPrefit));

  _workspace.import(*bkgFuncExp2Cand1Prefit);
  _workspace.import(*bkgFuncExp2Cand2Prefit);
  _workspace.import(*bkgFuncExp2Cand1);
  _workspace.import(*bkgFuncExp2Cand2);
  _workspace.import(*bkgFuncExp2Cand1Corr);
  _workspace.import(*bkgFuncExp2Cand2Corr);
  _workspace.import(*bkgFuncExp2Cand1CorrPrefit);
  _workspace.import(*bkgFuncExp2Cand2CorrPrefit);

  /// | ------------------------------------------------------------------ |
  /// | ----------------------- SIGNAL FUNCTIONS ------------------------- |
  /// | ------------------------------------------------------------------ |
  // signal pdf
  RooRealVar meanPrefit("meanPrefit", "Prefit mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigmaPrefit("sigmaPrefit", "Prefit sigma for signal", 0.02, 0.01, 0.04);

  RooRealVar mean("mean", "mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigma("sigma", "sigma for signal", 0.02, 0.01, 0.04);

  RooRealVar meanCorr("meanCorr", "Corr mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigmaCorr("sigmaCorr", "Corr sigma for signal", 0.02, 0.01, 0.04);

  RooRealVar meanCorrPrefit("meanCorrPrefit", "Corr mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigmaCorrPrefit("sigmaCorrPrefit", "Corr sigma for signal", 0.02, 0.01, 0.04);

  RooRealVar meanCand2Prefit("meanCand2Prefit", "Prefit mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigmaCand2Prefit("sigmaCand2Prefit", "Prefit sigma for signal", 0.02, 0.01, 0.04);

  RooRealVar meanCand2("meanCand2", "mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigmaCand2("sigmaCand2", "sigma for signal", 0.02, 0.01, 0.04);

  RooRealVar meanCand2Corr("meanCand2Corr", "Corr mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigmaCand2Corr("sigmaCand2Corr", "Corr sigma for signal", 0.02, 0.01, 0.04);

  RooRealVar meanCand2CorrPrefit("meanCand2CorrPrefit", "Corr mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigmaCand2CorrPrefit("sigmaCand2CorrPrefit", "Corr sigma for signal", 0.02, 0.01, 0.04);

  RooAbsPdf *sgnFuncGausCand1Prefit = new RooGaussian("sgnFuncGausCand1Prefit", "Prefit signal pdf of candidate 1", *massCand1, meanPrefit, sigmaPrefit);
  RooAbsPdf *sgnFuncGausCand2Prefit = new RooGaussian("sgnFuncGausCand2Prefit", "Prefit signal pdf of candidate 2", *massCand2, meanCand2Prefit, sigmaCand2Prefit);

  RooAbsPdf *sgnFuncGausCand1 = new RooGaussian("sgnFuncGausCand1", "signal pdf of candidate 1", *massCand1, mean, sigma);
  RooAbsPdf *sgnFuncGausCand2 = new RooGaussian("sgnFuncGausCand2", "signal pdf of candidate 2", *massCand2, meanCand2, sigmaCand2);

  RooAbsPdf *sgnFuncGausCand1Corr = new RooGaussian("sgnFuncGausCand1Corr", "Corr signal pdf of candidate 1", *massCand1, meanCorr, sigmaCorr);
  RooAbsPdf *sgnFuncGausCand2Corr = new RooGaussian("sgnFuncGausCand2Corr", "Corr signal pdf of candidate 2", *massCand2, meanCand2Corr, sigmaCand2Corr);

  RooAbsPdf *sgnFuncGausCand1CorrPrefit = new RooGaussian("sgnFuncGausCand1CorrPrefit", "Corr signal pdf of candidate 1", *massCand1, meanCorrPrefit, sigmaCorrPrefit);
  RooAbsPdf *sgnFuncGausCand2CorrPrefit = new RooGaussian("sgnFuncGausCand2CorrPrefit", "Corr signal pdf of candidate 2", *massCand2, meanCand2CorrPrefit, sigmaCand2CorrPrefit);

  _workspace.import(*sgnFuncGausCand1Prefit);
  _workspace.import(*sgnFuncGausCand2Prefit);
  _workspace.import(*sgnFuncGausCand1);
  _workspace.import(*sgnFuncGausCand2);
  _workspace.import(*sgnFuncGausCand1Corr);
  _workspace.import(*sgnFuncGausCand2Corr);
  _workspace.import(*sgnFuncGausCand1CorrPrefit);
  _workspace.import(*sgnFuncGausCand2CorrPrefit);

  RooRealVar alpha1("alpha1", "Alpha (transition)", 1.595, 1.4, 1.8);
  RooRealVar n1("n1", "n (steepness of tail)", 2.999, 2.8, 3.2);
  RooRealVar alpha2("alpha2", "Alpha (transition)", 1.4262, 1.2, 1.6);
  RooRealVar n2("n2", "n (steepness of tail)", 3.8174, 3.6, 4.0);

  alpha1.setConstant(kTRUE);
  n1.setConstant(kTRUE);
  alpha2.setConstant(kTRUE);
  n2.setConstant(kTRUE);

  RooAbsPdf *sgnFuncCBCand1Prefit = new BifurcatedCB("sgnFuncCBCand1Prefit", "Prefit signal pdf of candidate 1", *massCand1, meanPrefit, sigmaPrefit, alpha1, n1, alpha2, n2);
  RooAbsPdf *sgnFuncCBCand2Prefit = new BifurcatedCB("sgnFuncCBCand2Prefit", "Prefit signal pdf of candidate 2", *massCand2, meanPrefit, sigmaPrefit, alpha1, n1, alpha2, n2);

  RooAbsPdf *sgnFuncCBCand1 = new BifurcatedCB("sgnFuncCBCand1", "signal pdf of candidate 1", *massCand1, mean, sigma, alpha1, n1, alpha2, n2);
  RooAbsPdf *sgnFuncCBCand2 = new BifurcatedCB("sgnFuncCBCand2", "signal pdf of candidate 2", *massCand2, mean, sigma, alpha1, n1, alpha2, n2);

  RooAbsPdf *sgnFuncCBCand1Corr = new BifurcatedCB("sgnFuncCBCand1Corr", "Corr signal pdf of candidate 1", *massCand1, meanCorr, sigmaCorr, alpha1, n1, alpha2, n2);
  RooAbsPdf *sgnFuncCBCand2Corr = new BifurcatedCB("sgnFuncCBCand2Corr", "Corr signal pdf of candidate 2", *massCand2, meanCorr, sigmaCorr, alpha1, n1, alpha2, n2);

  RooAbsPdf *sgnFuncCBCand1CorrPrefit = new BifurcatedCB("sgnFuncCBCand1CorrPrefit", "Corr signal pdf of candidate 1", *massCand1, meanCorrPrefit, sigmaCorrPrefit, alpha1, n1, alpha2, n2);
  RooAbsPdf *sgnFuncCBCand2CorrPrefit = new BifurcatedCB("sgnFuncCBCand2CorrPrefit", "Corr signal pdf of candidate 2", *massCand2, meanCorrPrefit, sigmaCorrPrefit, alpha1, n1, alpha2, n2);

  _workspace.import(*sgnFuncCBCand1Prefit);
  _workspace.import(*sgnFuncCBCand2Prefit);
  _workspace.import(*sgnFuncCBCand1);
  _workspace.import(*sgnFuncCBCand2);
  _workspace.import(*sgnFuncCBCand1Corr);
  _workspace.import(*sgnFuncCBCand2Corr);
  _workspace.import(*sgnFuncCBCand1CorrPrefit);
  _workspace.import(*sgnFuncCBCand2CorrPrefit);

  /// | ------------------------------------------------------------------ |
  /// | ------------------- REFLECTION FUNCTIONS ------------------------- |
  /// | ------------------------------------------------------------------ |
  // reflection Gaussian
  RooRealVar meanReflPrefit("meanReflPrefit", "Prefit mean for reflections", 1.85, 0.0, 2.15);
  RooRealVar sigmaReflPrefit("sigmaReflPrefit", "Prefit sigma for reflection", 0.012, 0.0001, 0.3);

  RooRealVar meanRefl("meanRefl", "mean for reflections", 1.85, 0.0, 2.15);
  RooRealVar sigmaRefl("sigmaRefl", "sigma for reflection", 0.012, 0.0001, 0.3);

  RooRealVar meanReflCorr("meanReflCorr", "Corr mean for reflections", 1.85, 0.0, 2.15);
  RooRealVar sigmaReflCorr("sigmaReflCorr", "Corr sigma for reflection", 0.012, 0.0001, 0.3);

  RooRealVar meanReflCorrPrefit("meanReflCorrPrefit", "Prefit Corr mean for reflections", 1.85, 0.0, 2.15);
  RooRealVar sigmaReflCorrPrefit("sigmaReflCorrPrefit", "Prefit Corr sigma for reflection", 0.012, 0.0001, 0.3);

  RooAbsPdf *reflFuncGausCand1Prefit = new RooGaussian("reflFuncGausCand1Prefit", "Prefit reflection pdf of candidate 1", *massCand1, meanReflPrefit, sigmaReflPrefit);
  RooAbsPdf *reflFuncGausCand2Prefit = new RooGaussian("reflFuncGausCand2Prefit", "Prefit reflection pdf of candidate 2", *massCand2, meanReflPrefit, sigmaReflPrefit);

  RooAbsPdf *reflFuncGausCand1 = new RooGaussian("reflFuncGausCand1", "reflection pdf of candidate 1", *massCand1, meanRefl, sigmaRefl);
  RooAbsPdf *reflFuncGausCand2 = new RooGaussian("reflFuncGausCand2", "reflection pdf of candidate 2", *massCand2, meanRefl, sigmaRefl);

  RooAbsPdf *reflFuncGausCand1Corr = new RooGaussian("reflFuncGausCand1Corr", "Corr reflection pdf of candidate 1", *massCand1, meanReflCorr, sigmaReflCorr);
  RooAbsPdf *reflFuncGausCand2Corr = new RooGaussian("reflFuncGausCand2Corr", "Corr reflection pdf of candidate 2", *massCand2, meanReflCorr, sigmaReflCorr);

  RooAbsPdf *reflFuncGausCand1CorrPrefit = new RooGaussian("reflFuncGausCand1CorrPrefit", "Corr reflection pdf of candidate 1", *massCand1, meanReflCorrPrefit, sigmaReflCorrPrefit);
  RooAbsPdf *reflFuncGausCand2CorrPrefit = new RooGaussian("reflFuncGausCand2CorrPrefit", "Corr reflection pdf of candidate 2", *massCand2, meanReflCorrPrefit, sigmaReflCorrPrefit);

  _workspace.import(*reflFuncGausCand1Prefit);
  _workspace.import(*reflFuncGausCand2Prefit);
  _workspace.import(*reflFuncGausCand1);
  _workspace.import(*reflFuncGausCand2);
  _workspace.import(*reflFuncGausCand1Corr);
  _workspace.import(*reflFuncGausCand2Corr);
  _workspace.import(*reflFuncGausCand1CorrPrefit);
  _workspace.import(*reflFuncGausCand2CorrPrefit);

  // reflection double gaussian
  RooRealVar meanReflDoubleGausPrefit("meanReflDoubleGausPrefit", "Prefit mean for reflection double gaussian", 1.85, 0.0, 1.90);
  RooRealVar sigmaReflDoubleGausPrefit("sigmaReflDoubleGausPrefit", "Prefit sigmaReflDoubleGaus", 0.012, 0.0001, 0.2);
  RooRealVar fracReflPrefit("fracReflPrefit", "Prefit frac of the two reflected gaussians", 0.5, 0, 1.);

  RooRealVar meanReflDoubleGaus("meanReflDoubleGaus", "mean for reflection double gaussian", 1.85, 0.0, 1.90);
  RooRealVar sigmaReflDoubleGaus("sigmaReflDoubleGaus", "sigmaReflDoubleGaus", 0.012, 0.0001, 0.2);
  RooRealVar fracRefl("fracRefl", "frac of the two reflected gaussians", 0.5, 0, 1.);

  RooRealVar meanReflDoubleGausCorr("meanReflDoubleGausCorr", "Corr mean for reflection double gaussian", 1.85, 0.0, 1.90);
  RooRealVar sigmaReflDoubleGausCorr("sigmaReflDoubleGausCorr", "Corr sigmaReflDoubleGaus", 0.012, 0.0001, 0.2);
  RooRealVar fracReflCorr("fracReflCorr", "Corr frac of the two reflected gaussians", 0.5, 0, 1.);

  RooRealVar meanReflDoubleGausCorrPrefit("meanReflDoubleGausCorrPrefit", "Prefit Corr mean for reflection double gaussian", 1.85, 0.0, 1.90);
  RooRealVar sigmaReflDoubleGausCorrPrefit("sigmaReflDoubleGausCorrPrefit", "Prefit Corr sigmaReflDoubleGaus", 0.012, 0.0001, 0.2);
  RooRealVar fracReflCorrPrefit("fracReflCorrPrefit", "Prefit Corr frac of the two reflected gaussians", 0.5, 0, 1.);

  //    Second gaussian definition
  RooAbsPdf *relfFuncSecondGausCand1Prefit = new RooGaussian("relfFuncSecondGausCand1Prefit", "relfFuncSecondGausCand1Prefit", *massCand1, meanReflDoubleGausPrefit, sigmaReflDoubleGausPrefit);
  RooAbsPdf *relfFuncSecondGausCand2Prefit = new RooGaussian("relfFuncSecondGausCand2Prefit", "relfFuncSecondGausCand2Prefit", *massCand2, meanReflDoubleGausPrefit, sigmaReflDoubleGausPrefit);

  RooAbsPdf *relfFuncSecondGausCand1 = new RooGaussian("relfFuncSecondGausCand1", "relfFuncSecondGausCand1", *massCand1, meanReflDoubleGaus, sigmaReflDoubleGaus);
  RooAbsPdf *relfFuncSecondGausCand2 = new RooGaussian("relfFuncSecondGausCand2", "relfFuncSecondGausCand2", *massCand2, meanReflDoubleGaus, sigmaReflDoubleGaus);

  RooAbsPdf *relfFuncSecondGausCand1Corr = new RooGaussian("relfFuncSecondGausCand1Corr", "relfFuncSecondGausCand1Corr", *massCand1, meanReflDoubleGausCorr, sigmaReflDoubleGausCorr);
  RooAbsPdf *relfFuncSecondGausCand2Corr = new RooGaussian("relfFuncSecondGausCand2Corr", "relfFuncSecondGausCand2Corr", *massCand2, meanReflDoubleGausCorr, sigmaReflDoubleGausCorr);

  RooAbsPdf *relfFuncSecondGausCand1CorrPrefit = new RooGaussian("relfFuncSecondGausCand1CorrPrefit", "relfFuncSecondGausCand1CorrPrefit", *massCand1, meanReflDoubleGausCorrPrefit, sigmaReflDoubleGausCorrPrefit);
  RooAbsPdf *relfFuncSecondGausCand2CorrPrefit = new RooGaussian("relfFuncSecondGausCand2CorrPrefit", "relfFuncSecondGausCand2CorrPrefit", *massCand2, meanReflDoubleGausCorrPrefit, sigmaReflDoubleGausCorrPrefit);

  //    Compose the final double gaussian
  RooAbsPdf *reflFuncDoubleGausCand1Prefit = new RooAddPdf("reflFuncDoubleGausCand1Prefit", "Prefit reflection pdf of candidate 1", RooArgList(*reflFuncGausCand1Prefit, *relfFuncSecondGausCand1Prefit), fracReflPrefit);
  RooAbsPdf *reflFuncDoubleGausCand2Prefit = new RooAddPdf("reflFuncDoubleGausCand2Prefit", "Prefit reflection pdf of candidate 2", RooArgList(*reflFuncGausCand2Prefit, *relfFuncSecondGausCand2Prefit), fracReflPrefit);

  RooAbsPdf *reflFuncDoubleGausCand1 = new RooAddPdf("reflFuncDoubleGausCand1", "reflection pdf of candidate 1", RooArgList(*reflFuncGausCand1, *relfFuncSecondGausCand1), fracRefl);
  RooAbsPdf *reflFuncDoubleGausCand2 = new RooAddPdf("reflFuncDoubleGausCand2", "reflection pdf of candidate 2", RooArgList(*reflFuncGausCand2, *relfFuncSecondGausCand2), fracRefl);

  RooAbsPdf *reflFuncDoubleGausCand1Corr = new RooAddPdf("reflFuncDoubleGausCand1Corr", "Corr reflection pdf of candidate 1", RooArgList(*reflFuncGausCand1Corr, *relfFuncSecondGausCand1Corr), fracReflCorr);
  RooAbsPdf *reflFuncDoubleGausCand2Corr = new RooAddPdf("reflFuncDoubleGausCand2Corr", "Corr reflection pdf of candidate 2", RooArgList(*reflFuncGausCand2Corr, *relfFuncSecondGausCand2Corr), fracReflCorr);

  RooAbsPdf *reflFuncDoubleGausCand1CorrPrefit = new RooAddPdf("reflFuncDoubleGausCand1CorrPrefit", "Prefit Corr reflection pdf of candidate 1", RooArgList(*reflFuncGausCand1CorrPrefit, *relfFuncSecondGausCand1CorrPrefit), fracReflCorrPrefit);
  RooAbsPdf *reflFuncDoubleGausCand2CorrPrefit = new RooAddPdf("reflFuncDoubleGausCand2CorrPrefit", "Prefit Corr reflection pdf of candidate 2", RooArgList(*reflFuncGausCand2CorrPrefit, *relfFuncSecondGausCand2CorrPrefit), fracReflCorrPrefit);

  _workspace.import(*reflFuncDoubleGausCand1Prefit, RooFit::RecycleConflictNodes());
  _workspace.import(*reflFuncDoubleGausCand2Prefit, RooFit::RecycleConflictNodes());
  _workspace.import(*reflFuncDoubleGausCand1, RooFit::RecycleConflictNodes());
  _workspace.import(*reflFuncDoubleGausCand2, RooFit::RecycleConflictNodes());
  _workspace.import(*reflFuncDoubleGausCand1Corr, RooFit::RecycleConflictNodes());
  _workspace.import(*reflFuncDoubleGausCand2Corr, RooFit::RecycleConflictNodes());
  _workspace.import(*reflFuncDoubleGausCand1CorrPrefit, RooFit::RecycleConflictNodes());
  _workspace.import(*reflFuncDoubleGausCand2CorrPrefit, RooFit::RecycleConflictNodes());

  cout << "Workspace filled with functions" << endl;
}

void InvMassFitter2D::do2DFit(Bool_t draw, Bool_t doReflections, Bool_t isMc, TFile *fout)
{
  createDataset();
  // Declare observable variable
  _workspace.Print("v");
  RooDataSet *dataset = dynamic_cast<RooDataSet *>(_workspace.data("data"));
  RooDataSet *weightedDataset = dynamic_cast<RooDataSet *>(_workspace.data("weightedData"));

  if (!dataset) {
    cerr << "ERROR: dataset not found!" << endl;
    return;
  }
  double numEntries = dataset->numEntries();

  RooRealVar *massCand1 = _workspace.var("fMCand1");
  RooRealVar *massCand2 = _workspace.var("fMCand2");
  RooRealVar *ptCand1 = _workspace.var("fPtCand1");
  RooRealVar *ptCand2 = _workspace.var("fPtCand2");
  RooRealVar *yCand1 = _workspace.var("fYCand1");
  RooRealVar *yCand2 = _workspace.var("fYCand2");
  RooRealVar *phiCand1 = _workspace.var("fPhiCand1");
  RooRealVar *phiCand2 = _workspace.var("fPhiCand2");

  if (!massCand1 || !massCand2)
  {
    cerr << "Error: Variables not found in workspace!" << endl;
    return;
  }
  fillWorkspace(dataset);

  massCand1->setRange(_massMin, _massMax);
  massCand2->setRange(_massMin, _massMax);

  cout << "\n\n                  Integrated efficiency: " << _integratedEfficiency << endl;

  /// ----------------------------------------------------------------
  /// ------- LOAD PDF FUNCTIONS (SEVERAL OPTIONS POSSIBLE) ----------
  /// ----------------------------------------------------------------

  RooAbsPdf *bkgPdfCand1Prefit;
  RooAbsPdf *bkgPdfCand2Prefit;
  RooAbsPdf *sgnPdfCand1Prefit;
  RooAbsPdf *sgnPdfCand2Prefit;

  RooAbsPdf *bkgPdfCand1PrefitCorr;
  RooAbsPdf *bkgPdfCand2PrefitCorr;
  RooAbsPdf *sgnPdfCand1PrefitCorr;
  RooAbsPdf *sgnPdfCand2PrefitCorr;

  RooAbsPdf *reflPdfCand1Prefit;
  RooAbsPdf *reflPdfCand2Prefit;
  RooAbsPdf *reflPdfCand1PrefitCorr;
  RooAbsPdf *reflPdfCand2PrefitCorr;

  selectFitFunctions(sgnPdfCand1Prefit, sgnPdfCand2Prefit, bkgPdfCand1Prefit, bkgPdfCand2Prefit, reflPdfCand1Prefit, reflPdfCand2Prefit, Prefit);

  selectFitFunctions(sgnPdfCand1PrefitCorr, sgnPdfCand2PrefitCorr, bkgPdfCand1PrefitCorr, bkgPdfCand2PrefitCorr, reflPdfCand1PrefitCorr, reflPdfCand2PrefitCorr, PrefitCorr);

  RooFitResult *prefitResultCand1 = nullptr;
  RooFitResult *prefitResultCand2 = nullptr;
  RooFitResult *prefitCorrResultCand1 = nullptr;
  RooFitResult *prefitCorrResultCand2 = nullptr;
  if (doReflections) {
    prefitResultCand1 = fitAndPlot1DCandidate(sgnPdfCand1Prefit, bkgPdfCand1Prefit, reflPdfCand1Prefit, *massCand1, dataset, "Cand1", Form("fit_result_cand1_%s", "prefit"));
    prefitResultCand2 = fitAndPlot1DCandidate(sgnPdfCand2Prefit, bkgPdfCand2Prefit, reflPdfCand2Prefit, *massCand2, dataset, "Cand2", Form("fit_result_cand2_%s", "prefit"));

    prefitCorrResultCand1 = fitAndPlot1DCandidate(sgnPdfCand1PrefitCorr, bkgPdfCand1PrefitCorr, reflPdfCand1PrefitCorr, *massCand1, weightedDataset, "Cand1", Form("fit_result_weighted_cand1_%s", "prefitCorr"));
    prefitCorrResultCand2 = fitAndPlot1DCandidate(sgnPdfCand2PrefitCorr, bkgPdfCand2PrefitCorr, reflPdfCand2PrefitCorr, *massCand2, weightedDataset, "Cand2", Form("fit_result_weighted_cand2_%s", "prefitCorr"));
  } else {
    prefitResultCand1 = fitAndPlot1DCandidate(sgnPdfCand1Prefit, bkgPdfCand1Prefit, *massCand1, dataset, "Cand1", Form("fit_result_cand1_%s", "prefit"));
    prefitResultCand2 = fitAndPlot1DCandidate(sgnPdfCand2Prefit, bkgPdfCand2Prefit, *massCand2, dataset, "Cand2", Form("fit_result_cand2_%s", "prefit"));

    prefitCorrResultCand1 = fitAndPlot1DCandidate(sgnPdfCand1PrefitCorr, bkgPdfCand1PrefitCorr, *massCand1, weightedDataset, "Cand1", Form("fit_result_weighted_cand1_%s", "prefitCorr"));
    prefitCorrResultCand2 = fitAndPlot1DCandidate(sgnPdfCand2PrefitCorr, bkgPdfCand2PrefitCorr, *massCand2, weightedDataset, "Cand2", Form("fit_result_weighted_cand2_%s", "prefitCorr"));
  }

  /// ----------------------------------------------------------
  /// ------------ LOAD AND FIX PARAMETERS ---------------------
  /// ----------------------------------------------------------
  //  Most parameters are just set as initial values
  // Assuming fitResult is a valid RooFitResult* from a fit
  RooAbsPdf *bkgPdfCand1;
  RooAbsPdf *bkgPdfCand2;
  RooAbsPdf *sgnPdfCand1;
  RooAbsPdf *sgnPdfCand2;
  RooAbsPdf *reflPdfCand1;
  RooAbsPdf *reflPdfCand2;

  selectFitFunctions(sgnPdfCand1, sgnPdfCand2, bkgPdfCand1, bkgPdfCand2, reflPdfCand1, reflPdfCand2, Rawfit);

  RooAbsPdf *bkgPdfCand1Corr;
  RooAbsPdf *bkgPdfCand2Corr;
  RooAbsPdf *sgnPdfCand1Corr;
  RooAbsPdf *sgnPdfCand2Corr;
  RooAbsPdf *reflPdfCand1Corr;
  RooAbsPdf *reflPdfCand2Corr;

  selectFitFunctions(sgnPdfCand1Corr, sgnPdfCand2Corr, bkgPdfCand1Corr, bkgPdfCand2Corr, reflPdfCand1Corr, reflPdfCand2Corr, Corrfit);

  const RooArgList& paramsPrefitCand1 = prefitResultCand1->floatParsFinal();
  const RooArgList& paramsPrefitCand2 = prefitResultCand2->floatParsFinal();
  const RooArgList& paramsPrefitCorrCand1 = prefitCorrResultCand1->floatParsFinal();
  const RooArgList& paramsPrefitCorrCand2 = prefitCorrResultCand2->floatParsFinal();

  setPrefitParameters(paramsPrefitCand1);
  setPrefitParameters(paramsPrefitCand2);

  setPrefitParameters(paramsPrefitCorrCand1);
  setPrefitParameters(paramsPrefitCorrCand2);

  RooRealVar *mean = _workspace.var("mean");
  RooRealVar *sigma = _workspace.var("sigma");

  cout << "Mean, sigma and tau values set" << endl;

  RooRealVar *sigmaRefl = _workspace.var("sigmaRefl");
  RooRealVar *meanRefl = _workspace.var("meanRefl");
  RooRealVar *sigmaReflDoubleGaus = _workspace.var("sigmaReflDoubleGaus");
  RooRealVar *meanReflDoubleGaus = _workspace.var("meanReflDoubleGaus");
  RooRealVar *fracRefl = _workspace.var("fracRefl");

  if (!reflPdfCand1 || !reflPdfCand2)
  {
    cerr << "ERROR: attempted to access an empty pointer (reflPdf functions)!" << endl;
    return;
  }

  /// -----------------------------------------------------------------
  /// --------------- DEFINE YIELD PARAMETERS FOR FIT -----------------
  /// -----------------------------------------------------------------

  RooArgSet normSet(*massCand1, *massCand2);
  RooArgSet weightedNormSet(*massCand1, *massCand2);

  // RAW PARAMETERS
  RooRealVar *nBkg1 = new RooRealVar("nBkg1", "background yield of cand 1", 300.0, 50.0, 1000.0);
  RooRealVar *nSgn1 = new RooRealVar("nSgn1", "signal yield of cand1", 70.0, 0.0, 300.0);
  RooFormulaVar *nRefl1 = new RooFormulaVar("nRefl1", "reflected signal yield of cand1", "@0 * @1", RooArgList(_reflOverSgn, *nSgn1));

  RooRealVar *nBkg2 = new RooRealVar("nBkg2", "background yield of cand2", 300.0, 50.0, 1000.0);
  RooRealVar *nSgn2 = new RooRealVar("nSgn2", "signal yield of cand2", 70.0, 0.0, 300.0);
  RooFormulaVar *nRefl2 = new RooFormulaVar("nRefl2", "reflected signal yield of cand 2", "@0 * @1", RooArgList(_reflOverSgn, *nSgn2));

  RooFormulaVar nBkgBkg("nBkgBkg", "fraction of bkg x bkg", "@0 * @1", RooArgList(*nBkg1, *nBkg2));
  RooFormulaVar nBkgSgn("nBkgSgn", "fraction of bkg x sgn", "@0 * @1", RooArgList(*nBkg1, *nSgn2));
  RooFormulaVar nSgnBkg("nSgnBkg", "fraction of sgn x bkg", "@0 * @1", RooArgList(*nSgn1, *nBkg2));
  RooFormulaVar nBkgRefl("nBkgRefl", "fraction of bkg x refl", "@0 * @1", RooArgList(*nBkg1, *nRefl2));
  RooFormulaVar nReflBkg("nReflBkg", "fraction of refl x bkg", "@0 * @1", RooArgList(*nRefl1, *nBkg2));
  RooFormulaVar nSgnSgn("nSgnSgn", "fraction of sgn x sgn", "@0 * @1", RooArgList(*nSgn1, *nSgn2));
  RooFormulaVar nReflRefl("nReflRefl", "fraction of refl x refl", "@0 * @1", RooArgList(*nRefl1, *nRefl2));
  RooFormulaVar nSgnRefl("nSgnRefl", "fraction of sgn x refl", "@0 * @1", RooArgList(*nSgn1, *nRefl2));
  RooFormulaVar nReflSgn("nReflSgn", "fraction of refl x sgn", "@0 * @1", RooArgList(*nRefl1, *nSgn2));

  // WEIGHTED PARAMETERS
  RooRealVar *nBkgCorr1 = new RooRealVar("nBkgCorr1", "Corrected background yield of cand 1", 10000, 100, 200000);
  RooRealVar *nSgnCorr1 = new RooRealVar("nSgnCorr1", "Corrected signal yield of cand1", 2500, 10.0, 50000.0);
  RooFormulaVar *nReflCorr1 = new RooFormulaVar("nReflCorr1", "Corrected reflected signal yield of cand1", "@0 * @1", RooArgList(_reflOverSgn, *nSgnCorr1));

  RooRealVar *nBkgCorr2 = new RooRealVar("nBkgCorr2", "Corrected background yield of cand2", 10000, 100, 200000);
  RooRealVar *nSgnCorr2 = new RooRealVar("nSgnCorr2", "Corrected signal yield of cand2", 2500, 10.0, 50000.0);
  RooFormulaVar *nReflCorr2 = new RooFormulaVar("nReflCorr2", "Corrected reflected signal yield of cand 2", "@0 * @1", RooArgList(_reflOverSgn, *nSgnCorr2));

  RooFormulaVar nBkgBkgCorr("nBkgBkgCorr", "Corrected fraction of bkg x bkg", "@0 * @1", RooArgList(*nBkgCorr1, *nBkgCorr2));
  RooFormulaVar nBkgSgnCorr("nBkgSgnCorr", "Corrected fraction of bkg x sgn", "@0 * @1", RooArgList(*nBkgCorr1, *nSgnCorr2));
  RooFormulaVar nSgnBkgCorr("nSgnBkgCorr", "Corrected fraction of sgn x bkg", "@0 * @1", RooArgList(*nSgnCorr1, *nBkgCorr2));
  RooFormulaVar nBkgReflCorr("nBkgReflCorr", "Corrected fraction of bkg x refl", "@0 * @1", RooArgList(*nBkgCorr1, *nReflCorr2));
  RooFormulaVar nReflBkgCorr("nReflBkgCorr", "Corrected fraction of refl x bkg", "@0 * @1", RooArgList(*nReflCorr1, *nBkgCorr2));
  RooFormulaVar nSgnSgnCorr("nSgnSgnCorr", "Corrected fraction of sgn x sgn", "@0 * @1", RooArgList(*nSgnCorr1, *nSgnCorr2));
  RooFormulaVar nReflReflCorr("nReflReflCorr", "Corrected fraction of refl x refl", "@0 * @1", RooArgList(*nReflCorr1, *nReflCorr2));
  RooFormulaVar nSgnReflCorr("nSgnReflCorr", "Corrected fraction of sgn x refl", "@0 * @1", RooArgList(*nSgnCorr1, *nReflCorr2));
  RooFormulaVar nReflSgnCorr("nReflSgnCorr", "Corrected fraction of refl x sgn", "@0 * @1", RooArgList(*nReflCorr1, *nSgnCorr2));

  RooProdPdf sgnSgnFunc2D("sgnSgnFunc2D", "sgn * sgn 2D term", RooArgList(*sgnPdfCand1, *sgnPdfCand2));
  RooProdPdf bkgBkgFunc2D("bkgBkgFunc2D", "bkg * bkg 2D term", RooArgList(*bkgPdfCand1, *bkgPdfCand2));
  RooProdPdf sgnBkgFunc2D("sgnBkgFunc2D", "sgn * bkg cross-term", RooArgList(*sgnPdfCand1, *bkgPdfCand2));
  RooProdPdf bkgSgnFunc2D("bkgSgnFunc2D", "bkg * sgn cross-term", RooArgList(*bkgPdfCand1, *sgnPdfCand2));

  RooProdPdf sgnReflFunc2D("sgnReflFunc2D", "sgn * refl 2D term", RooArgList(*sgnPdfCand1, *reflPdfCand2));
  RooProdPdf bkgReflFunc2D("bkgReflFunc2D", "bkg * refl 2D term", RooArgList(*bkgPdfCand1, *reflPdfCand2));
  RooProdPdf reflBkgFunc2D("reflBkgFunc2D", "refl * bkg cross-term", RooArgList(*reflPdfCand1, *bkgPdfCand2));
  RooProdPdf reflSgnFunc2D("reflSgnFunc2D", "refl * sgn cross-term", RooArgList(*reflPdfCand1, *sgnPdfCand2));
  RooProdPdf reflReflFunc2D("reflReflFunc2D", "refl * refl cross-term", RooArgList(*reflPdfCand1, *reflPdfCand2));

  RooProdPdf sgnSgnFunc2DCorr("sgnSgnFunc2DCorr", "sgn * sgn 2D term", RooArgList(*sgnPdfCand1Corr, *sgnPdfCand2Corr));
  RooProdPdf bkgBkgFunc2DCorr("bkgBkgFunc2DCorr", "bkg * bkg 2D term", RooArgList(*bkgPdfCand1Corr, *bkgPdfCand2Corr));
  RooProdPdf sgnBkgFunc2DCorr("sgnBkgFunc2DCorr", "sgn * bkg cross-term", RooArgList(*sgnPdfCand1Corr, *bkgPdfCand2Corr));
  RooProdPdf bkgSgnFunc2DCorr("bkgSgnFunc2DCorr", "bkg * sgn cross-term", RooArgList(*bkgPdfCand1Corr, *sgnPdfCand2Corr));

  RooProdPdf sgnReflFunc2DCorr("sgnReflFunc2DCorr", "sgn * refl 2D term", RooArgList(*sgnPdfCand1Corr, *reflPdfCand2Corr));
  RooProdPdf bkgReflFunc2DCorr("bkgReflFunc2DCorr", "bkg * refl 2D term", RooArgList(*bkgPdfCand1Corr, *reflPdfCand2Corr));
  RooProdPdf reflBkgFunc2DCorr("reflBkgFunc2DCorr", "refl * bkg cross-term", RooArgList(*reflPdfCand1Corr, *bkgPdfCand2Corr));
  RooProdPdf reflSgnFunc2DCorr("reflSgnFunc2DCorr", "refl * sgn cross-term", RooArgList(*reflPdfCand1Corr, *sgnPdfCand2Corr));
  RooProdPdf reflReflFunc2DCorr("reflReflFunc2DCorr", "refl * refl cross-term", RooArgList(*reflPdfCand1Corr, *reflPdfCand2Corr));

  // There isn't a direct equivalent of RooPlot for >1 dimensions
  TH2D *hMassCorrelations = new TH2D("hMassCorrelations", "data of mass correlations", 40, _massMin, _massMax, 40, _massMin, _massMax);
  dataset->fillHistogram(hMassCorrelations, RooArgList(*massCand1, *massCand2));

  TH2D *hMassWeightedCorrelations = new TH2D("hMassWeightedCorrelations", "data of weighted mass correlations", 40, _massMin, _massMax, 40, _massMin, _massMax);
  weightedDataset->fillHistogram(hMassWeightedCorrelations, RooArgList(*massCand1, *massCand2));

  // Define the difference between nSgn1 and nSgn2
  RooFormulaVar diffSgn("diffSgn", "difference between nSgn1 and nSgn2", "@0 - @1", RooArgList(*nSgn1, *nSgn2));
  RooFormulaVar diffBkg("diffBkg", "difference between nBkg1 and nBkg2", "@0 - @1", RooArgList(*nBkg1, *nBkg2));
  RooFormulaVar diffSgnCorr("diffSgnCorr", "difference between nSgnCorr1 and nSgnCorr2", "@0 - @1", RooArgList(*nSgnCorr1, *nSgnCorr2));
  RooFormulaVar diffBkgCorr("diffBkgCorr", "difference between nBkgCorr1 and nBkgCorr2", "@0 - @1", RooArgList(*nBkgCorr1, *nBkgCorr2));

  // Define the Gaussian constraint to keep nSgn1 and nSgn2 symmetric
  RooRealVar sigmaConstraint("sigmaConstraint", "width of Gaussian constraint", 10.0, 1.0, 100000.0); // adjust this as needed
  sigmaConstraint.setConstant(kTRUE);
  RooRealVar sigmaConstraintCorr("sigmaConstraintCorr", "width of Gaussian constraint corr", 500.0, 1.0, 10000.0); // adjust this as needed
  sigmaConstraintCorr.setConstant(kTRUE);
  RooGaussian *constraintSgn = new RooGaussian("constraintSgn", "Gaussian constraint on nSgn1 - nSgn2", diffSgn, 0, sigmaConstraint);
  RooGaussian *constraintSgnCorr = new RooGaussian("constraintSgnCorr", "Gaussian constraint on nSgnCorr1 - nSgnCorr2", diffSgnCorr, 0, sigmaConstraintCorr);
  RooGaussian *constraintBkgCorr = new RooGaussian("constraintBkgCorr", "Gaussian constraint on nBkgCorr1 - nBkgCorr2", diffBkgCorr, 0, sigmaConstraintCorr);
  // RooArgSet constraintsCorr(*constraintSgnCorr, *constraintBkgCorr);

  // Final fit
  if (isMc) {
    if (doReflections) {
      _totPdf2D = new RooAddPdf("_totPdf2D", "signal + reflection pdf 2D",
        RooArgList(sgnSgnFunc2D, sgnReflFunc2D, reflSgnFunc2D, sgnBkgFunc2D,
                   bkgSgnFunc2D, bkgBkgFunc2D, bkgReflFunc2D, reflBkgFunc2D, reflReflFunc2D),
        RooArgList(nSgnSgn, nSgnRefl, nReflSgn, nSgnBkg, nBkgSgn, nBkgBkg, nBkgRefl, nReflBkg, nReflRefl));

      _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "Corrected signal + reflection pdf 2D",
                RooArgList(sgnSgnFunc2D, sgnReflFunc2D, reflSgnFunc2D, reflReflFunc2D),
                RooArgList(nSgnSgnCorr, nSgnReflCorr, nReflSgnCorr, nReflReflCorr));
    } else {
      _totPdf2D = new RooAddPdf("_totPdf2D", "signal pdf 2D",
                                RooArgList(sgnSgnFunc2D),
                                RooArgList(nSgnSgn));

      _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "Corrected signal pdf 2D",
                                        RooArgList(sgnSgnFunc2DCorr),
                                        RooArgList(nSgnSgnCorr));
    }
  } else {
    if (doReflections) {
      _totPdf2D = new RooAddPdf("_totPdf2D", "background + signal + reflection pdf 2D",
                                RooArgList(sgnSgnFunc2D, sgnReflFunc2D, reflSgnFunc2D, sgnBkgFunc2D,
                                           bkgSgnFunc2D, bkgBkgFunc2D, bkgReflFunc2D, reflBkgFunc2D, reflReflFunc2D),
                                RooArgList(nSgnSgn, nSgnRefl, nReflSgn, nSgnBkg, nBkgSgn, nBkgBkg, nBkgRefl, nReflBkg, nReflRefl));

      _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "Corrected background + signal + reflection pdf 2D",
                                        RooArgList(sgnSgnFunc2DCorr, sgnReflFunc2DCorr, reflSgnFunc2DCorr, sgnBkgFunc2DCorr, bkgSgnFunc2DCorr,
                                                   bkgBkgFunc2DCorr, bkgReflFunc2DCorr, reflBkgFunc2DCorr, reflReflFunc2DCorr),
                                        RooArgList(nSgnSgnCorr, nSgnReflCorr, nReflSgnCorr, nSgnBkgCorr, nBkgSgnCorr, nBkgBkgCorr,
                                                   nBkgReflCorr, nReflBkgCorr, nReflReflCorr));
    } else {
      _totPdf2D = new RooAddPdf("_totPdf2D", "background + signal pdf 2D",
                                RooArgList(bkgBkgFunc2D, sgnSgnFunc2D, sgnBkgFunc2D, bkgSgnFunc2D),
                                RooArgList(nBkgBkg, nSgnSgn, nSgnBkg, nBkgSgn));

      _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "background + signal pdf 2D",
                                        RooArgList(bkgBkgFunc2DCorr, sgnSgnFunc2DCorr, sgnBkgFunc2DCorr, bkgSgnFunc2DCorr),
                                        RooArgList(nBkgBkgCorr, nSgnSgnCorr, nSgnBkgCorr, nBkgSgnCorr));
    }
  }

  _totPdf2D->fixCoefNormalization(normSet);
  _weightedTotPdf2D->fixCoefNormalization(weightedNormSet);

  // Include the constraint in the total model
  RooProdPdf *modelWithConstraint = new RooProdPdf("modelWithConstraint", "model with signal symmetry constraint",
                                                   RooArgList(*_totPdf2D, *constraintSgn));

  RooProdPdf *modelWithConstraintCorr = new RooProdPdf("modelWithConstraintCorr", "model with signal symmetry constraint Corr",
                                                       RooArgList(*_weightedTotPdf2D, *constraintSgnCorr));

  // Perform the fit
  RooFitResult *resultConstrained = modelWithConstraint->fitTo(*dataset, RooFit::Save(),
                                                               //RooFit::SumW2Error(kTRUE),
                                                               RooFit::Minimizer("Minuit2", "Migrad"),
                                                               RooFit::Range(_massMin, _massMax),
                                                               RooFit::Strategy(1),
                                                               //RooFit::Optimize(2),
                                                               RooFit::MaxCalls(1000), // Adjust max calls as needed
                                                               RooFit::PrintLevel(-1)   // Suppress verbose output
  );
  resultConstrained->correlationMatrix().Print();
  if (resultConstrained->status() == 0) {
    std::cout << "Fit Converged!" << std::endl;
  } else if (resultConstrained->status() == 1) {
    std::cout << "Fit did not fully converge (status=1), but covariance matrix is OK." << std::endl;
  } else {
    std::cout << "Fit failed (status=" << resultConstrained->status() << ")." << std::endl;
  }
  resultConstrained->Print("v");

  const RooArgList& paramsRawfit = resultConstrained->floatParsFinal();

  RooFitResult *resultConstrainedCorr = modelWithConstraintCorr->fitTo(*weightedDataset, RooFit::Save(),
                                                                       RooFit::SumW2Error(kTRUE),
                                                                       RooFit::Minimizer("Minuit2", "Migrad"),
                                                                       RooFit::Range(_massMin, _massMax),
                                                                       RooFit::Strategy(1),
                                                                       //RooFit::Optimize(2),
                                                                       RooFit::MaxCalls(1000), // Adjust max calls as needed
                                                                       RooFit::PrintLevel(-1)   // Suppress verbose output
  );

  if (resultConstrainedCorr->status() == 0) {
    std::cout << "Fit Converged!" << std::endl;
  } else if (resultConstrainedCorr->status() == 1) {
    std::cout << "Fit did not fully converge (status=1), but covariance matrix is OK." << std::endl;
  } else {
    std::cout << "Fit failed (status=" << resultConstrainedCorr->status() << ")." << std::endl;
  }

  resultConstrainedCorr->Print("v");

  // -------------------------------------
  // -- Get raw yields per contribution --
  // -------------------------------------

  // Set the range for massCand1 and massCand2
  massCand1->setRange("3sigmaRange", mean->getVal() - 3 * sigma->getVal(), mean->getVal() + 3 * sigma->getVal());
  massCand2->setRange("3sigmaRange", mean->getVal() - 3 * sigma->getVal(), mean->getVal() + 3 * sigma->getVal());

  // Assign the final fitted parameters to the model with constraints
  modelWithConstraint->getParameters(*massCand1)->assignValueOnly(resultConstrained->floatParsFinal());
  modelWithConstraint->getParameters(*massCand2)->assignValueOnly(resultConstrained->floatParsFinal());

  cout << "\n\n ||||| RAW INTEGRAL VALUES |||||" << endl;
  cout << "Signal x Signal ";
  RooRealVar *yieldSgnSgn = getYieldInRange(resultConstrained, massCand1, massCand2, sgnSgnFunc2D, nSgnSgn, "3sigmaRange");

  cout << "Signal x Background ";
  RooRealVar *yieldSgnBkg = getYieldInRange(resultConstrained, massCand1, massCand2, sgnBkgFunc2D, nSgnBkg, "3sigmaRange");

  cout << "Background x Signal ";
  RooRealVar *yieldBkgSgn = getYieldInRange(resultConstrained, massCand1, massCand2, bkgSgnFunc2D, nBkgSgn, "3sigmaRange");

  cout << "Background x Background ";
  RooRealVar *yieldBkgBkg = getYieldInRange(resultConstrained, massCand1, massCand2, bkgBkgFunc2D, nBkgBkg, "3sigmaRange");

  RooRealVar *yieldSgnRefl = nullptr;
  RooRealVar *yieldReflSgn = nullptr;
  RooRealVar *yieldBkgRefl = nullptr;
  RooRealVar *yieldReflBkg = nullptr;
  RooRealVar *yieldReflRefl = nullptr;
  if (doReflections)
  {
    cout << "Signal x Reflected ";
    yieldSgnRefl = getYieldInRange(resultConstrained, massCand1, massCand2, sgnReflFunc2D, nSgnRefl, "3sigmaRange");

    cout << "Reflected x Signal ";
    yieldReflSgn = getYieldInRange(resultConstrained, massCand1, massCand2, reflSgnFunc2D, nReflSgn, "3sigmaRange");

    cout << "Background x Reflected ";
    yieldBkgRefl = getYieldInRange(resultConstrained, massCand1, massCand2, bkgReflFunc2D, nBkgRefl, "3sigmaRange");

    cout << "Reflected x Background ";
    yieldReflBkg = getYieldInRange(resultConstrained, massCand1, massCand2, reflBkgFunc2D, nReflBkg, "3sigmaRange");

    cout << "Reflected x Reflected ";
    yieldReflRefl = getYieldInRange(resultConstrained, massCand1, massCand2, reflReflFunc2D, nReflRefl, "3sigmaRange");
  }
  // Display results
  cout << "\n\n ||||| RAW RESULTS |||||" << endl;
  cout << "Number of SgnSgn candidates " << nSgnSgn.getVal() << " +- " << nSgnSgn.getPropagatedError(*resultConstrained) << endl;
  cout << "Number of SgnBkg candidates " << nSgnBkg.getVal() << " +- " << nSgnBkg.getPropagatedError(*resultConstrained) << endl;
  cout << "Number of BkgSgn candidates " << nBkgSgn.getVal() << " +- " << nBkgSgn.getPropagatedError(*resultConstrained) << endl;
  cout << "Number of BkgBkg candidates " << nBkgBkg.getVal() << " +- " << nBkgBkg.getPropagatedError(*resultConstrained) << endl;
  if (doReflections)
  {
    cout << "Number of ReflSgn candidates " << nReflSgn.getVal() << " +- " << nReflSgn.getPropagatedError(*resultConstrained) << endl;
    cout << "Number of SgnRefl candidates " << nSgnRefl.getVal() << " +- " << nSgnRefl.getPropagatedError(*resultConstrained) << endl;
    cout << "Number of BkgRefl candidates " << nBkgRefl.getVal() << " +- " << nBkgRefl.getPropagatedError(*resultConstrained) << endl;
    cout << "Number of ReflBkg candidates " << nReflBkg.getVal() << " +- " << nReflBkg.getPropagatedError(*resultConstrained) << endl;
    cout << "Number of ReflRefl candidates " << nReflRefl.getVal() << " +- " << nReflRefl.getPropagatedError(*resultConstrained) << endl;
  }
  cout << "\n ---------------------------------------------- " << endl;
  cout << " ---------------------------------------------- \n"
       << endl;

  cout << " SgnSgn yield " << yieldSgnSgn->getVal() << " +- " << yieldSgnSgn->getError() << endl;
  cout << " SgnBkg yield " << yieldSgnBkg->getVal() << " +- " << yieldSgnBkg->getError() << endl;
  cout << " BkgSgn yield " << yieldBkgSgn->getVal() << " +- " << yieldBkgSgn->getError() << endl;
  cout << " BkgBkg yield " << yieldBkgBkg->getVal() << " +- " << yieldBkgBkg->getError() << endl;

  if (doReflections)
  {
    cout << " SgnRefl yield " << yieldSgnRefl->getVal() << " +- " << yieldSgnRefl->getError() << endl;
    cout << " ReflSgn yield " << yieldReflSgn->getVal() << " +- " << yieldReflSgn->getError() << endl;
    cout << " ReflBkg yield " << yieldBkgRefl->getVal() << " +- " << yieldBkgRefl->getError() << endl;
    cout << " BkgRefl yield " << yieldReflBkg->getVal() << " +- " << yieldReflBkg->getError() << endl;
    cout << " ReflRefl yield " << yieldReflRefl->getVal() << " +- " << yieldReflRefl->getError() << endl;
  }
  cout << "\n ---------------------------------------------- " << endl;
  cout << " ---------------------------------------------- \n\n\n"
       << endl;
 cout << "joder" << endl;
  // Save these results in a histogram
  TH1D *hYields = new TH1D("hYields", "Raw Yields", 4, 0, 4); // 4 bins for SgnSgn, SgnBkg, BkgSgn, BkgBkg
  hYields->GetXaxis()->SetBinLabel(1, "SgnSgn");
  hYields->GetXaxis()->SetBinLabel(2, "SgnBkg");
  hYields->GetXaxis()->SetBinLabel(3, "BkgSgn");
  hYields->GetXaxis()->SetBinLabel(4, "BkgBkg");

  hYields->SetBinContent(1, yieldSgnSgn->getVal());
  hYields->SetBinError(1, yieldSgnSgn->getError());
  hYields->SetBinContent(2, yieldSgnBkg->getVal());
  hYields->SetBinError(2, yieldSgnBkg->getError());
  hYields->SetBinContent(3, yieldBkgSgn->getVal());
  hYields->SetBinError(3, yieldBkgSgn->getError());
  hYields->SetBinContent(4, yieldBkgBkg->getVal());
  hYields->SetBinError(4, yieldBkgBkg->getError());
  if (doReflections)
  {
    hYields->SetBins(9, 0, 9);
    hYields->GetXaxis()->SetBinLabel(5, "SgnRefl");
    hYields->GetXaxis()->SetBinLabel(6, "ReflSgn");
    hYields->GetXaxis()->SetBinLabel(7, "BkgRefl");
    hYields->GetXaxis()->SetBinLabel(8, "ReflBkg");
    hYields->GetXaxis()->SetBinLabel(9, "ReflRefl");

    hYields->SetBinContent(5, yieldSgnRefl->getVal());
    hYields->SetBinError(5, yieldSgnRefl->getError());
    hYields->SetBinContent(6, yieldReflSgn->getVal());
    hYields->SetBinError(6, yieldReflSgn->getError());
    hYields->SetBinContent(7, yieldBkgRefl->getVal());
    hYields->SetBinError(7, yieldBkgRefl->getError());
    hYields->SetBinContent(8, yieldReflBkg->getVal());
    hYields->SetBinError(8, yieldReflBkg->getError());
    hYields->SetBinContent(9, yieldReflRefl->getVal());
    hYields->SetBinError(9, yieldReflRefl->getError());
  }

  // ------------------------------------------------------
  // -- Get efficiency corrected yields per contribution --
  // ------------------------------------------------------

  modelWithConstraintCorr->getParameters(*massCand1)->assignValueOnly(resultConstrainedCorr->floatParsFinal());
  modelWithConstraintCorr->getParameters(*massCand2)->assignValueOnly(resultConstrainedCorr->floatParsFinal());

  cout << "\n\n ||||| EFFICIENCY CORRECTED INTEGRAL VALUES |||||" << endl;
  cout << "Signal x Signal ";
  RooRealVar *yieldSgnSgnCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, sgnSgnFunc2DCorr, nSgnSgnCorr, "3sigmaRange");

  cout << "Signal x Background ";
  RooRealVar *yieldSgnBkgCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, sgnBkgFunc2DCorr, nSgnBkgCorr, "3sigmaRange");

  cout << "Background x Signal ";
  RooRealVar *yieldBkgSgnCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, bkgSgnFunc2DCorr, nBkgSgnCorr, "3sigmaRange");

  cout << "Background x Background ";
  RooRealVar *yieldBkgBkgCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, bkgBkgFunc2DCorr, nBkgBkgCorr, "3sigmaRange");

  RooRealVar *yieldSgnReflCorr = nullptr;
  RooRealVar *yieldReflSgnCorr = nullptr;
  RooRealVar *yieldBkgReflCorr = nullptr;
  RooRealVar *yieldReflBkgCorr = nullptr;
  RooRealVar *yieldReflReflCorr = nullptr;
  if (doReflections)
  {

    cout << "Signal x Reflected ";
    yieldSgnReflCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, sgnReflFunc2DCorr, nSgnReflCorr, "3sigmaRange");

    cout << "Reflected x Signal ";
    yieldReflSgnCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, reflSgnFunc2DCorr, nReflSgnCorr, "3sigmaRange");

    cout << "Background x Reflected ";
    yieldBkgReflCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, bkgReflFunc2DCorr, nBkgReflCorr, "3sigmaRange");

    cout << "Reflected x Background ";
    yieldReflBkgCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, reflBkgFunc2DCorr, nReflBkgCorr, "3sigmaRange");

    cout << "Reflected x Reflected ";
    yieldReflReflCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, reflReflFunc2DCorr, nReflReflCorr, "3sigmaRange");
  }
  // Display results
  cout << "\n\n ||||| EFFICIENCY CORRECTED RESULTS ||||| " << endl;
  cout << "Number of SgnSgn candidates " << nSgnSgnCorr.getVal() << " +- " << nSgnSgnCorr.getPropagatedError(*resultConstrainedCorr) << endl;
  cout << "Number of SgnBkg candidates " << nSgnBkgCorr.getVal() << " +- " << nSgnBkgCorr.getPropagatedError(*resultConstrainedCorr) << endl;
  cout << "Number of BkgSgn candidates " << nBkgSgnCorr.getVal() << " +- " << nBkgSgnCorr.getPropagatedError(*resultConstrainedCorr) << endl;
  cout << "Number of BkgBkg candidates " << nBkgBkgCorr.getVal() << " +- " << nBkgBkgCorr.getPropagatedError(*resultConstrainedCorr) << endl;
  if (doReflections)
  {
    cout << "Number of ReflSgn candidates " << nReflSgnCorr.getVal() << " +- " << nReflSgnCorr.getPropagatedError(*resultConstrainedCorr) << endl;
    cout << "Number of SgnRefl candidates " << nSgnReflCorr.getVal() << " +- " << nSgnReflCorr.getPropagatedError(*resultConstrainedCorr) << endl;
    cout << "Number of BkgRefl candidates " << nBkgReflCorr.getVal() << " +- " << nBkgReflCorr.getPropagatedError(*resultConstrainedCorr) << endl;
    cout << "Number of ReflBkg candidates " << nReflBkgCorr.getVal() << " +- " << nReflBkgCorr.getPropagatedError(*resultConstrainedCorr) << endl;
    cout << "Number of ReflRefl candidates " << nReflReflCorr.getVal() << " +- " << nReflReflCorr.getPropagatedError(*resultConstrainedCorr) << endl;
  }
  cout << " ---------------------------------------------- " << endl;
  cout << " ---------------------------------------------- \n"
       << endl;

  cout << " SgnSgn yield " << yieldSgnSgnCorr->getVal() << " +- " << yieldSgnSgnCorr->getError() << endl;
  cout << " SgnBkg yield " << yieldSgnBkgCorr->getVal() << " +- " << yieldSgnBkgCorr->getError() << endl;
  cout << " BkgSgn yield " << yieldBkgSgnCorr->getVal() << " +- " << yieldBkgSgnCorr->getError() << endl;
  cout << " BkgBkg yield " << yieldBkgBkgCorr->getVal() << " +- " << yieldBkgBkgCorr->getError() << endl;
  if (doReflections)
  {
    cout << " SgnRefl yield " << yieldSgnReflCorr->getVal() << " +- " << yieldSgnReflCorr->getError() << endl;
    cout << " ReflSgn yield " << yieldReflSgnCorr->getVal() << " +- " << yieldReflSgnCorr->getError() << endl;
    cout << " ReflBkg yield " << yieldBkgReflCorr->getVal() << " +- " << yieldBkgReflCorr->getError() << endl;
    cout << " BkgRefl yield " << yieldReflBkgCorr->getVal() << " +- " << yieldReflBkgCorr->getError() << endl;
    cout << " ReflRefl yield " << yieldReflReflCorr->getVal() << " +- " << yieldReflReflCorr->getError() << endl;
  }
  cout << " ---------------------------------------------- " << endl;
  cout << " ---------------------------------------------- \n\n\n"
       << endl;

  // Save the results in a histogram
  TH1D *hYieldsCorr = new TH1D("hYieldsCorr", "Raw Yields Corrected by efficiencies", 4, 0, 4); // 4 bins for SgnSgn, SgnBkg, BkgSgn, BkgBkg
  hYieldsCorr->GetXaxis()->SetBinLabel(1, "SgnSgn");
  hYieldsCorr->GetXaxis()->SetBinLabel(2, "SgnBkg");
  hYieldsCorr->GetXaxis()->SetBinLabel(3, "BkgSgn");
  hYieldsCorr->GetXaxis()->SetBinLabel(4, "BkgBkg");

  hYieldsCorr->SetBinContent(1, yieldSgnSgnCorr->getVal());
  hYieldsCorr->SetBinError(1, yieldSgnSgnCorr->getError());
  hYieldsCorr->SetBinContent(2, yieldSgnBkgCorr->getVal());
  hYieldsCorr->SetBinError(2, yieldSgnBkgCorr->getError());
  hYieldsCorr->SetBinContent(3, yieldBkgSgnCorr->getVal());
  hYieldsCorr->SetBinError(3, yieldBkgSgnCorr->getError());
  hYieldsCorr->SetBinContent(4, yieldBkgBkgCorr->getVal());
  hYieldsCorr->SetBinError(4, yieldBkgBkgCorr->getError());

  if (doReflections)
  {
    hYieldsCorr->SetBins(9, 0, 9);
    hYieldsCorr->GetXaxis()->SetBinLabel(5, "SgnRefl");
    hYieldsCorr->GetXaxis()->SetBinLabel(6, "ReflSgn");
    hYieldsCorr->GetXaxis()->SetBinLabel(7, "BkgRefl");
    hYieldsCorr->GetXaxis()->SetBinLabel(8, "ReflBkg");
    hYieldsCorr->GetXaxis()->SetBinLabel(9, "ReflRefl");

    hYieldsCorr->SetBinContent(5, yieldSgnReflCorr->getVal());
    hYieldsCorr->SetBinError(5, yieldSgnReflCorr->getError());
    hYieldsCorr->SetBinContent(6, yieldReflSgnCorr->getVal());
    hYieldsCorr->SetBinError(6, yieldReflSgnCorr->getError());
    hYieldsCorr->SetBinContent(7, yieldBkgReflCorr->getVal());
    hYieldsCorr->SetBinError(7, yieldBkgReflCorr->getError());
    hYieldsCorr->SetBinContent(8, yieldReflBkgCorr->getVal());
    hYieldsCorr->SetBinError(8, yieldReflBkgCorr->getError());
    hYieldsCorr->SetBinContent(9, yieldReflReflCorr->getVal());
    hYieldsCorr->SetBinError(9, yieldReflReflCorr->getError());
  }

  TH1F* hWidth = new TH1F("hWidth", "Width", 1, 0, 1);
  hWidth->SetBinContent(1, sigma->getVal());
  hWidth->SetBinError(1, sigma->getError());

  fout->cd();
  hWidth->Write();
  hYields->Write();
  hYieldsCorr->Write();

  // Check by dividing each component by the integrated eff^2
  cout << " ---------------------------------------------- " << endl;
  cout << " ---------------------------------------------- \n"
       << endl;

  cout << " NUMBER OF CANDIDATES / INTEGRATED EFF^2" << endl;

  cout << " SgnSgn yield " << yieldSgnSgn->getVal() / (_integratedEfficiency * _integratedEfficiency) << " +- " << yieldSgnSgn->getError() / (_integratedEfficiency * _integratedEfficiency) << endl;
  cout << " SgnBkg yield " << yieldSgnBkg->getVal() / (_integratedEfficiency * _integratedEfficiency) << " +- " << yieldSgnBkg->getError() / (_integratedEfficiency * _integratedEfficiency) << endl;
  cout << " BkgSgn yield " << yieldBkgSgn->getVal() / (_integratedEfficiency * _integratedEfficiency) << " +- " << yieldBkgSgn->getError() / (_integratedEfficiency * _integratedEfficiency) << endl;
  cout << " BkgBkg yield " << yieldBkgBkg->getVal() / (_integratedEfficiency * _integratedEfficiency) << " +- " << yieldBkgBkg->getError() / (_integratedEfficiency * _integratedEfficiency) << endl;
  if (doReflections)
  {
    cout << " SgnRefl yield " << yieldSgnRefl->getVal() / (_integratedEfficiency * _integratedEfficiency) << " +- " << yieldSgnRefl->getError() / (_integratedEfficiency * _integratedEfficiency) << endl;
    cout << " ReflSgn yield " << yieldReflSgn->getVal() / (_integratedEfficiency * _integratedEfficiency) << " +- " << yieldReflSgn->getError() / (_integratedEfficiency * _integratedEfficiency) << endl;
    cout << " ReflBkg yield " << yieldBkgRefl->getVal() / (_integratedEfficiency * _integratedEfficiency) << " +- " << yieldBkgRefl->getError() / (_integratedEfficiency * _integratedEfficiency) << endl;
    cout << " BkgRefl yield " << yieldReflBkg->getVal() / (_integratedEfficiency * _integratedEfficiency) << " +- " << yieldReflBkg->getError() / (_integratedEfficiency * _integratedEfficiency) << endl;
    cout << " ReflRefl yield " << yieldReflRefl->getVal() / (_integratedEfficiency * _integratedEfficiency) << " +- " << yieldReflRefl->getError() / (_integratedEfficiency * _integratedEfficiency) << endl;
  }
  cout << " ---------------------------------------------- " << endl;
  cout << " ---------------------------------------------- \n\n\n"
       << endl;

  // Assuming your fit PDF is called pdf2D and it is defined with RooRealVars x and y
  double xVal = 1.85; // Your desired x coordinate
  double yVal = 1.85; // Your desired y coordinate

  // Retrieve the list of floating parameters from the fit result
  RooArgSet *params = modelWithConstraint->getParameters(RooArgSet(*massCand1, *massCand2));
  params->assignValueOnly(resultConstrained->floatParsFinal()); // Synchronize values

  // Set x and y to the specific point of interest
  massCand1->setVal(xVal);
  massCand2->setVal(yVal);

  // Evaluate the PDF value at (x, y)
  double value = modelWithConstraint->getVal(RooArgSet(*massCand1, *massCand2));

  // Get the propagated error at this point
  double error = modelWithConstraint->getPropagatedError(*resultConstrained, RooArgSet(*massCand1, *massCand2));

  // Calculate the relative error
  double relativeError = (value != 0) ? error / value : 0;

  std::cout << "Value at (" << xVal << ", " << yVal << "): " << value << std::endl;
  std::cout << "Error at (" << xVal << ", " << yVal << "): " << error << std::endl;
  std::cout << "Relative Error: " << relativeError << std::endl;

  RooArgSet *paramsCorr = modelWithConstraintCorr->getParameters(RooArgSet(*massCand1, *massCand2));
  paramsCorr->assignValueOnly(resultConstrainedCorr->floatParsFinal()); // Synchronize values

  double valueCorr = modelWithConstraintCorr->getVal(RooArgSet(*massCand1, *massCand2));
  double errorCorr = modelWithConstraintCorr->getPropagatedError(*resultConstrainedCorr, RooArgSet(*massCand1, *massCand2));
  double relativeErrorCorr = (value != 0) ? errorCorr / valueCorr : 0;

  std::cout << "ValueCorr at (" << xVal << ", " << yVal << "): " << valueCorr << std::endl;
  std::cout << "ErrorCorr at (" << xVal << ", " << yVal << "): " << errorCorr << std::endl;
  std::cout << "Relative Error Corr: " << relativeErrorCorr << std::endl;
  cout << "\n\n\n";

  // Fill histogram with fit results
  double binWidth = 0.005; // Desired bin width
  int nBins = (int)((_massMax - _massMin) / binWidth);

  // Set the calculated number of bins
  // massCand1->setBins(nBins);
  // massCand2->setBins(nBins);

  TH2D *histFit = new TH2D("histFit", "Fit Results;m(D_{1}^{0}) (GeV/#it{c}^2);#it{M}(D_{2}^{0}) (GeV/#it{c}^2)",
                           nBins, _massMin, _massMax, nBins, _massMin, _massMax);

  TH2D *histFitWeighted = new TH2D("histFitWeighted", "Efficiency corrected fit Results;m(D_{1}^{0}) (GeV/#it{c}^2);#it{M}(D_{2}^{0}) (GeV/#it{c}^2)",
                                   nBins, _massMin, _massMax, nBins, _massMin, _massMax);

  plot2DFit(hMassCorrelations, histFit, modelWithConstraint, draw, fout, "cRawMassCorrelations");
  plot2DFit(hMassWeightedCorrelations, histFitWeighted, modelWithConstraintCorr, draw, fout, "cWeightedMassCorrelations");

  cout << " Integral of hMassCorrelations = " << hMassCorrelations->Integral() << endl;
  cout << " Integral of hMassWeightedCorrelations = " << hMassWeightedCorrelations->Integral() << endl;

  if (draw)
  {
    // Plot projections
    plotProjectionsAfterFit(resultConstrained, modelWithConstraint, dataset, "cRawMassProjections", fout, doReflections, "");
    plotProjectionsAfterFit(resultConstrainedCorr, modelWithConstraintCorr, weightedDataset, "cWeightedMassProjections", fout, doReflections, "Corr");
  }

  // ----------------------------------------------------------
  // -- Make checks comparing raw and weighted distributions --
  // ----------------------------------------------------------

  // Normalise by dividing by the total number of candidates in a 3sgima range from the peak
  // Check: ratio of histFitScaled / histFitWeightdScaled (using fit functions)

  double integral3SigmaHistFit = histFit->Integral(
      histFit->GetXaxis()->FindBin(mean->getVal() - 3 * sigma->getVal()),
      histFit->GetXaxis()->FindBin(mean->getVal() + 3 * sigma->getVal()),
      histFit->GetYaxis()->FindBin(mean->getVal() - 3 * sigma->getVal()),
      histFit->GetYaxis()->FindBin(mean->getVal() + 3 * sigma->getVal()));

  double integral3SigmaHistFitWeighted = histFitWeighted->Integral(
      histFitWeighted->GetXaxis()->FindBin(mean->getVal() - 3 * sigma->getVal()),
      histFitWeighted->GetXaxis()->FindBin(mean->getVal() + 3 * sigma->getVal()),
      histFitWeighted->GetYaxis()->FindBin(mean->getVal() - 3 * sigma->getVal()),
      histFitWeighted->GetYaxis()->FindBin(mean->getVal() + 3 * sigma->getVal()));

  TH2D *histFitScaled = (TH2D *)histFit->Clone("histFitScaled");
  TH2D *histFitWeightedScaled = (TH2D *)histFitWeighted->Clone("histFitWeightedScaled");
  histFitScaled->Rebin2D(2, 2);
  histFitWeightedScaled->Rebin2D(2, 2);

  histFitScaled->Scale(1.0 / integral3SigmaHistFit);
  histFitWeightedScaled->Scale(1.0 / integral3SigmaHistFitWeighted);
  TH2D *ratioHist2D = (TH2D *)histFitWeightedScaled->Clone("ratioHist2D");
  ratioHist2D->Divide(histFitScaled); // Perform bin-by-bin division

  // Check: Señal(raw) / AE(average)^2 = Señal(plot eff corrected)
  TH2D *histFitIntEff = (TH2D *)histFit->Clone("histFitIntEff");
  histFitIntEff->Scale(1.0 / (_integratedEfficiency * _integratedEfficiency));
  TH2D *ratioHist2DIntEff = (TH2D *)histFitIntEff->Clone("ratioHist2DIntEff");
  ratioHist2DIntEff->Divide(histFitWeighted);

  fout->cd();
  TCanvas *cRatioIntEff = new TCanvas("cRatioIntEff", "cRatioIntEff", 1000, 700);
  ratioHist2DIntEff->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
  ratioHist2DIntEff->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
  ratioHist2DIntEff->SetTitle("(raw yield/(integrated efficiency)^{2}) / efficiency-corrected distr.");
  ratioHist2DIntEff->Draw("colz");
  cRatioIntEff->Write();
  cRatioIntEff->Close();

  TCanvas *cRatioFunctions = new TCanvas("cRatioFunctions", "cRatioFunctions", 1000, 700);
  ratioHist2D->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
  ratioHist2D->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
  ratioHist2D->SetTitle("Normalised Efficiency-corrected / non-corrected fit functions");
  ratioHist2D->Draw("colz");
  cRatioFunctions->Write();
  cRatioFunctions->Close();

  // Normalise by dividing by the total number of candidates in a 3sigma range from the peak
  // Check: ratio of hMassScaled / hMassWeightdScaled (using datapoints this time)

  double integral3SigmaHistCorr = hMassCorrelations->Integral(
      hMassCorrelations->GetXaxis()->FindBin(mean->getVal() - 3 * sigma->getVal()),
      hMassCorrelations->GetXaxis()->FindBin(mean->getVal() + 3 * sigma->getVal()),
      hMassCorrelations->GetYaxis()->FindBin(mean->getVal() - 3 * sigma->getVal()),
      hMassCorrelations->GetYaxis()->FindBin(mean->getVal() + 3 * sigma->getVal()));

  double integral3SigmaHistCorrWeighted = hMassWeightedCorrelations->Integral(
      hMassWeightedCorrelations->GetXaxis()->FindBin(mean->getVal() - 3 * sigma->getVal()),
      hMassWeightedCorrelations->GetXaxis()->FindBin(mean->getVal() + 3 * sigma->getVal()),
      hMassWeightedCorrelations->GetYaxis()->FindBin(mean->getVal() - 3 * sigma->getVal()),
      hMassWeightedCorrelations->GetYaxis()->FindBin(mean->getVal() + 3 * sigma->getVal()));

  TH2D *hMassCorrelationsScaled = (TH2D *)hMassCorrelations->Clone("hMassCorrelationsScaled");
  TH2D *hMassWeightedCorrelationsScaled = (TH2D *)hMassWeightedCorrelations->Clone("hMassWeightedCorrelationsScaled");
  hMassCorrelationsScaled->Scale(1.0 / integral3SigmaHistCorr);
  hMassWeightedCorrelationsScaled->Scale(1.0 / integral3SigmaHistCorrWeighted);
  TH2D *ratioHistMass2D = (TH2D *)hMassWeightedCorrelationsScaled->Clone("ratioHistMass2D");
  ratioHistMass2D->Divide(hMassCorrelationsScaled); // Perform bin-by-bin division

  // Check: Señal(raw) / AE(average)^2 = Señal(plot eff corrected)
  TH2D *hMassCorrelationsIntEff = (TH2D *)hMassCorrelations->Clone("hMassWeightedCorrelationsIntEff");
  hMassCorrelationsIntEff->Scale(1.0 / (_integratedEfficiency * _integratedEfficiency));
  TH2D *ratioHistMass2DIntEff = (TH2D *)hMassCorrelationsIntEff->Clone("ratioHistMass2DIntEff");
  ratioHistMass2DIntEff->Divide(hMassWeightedCorrelations); // Perform bin-by-bin division

  TCanvas *cRatioMass = new TCanvas("cRatioMass", "cRatioMass", 1000, 700);
  ratioHistMass2D->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
  ratioHistMass2D->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
  ratioHistMass2D->SetTitle("Normalised Efficiency-corrected / non-corrected hMass");
  ratioHistMass2D->Draw("colz");
  cRatioMass->Write();
  cRatioMass->Close();

  TCanvas *cRatioMassIntEff = new TCanvas("cRatioMassIntEff", "cRatioMassIntEff", 1000, 700);
  ratioHistMass2DIntEff->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
  ratioHistMass2DIntEff->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
  ratioHistMass2DIntEff->SetTitle("(raw yield/(integrated efficiency)^{2}) / efficiency-corrected distr.");
  ratioHistMass2DIntEff->Draw("colz");
  cRatioMassIntEff->Write();
  cRatioMassIntEff->Close();

  // -----------------------------------------------------------------------
  // -- Plot pT distributions before and after correcting by efficiencies --
  // -----------------------------------------------------------------------


  resultConstrained->correlationMatrix().Print();
  if (resultConstrained->status() == 0)
  {
    std::cout << "Fit Converged!" << std::endl;
  }
  else
  {
    std::cout << "Fit did not converge." << std::endl;
  }
  resultConstrained->Print("v");

  resultConstrainedCorr->correlationMatrix().Print();
  if (resultConstrainedCorr->status() == 0)
  {
    std::cout << "Fit Converged!" << std::endl;
  }
  else
  {
    std::cout << "Fit did not converge." << std::endl;
  }
  resultConstrainedCorr->Print("v");
} // do2Dfit

void InvMassFitter2D::plot2DFit(TH2D *hMassCorrelations, TH2D *histFit, RooProdPdf *model, Bool_t draw, TFile *fout, TString const &cName)
{
  RooRealVar *massCand1 = _workspace.var("fMCand1");
  RooRealVar *massCand2 = _workspace.var("fMCand2");

  massCand1->setRange(_massMin, _massMax);
  massCand2->setRange(_massMin, _massMax);

  // Fill histogram with fit results
  Int_t binMin = hMassCorrelations->GetXaxis()->FindBin(_massMin);
  Int_t binMax = hMassCorrelations->GetXaxis()->FindBin(_massMax);
  for (Int_t binx = 1; binx <= histFit->GetNbinsX(); ++binx)
  {
    for (Int_t biny = 1; biny <= histFit->GetNbinsY(); ++biny)
    {
      Double_t xval = histFit->GetXaxis()->GetBinCenter(binx);
      Double_t yval = histFit->GetYaxis()->GetBinCenter(biny);
      massCand1->setVal(xval);
      massCand2->setVal(yval);
      histFit->SetBinContent(binx, biny, model->getVal());
    }
  }
  if (draw)
  {
    TCanvas *c2d = new TCanvas(cName, cName, 1000, 700);
    c2d->SetLeftMargin(0.1);
    c2d->SetRightMargin(0.1);
    c2d->SetTopMargin(0.2);
    c2d->SetBottomMargin(0.1);

    histFit->SetLineColor(kRed);
    histFit->SetLineWidth(1);
    histFit->SetTitle("");
    histFit->GetXaxis()->SetTitle("");
    histFit->GetYaxis()->SetTitle("");
    histFit->GetZaxis()->SetTitle("");

    //histFit->GetXaxis()->SetDecimals(2);
    //histFit->GetYaxis()->SetDecimals(2);
    //histFit->GetZaxis()->SetMaxDigits(3);  // Forces scientific notation for large/small numbers

    /* histFit->GetXaxis()->SetTicks(0);
    histFit->GetYaxis()->SetTicks(0); */
    //histFit->GetZaxis()->SetTicks(0);

    hMassCorrelations->GetXaxis()->SetDecimals(2);
    hMassCorrelations->GetYaxis()->SetDecimals(2);
    hMassCorrelations->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
    hMassCorrelations->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
    hMassCorrelations->GetZaxis()->SetTitle("Counts per (3.5 MeV/#it{c}^{2})^{2}");
    hMassCorrelations->GetZaxis()->SetMaxDigits(3);  // Forces scientific notation for large/small numbers
    hMassCorrelations->GetXaxis()->CenterTitle();
    hMassCorrelations->GetYaxis()->CenterTitle();
    hMassCorrelations->GetXaxis()->SetTitleOffset(2);    // Adjust X-axis title offset
    hMassCorrelations->GetYaxis()->SetTitleOffset(2.2);    // Adjust Y-axis title offset
    hMassCorrelations->GetZaxis()->SetTitleOffset(1.45); // Adjust Z-axis title offset
    hMassCorrelations->SetTitle("");

    TPaveText *ptTitle = new TPaveText(0.12, 0.8, 0.4, 0.94, "brNDC");
    ptTitle->SetTextFont(42);
    ptTitle->SetTextSize(0.04);
    ptTitle->SetBorderSize(0);
    ptTitle->SetFillStyle(0);
    ptTitle->SetTextAlign(11);

    // Add text lines
    ptTitle->AddText("#scale[1.3]{ALICE Work in Progress}");
    ptTitle->AddText("#scale[1.15]{pp, #sqrt{#it{s}} = 13.6 TeV}");

    TPaveText *pt = new TPaveText(0.5, 0.69, 0.85, 0.93, "brNDC");
    pt->SetTextFont(42);
    pt->SetTextSize(0.04);
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(11);

    // Add text lines
    pt->AddText("D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.");
    if (strcmp(_pairType, "OS") == 0)
    {
      pt->AddText("D^{0}#bar{D}^{0} + #bar{D}^{0}D^{0} pairs");
    }
    else if (strcmp(_pairType, "LS") == 0)
    {
      pt->AddText("D^{0}D^{0} + #bar{D}^{0}#bar{D}^{0} pairs");
    }
    pt->AddText(Form("%.1f < #it{p}^{D}_{T} < %.1f GeV/#it{c}", _ptMin, _ptMax));

    gStyle->SetOptStat(0);
    hMassCorrelations->Draw("LEGO2");
    histFit->Draw("sameLEGO");

    ptTitle->Draw("same");
    pt->Draw("same");
    TString c2dName_pdf = cName + ".pdf";
    c2d->SaveAs(c2dName_pdf);
    TString c2dName_eps = cName + ".eps";
    c2d->SaveAs(c2dName_eps);
    fout->cd();
    c2d->Write();

    TString c1dName = cName + "_colz";
    TCanvas *c1d = new TCanvas(c1dName, c1dName, 1000, 700);
    hMassCorrelations->Draw("colz");
    c1d->Write();
    c1d->Close();
  }
}

void InvMassFitter2D::plotProjectionsAfterFit(RooFitResult *fitResult, RooProdPdf *model, RooDataSet *dataset, TString saveName, TFile *fout, bool doReflections, const char* suffix)
{

  RooRealVar *massCand1 = _workspace.var("fMCand1");
  RooRealVar *massCand2 = _workspace.var("fMCand2");

  massCand1->setRange(_massMin, _massMax);
  massCand2->setRange(_massMin, _massMax);

  double binWidth = 0.005; // Desired bin width
  int nBins = (int)((_massMax - _massMin) / binWidth);

  // Set the calculated number of bins
  massCand1->setBins(nBins);
  massCand2->setBins(nBins);

  RooPlot *frameCand1 = massCand1->frame(RooFit::Range(_massMin, _massMax), RooFit::Title("Invariant-mass of Candidate 1")); // define the frame to plot
  RooPlot *frameCand2 = massCand2->frame(RooFit::Range(_massMin, _massMax), RooFit::Title("Invariant-mass of Candidate 2")); // define the frame to plot
  frameCand1->GetXaxis()->SetLabelSize(0.05);  // Increase font size for X-axis numbers
  frameCand1->GetYaxis()->SetLabelSize(0.05);  // Increase font size for Y-axis numbers
  
  frameCand2->GetXaxis()->SetLabelSize(0.05);
  frameCand2->GetYaxis()->SetLabelSize(0.05);
  dataset->plotOn(frameCand1, RooFit::Name("data_cand1"), RooFit::DrawOption("pz"));
  dataset->plotOn(frameCand2, RooFit::Name("data_cand2"), RooFit::DrawOption("pz"));

  RooDataHist rooDsBinned("rooDsBinned", "Binned version of dataset", RooArgSet(*massCand1, *massCand2), *dataset);
  // Compute chi^2
  RooAbsReal* chi2 = model->createChi2(rooDsBinned, RooFit::DataError(RooAbsData::SumW2), RooFit::Extended(true));

  int nPar = fitResult->floatParsFinal().getSize(); // get the number of free parameters
  // Compute number of bins
  int nBinsDs = rooDsBinned.numEntries();
  int ndof = nBinsDs - nPar;

  double chi2_ndf = chi2->getVal() / ndof;
  cout << "chi2: " << chi2->getVal() << ", ndof: " << ndof << ", chi2_ndf: " << chi2_ndf << endl;

  std::string histochiname = std::string("hChi2_ndf_") + suffix;
  TH1F* hChi2_ndf = new TH1F(histochiname.c_str(), "#chi^{2}/ndf", 1, 0, 1);
  hChi2_ndf->SetBinContent(1, chi2_ndf);
  fout->cd();
  hChi2_ndf->Write();

  TLegend *legend = new TLegend(0.6, 0.55, 0.94, 0.79);
  legend->SetBorderSize(0);
  TLegend *legend_2 = new TLegend(0.55, 0.59, 0.89, 0.79);
  legend_2->SetBorderSize(0);

  RooCurve* sgnReflFunc2DCurve   = nullptr;
  RooCurve* bkgReflFunc2DCurve   = nullptr;
  RooCurve* reflBkgFunc2DCurve   = nullptr;
  RooCurve* reflSgnFunc2DCurve   = nullptr;
  RooCurve* reflReflFunc2DCurve  = nullptr;
  if (doReflections)
  {
    model->plotOn(frameCand1, Components((std::string("bkgReflFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kViolet), Name("bkgReflFunc2D_curve"));
    bkgReflFunc2DCurve = (RooCurve *)frameCand1->findObject("bkgReflFunc2D_curve", RooCurve::Class());
    model->plotOn(frameCand1, Components((std::string("sgnReflFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kMagenta - 7), Name("sgnReflFunc2D_curve"));
    sgnReflFunc2DCurve = (RooCurve *)frameCand1->findObject("sgnReflFunc2D_curve", RooCurve::Class());
    model->plotOn(frameCand1, Components((std::string("reflBkgFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kCyan - 3), Name("reflBkgFunc2D_curve"));
    reflBkgFunc2DCurve = (RooCurve *)frameCand1->findObject("reflBkgFunc2D_curve", RooCurve::Class());
    model->plotOn(frameCand1, Components((std::string("reflSgnFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kMagenta - 5), Name("reflSgnFunc2D_curve"));
    reflSgnFunc2DCurve = (RooCurve *)frameCand1->findObject("reflSgnFunc2D_curve", RooCurve::Class());
    model->plotOn(frameCand1, Components((std::string("reflReflFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kAzure + 3), Name("reflReflFunc2D_curve"));
    reflReflFunc2DCurve = (RooCurve *)frameCand1->findObject("reflReflFunc2D_curve", RooCurve::Class());

    model->plotOn(frameCand2, Components((std::string("bkgReflFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kViolet), Name("bkgReflFunc2D_curve"));
    model->plotOn(frameCand2, Components((std::string("sgnReflFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kMagenta - 7), Name("sgnReflFunc2D_curve"));
    model->plotOn(frameCand2, Components((std::string("reflBkgFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kCyan - 3), Name("reflBkgFunc2D_curve"));
    model->plotOn(frameCand2, Components((std::string("reflSgnFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kMagenta - 5), Name("reflSgnFunc2D_curve"));
    model->plotOn(frameCand2, Components((std::string("reflReflFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(2), LineColor(kAzure + 3), Name("reflReflFunc2D_curve"));
  }

  // Plot the components of the model and keep references to the returned RooCurve objects
  model->plotOn(frameCand1, Components((std::string("bkgBkgFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(kDashed), LineColor(kRed), Name("bkgBkgFunc2D_curve"));
  RooCurve *bkgBkgFunc2DCurve = (RooCurve *)frameCand1->findObject("bkgBkgFunc2D_curve", RooCurve::Class());
  model->plotOn(frameCand1, Components((std::string("sgnBkgFunc2D") + suffix).c_str()), LineStyle(3), LineWidth(2), LineColor(kOrange + 6), Name("sgnBkgFunc2D_curve"));
  RooCurve *sgnBkgFunc2DCurve = (RooCurve *)frameCand1->findObject("sgnBkgFunc2D_curve", RooCurve::Class());
  model->plotOn(frameCand1, Components((std::string("bkgSgnFunc2D") + suffix).c_str()), LineStyle(4), LineWidth(2), LineColor(kOrange), Name("bkgSgnFunc2D_curve"));
  RooCurve *bkgSgnFunc2DCurve = (RooCurve *)frameCand1->findObject("bkgSgnFunc2D_curve", RooCurve::Class());
  model->plotOn(frameCand1, DrawOption("F"), Components((std::string("sgnSgnFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(1), LineColor(kGreen +2), FillColor(kGreen -2), Name("sgnSgnFunc2D_curve"));
  RooCurve *sgnSgnFunc2DCurve = (RooCurve *)frameCand1->findObject("sgnSgnFunc2D_curve", RooCurve::Class());

  // Now cand 2
  model->plotOn(frameCand2, Components((std::string("bkgBkgFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(kDashed), LineColor(kRed), Name("bkgBkgFunc2D_curve"));
  model->plotOn(frameCand2, Components((std::string("sgnBkgFunc2D") + suffix).c_str()), LineStyle(3), LineWidth(2), LineColor(kOrange + 6), Name("sgnBkgFunc2D_curve"));
  model->plotOn(frameCand2, Components((std::string("bkgSgnFunc2D") + suffix).c_str()), LineStyle(4), LineWidth(2), LineColor(kOrange), Name("bkgSgnFunc2D_curve"));
  model->plotOn(frameCand2, DrawOption("F"), Components((std::string("sgnSgnFunc2D") + suffix).c_str()), LineWidth(2), LineStyle(1), LineColor(kGreen + 2), FillColor(kGreen -2), Name("sgnSgnFunc2D_curve"));

  model->plotOn(frameCand1, Name("model_curve")); // Plot the full model
  RooCurve *modelMCCurve = (RooCurve *)frameCand1->findObject("model_curve", RooCurve::Class());
  model->plotOn(frameCand2, Name("model_curve")); // Plot the full model

  // Add entries to the legend
  legend->AddEntry("data_cand1", "Data", "zp");              // Entry for data
  legend->AddEntry("model_curve", "Total fit", "l");              // Entry for data

  if (strcmp(_pairType, "LS") == 0) {
    legend->AddEntry(sgnSgnFunc2DCurve, "Sig. D^{0} - Sig. D^{0}", "l"); // Entry for signal
    legend->AddEntry(bkgBkgFunc2DCurve, "Bkg. D^{0} - Bkg. D^{0}", "l"); // Entry for background
    legend->AddEntry(sgnBkgFunc2DCurve, "Sig. D^{0} - Bkg. D^{0}", "l"); // Entry for cross term sgn x bkg
    legend->AddEntry(bkgSgnFunc2DCurve, "Bkg. D^{0} - Sig. D^{0}", "l"); // Entry for cross term bkg x sgn
  } else if (strcmp(_pairType, "OS") == 0) {
    legend->AddEntry(sgnSgnFunc2DCurve, "Sig. D^{0} - Sig. #bar{D^{0}}", "l"); // Entry for signal
    legend->AddEntry(bkgBkgFunc2DCurve, "Bkg. D^{0} - Bkg. #bar{D^{0}}", "l"); // Entry for background
    legend->AddEntry(sgnBkgFunc2DCurve, "Sig. D^{0} - Bkg. #bar{D^{0}}", "l"); // Entry for cross term sgn x bkg
    legend->AddEntry(bkgSgnFunc2DCurve, "Bkg. D^{0} - Sig. #bar{D^{0}}", "l"); // Entry for cross term bkg x sgn
  }

  if (doReflections)
  {
    if (strcmp(_pairType, "LS") == 0) {
      legend_2->AddEntry(sgnReflFunc2DCurve, "Sig D^{0} - Refl D^{0}", "l");   // Entry for cross term sgn x refl
      legend_2->AddEntry(bkgReflFunc2DCurve, "Bkg D^{0} - Refl D^{0}", "l");   // Entry for cross term bkg x refl
      legend_2->AddEntry(reflSgnFunc2DCurve, "Refl D^{0} - Sig D^{0}", "l");   // Entry for cross term refl x sgn
      legend_2->AddEntry(reflBkgFunc2DCurve, "Refl D^{0} - Bkg D^{0}", "l");   // Entry for cross term refl x bkg
      legend_2->AddEntry(reflReflFunc2DCurve, "Refl D^{0} - Refl D^{0}", "l"); // Entry for cross term refl x refl
    } else if (strcmp(_pairType, "OS") == 0) {
      legend_2->AddEntry(sgnReflFunc2DCurve, "Sig D^{0} - Refl #bar{D^{0}}", "l");   // Entry for cross term sgn x refl
      legend_2->AddEntry(bkgReflFunc2DCurve, "Bkg D^{0} - Refl #bar{D^{0}}", "l");   // Entry for cross term bkg x refl
      legend_2->AddEntry(reflSgnFunc2DCurve, "Refl D^{0} - Sig #bar{D^{0}}", "l");   // Entry for cross term refl x sgn
      legend_2->AddEntry(reflBkgFunc2DCurve, "Refl D^{0} - Bkg #bar{D^{0}}", "l");   // Entry for cross term refl x bkg
      legend_2->AddEntry(reflReflFunc2DCurve, "Refl D^{0} - Refl #bar{D^{0}}", "l"); // Entry for cross term refl x refl
    }
  }

  // Display chi2/ndf and parameter values
  TPaveText *chi2Text = new TPaveText(0.17, 0.6, 0.4, 0.67, "brNDC");
  chi2Text->SetTextFont(42);
  chi2Text->SetTextSize(0.04);
  chi2Text->SetBorderSize(0);
  chi2Text->SetFillStyle(0);
  chi2Text->SetTextAlign(11);

  // Add text lines
  if (strcmp(suffix, "") == 0 && chi2_ndf != 0) {
    chi2Text->AddText(Form("#chi^{2}/ndf = %.3f", chi2_ndf));
  }

  TPaveText *ptTitle = new TPaveText(0.17, 0.8, 0.45, 0.9, "brNDC");
  ptTitle->SetTextFont(42);
  ptTitle->SetTextSize(0.04);
  ptTitle->SetBorderSize(0);
  ptTitle->SetFillStyle(0);
  ptTitle->SetTextAlign(11);

    // Add text lines
  ptTitle->AddText("#scale[1.35]{ALICE Work in Progress}");
  ptTitle->AddText("#scale[1.2]{pp, #sqrt{#it{s}} = 13.6 TeV}");

  TPaveText *pt = new TPaveText(0.17, 0.69, 0.5, 0.8, "brNDC");
    pt->SetTextFont(42);
    pt->SetTextSize(0.04);
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetTextAlign(11);

    // Add text lines
    pt->AddText(Form("#scale[1.1]{%.1f < #it{p}_{T} < %.1f GeV/#it{c}}", _ptMin, _ptMax));
    pt->AddText("#scale[1.1]{|#it{y}| < 0.5}");

    TPaveText *pt_cand1 = new TPaveText(0.17, 0.75, 0.7, 0.9, "brNDC");
    pt_cand1->SetTextFont(42);
    pt_cand1->SetTextSize(0.04);
    pt_cand1->SetBorderSize(0);
    pt_cand1->SetFillStyle(0);
    pt_cand1->SetTextAlign(11);

    // Add text lines
    pt_cand1->AddText("#scale[1.1]{D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.}");
    if (strcmp(_pairType, "OS") == 0)
    {
      pt_cand1->AddText("#scale[1.1]{D^{0}#bar{D}^{0} + #bar{D}^{0}D^{0} pairs}");
    }
    else if (strcmp(_pairType, "LS") == 0)
    {
      pt_cand1->AddText("#scale[1.1]{D^{0}D^{0} + #bar{D}^{0}#bar{D}^{0} pairs}");
    }

  // Create a canvas to draw the projections
  TCanvas *c_proj = new TCanvas(saveName, saveName, 1250, 600);
  c_proj->Divide(2, 1);     // Divide canvas into two pads
  gPad->SetLeftMargin(0.2); // Increase the left margin (default is ~0.1)


  // Draw the X projection
  c_proj->cd(1); // Go to the first pad
  gPad->SetTicks(1, 1);  // Enable ticks
  gPad->SetLeftMargin(0.15);   // Adjust left margin (e.g., 15%)
  gPad->SetRightMargin(0.02);  // Adjust right margin
  gPad->SetTopMargin(0.05);    // Adjust top margin
  gPad->SetBottomMargin(0.15); // Adjust bottom margin

  double maxY1 = frameCand1->GetMaximum();
  frameCand1->SetMaximum(maxY1 + maxY1/4);

  frameCand1->SetTitle("");
  frameCand1->SetXTitle("#it{M}_{K#pi} (GeV/#it{c}^{2})");
  frameCand1->SetYTitle("Counts per 5 MeV/#it{c}^{2}");
  frameCand1->GetXaxis()->SetTitleSize(0.05); // Set X-axis title size
  frameCand1->GetYaxis()->SetTitleSize(0.05); // Set Y-axis title size

  frameCand1->Draw();
  ptTitle->Draw("same");
  legend->Draw();
  
  pt->Draw("same");

  // Draw the Y projection
  c_proj->cd(2); // Go to the second pad
  gPad->SetTicks(1, 1);  // Enable ticks
  gPad->SetLeftMargin(0.15);   // Adjust left margin (e.g., 15%)
  gPad->SetRightMargin(0.02);  // Adjust right margin
  gPad->SetTopMargin(0.05);    // Adjust top margin
  gPad->SetBottomMargin(0.15); // Adjust bottom margin

  frameCand2->SetTitle("");
  frameCand2->SetYTitle("Counts per 5 MeV/#it{c}^{2}");
  frameCand2->SetXTitle("#it{M}_{K#pi} (GeV/#it{c}^{2})");
  frameCand2->GetXaxis()->SetTitleSize(0.05); // Set X-axis title size
  frameCand2->GetYaxis()->SetTitleSize(0.05); // Set Y-axis title size

  double maxY2 = frameCand2->GetMaximum();
  frameCand2->SetMaximum(maxY2 + maxY2/4);

  frameCand2->Draw();
  chi2Text->Draw("same");
  pt_cand1->Draw("same");

  if (doReflections)
  {
    legend_2->Draw();
  }

  // Show the canvas
  c_proj->Draw();
  fout->cd();
  c_proj->Write();
  TString name = saveName + ".pdf";
  c_proj->SaveAs(name);

  // Define TH1F histograms with 50 bins in the range (_massMin, _massMax)
  TH1F *histoCand1 = new TH1F("histoCand1", "Invariant Mass of Candidate 1", 10, _massMin, _massMax);
  TH1F *histoCand2 = new TH1F("histoCand2", "Invariant Mass of Candidate 2", 10, _massMin, _massMax);

  // Fill the histograms with data from the RooPlot frames
  for (int i = 1; i <= histoCand1->GetNbinsX(); ++i)
  {
    double x = histoCand1->GetBinCenter(i);
    double y1 = ((RooCurve*)frameCand1->getHist("data_cand1"))->Eval(x);
    double y2 = ((RooCurve*)frameCand2->getHist("data_cand2"))->Eval(x);

    histoCand1->SetBinContent(i, y1);
    histoCand2->SetBinContent(i, y2);

    // Optionally set errors, assuming you have uncertainties
    histoCand1->SetBinError(i, std::sqrt(y1)); // for Poisson stats
    histoCand2->SetBinError(i, std::sqrt(y2));
  }

  histoCand1->Divide(histoCand2);
  TLine *line = new TLine(_massMin, 1, _massMax, 1);

  line->SetLineColor(kBlack);
  line->SetLineWidth(2);

  TCanvas *c_ratio = new TCanvas("c_ratioCands", "c_ratioCandsMasses", 800, 600);
  histoCand1->SetTitle("massCand1 / massCand2");
  histoCand1->GetXaxis()->SetTitle("Inv. mass");
  histoCand1->GetYaxis()->SetTitle("Ratio");
  histoCand1->SetMarkerStyle(20);
  histoCand1->SetMarkerColor(kBlue + 2);
  histoCand1->Draw("PE");
  line->Draw("same");
  c_ratio->Write();
  c_ratio->SaveAs("ratio_massCand1OverMassCand2.pdf");
  c_ratio->Close();

  // Plot residuals
  RooPlot *residFrameCand1 = massCand1->frame(RooFit::Title("Residuals of Candidate 1"));
  RooPlot *residFrameCand2 = massCand2->frame(RooFit::Title("Residuals of Candidate 2"));

  // Compute residuals (Data - Fit) and add to the frames
  RooHist* residualsCand1 = (RooHist*)frameCand1->residHist("data_cand1", "model_curve");
  RooHist* residualsCand2 = (RooHist*)frameCand2->residHist("data_cand2", "model_curve");

  // Add residuals to the frames
  residualsCand1->SetMarkerColor(kBlue);
  residualsCand2->SetMarkerColor(kBlue);
  residualsCand1->SetLineColor(kBlue);
  residualsCand2->SetLineColor(kBlue);
  residFrameCand1->addPlotable(residualsCand1, "P");
  residFrameCand2->addPlotable(residualsCand2, "P");

  TCanvas *cResiduals = new TCanvas("cResiduals", "Residuals", 800, 600);
  cResiduals->Divide(1, 2); // Divide into two pads for cand1 and cand2

  cResiduals->cd(1);
  residFrameCand1->SetYTitle("Counts per 5 MeV/#it{c}^{2}");
  residFrameCand1->SetXTitle("M(K#pi) (GeV/#it{c}^{2})");
  residFrameCand1->GetXaxis()->SetTitleSize(0.05); // Set X-axis title size
  residFrameCand1->GetYaxis()->SetTitleSize(0.05); // Set Y-axis title size
  residFrameCand1->Draw();
  line->Draw("same");

  cResiduals->cd(2);
  residFrameCand2->SetYTitle("Counts per 5 MeV/#it{c}^{2}");
  residFrameCand2->SetXTitle("M(K#pi) (GeV/#it{c}^{2})");
  residFrameCand2->GetXaxis()->SetTitleSize(0.05); // Set X-axis title size
  residFrameCand2->GetYaxis()->SetTitleSize(0.05); // Set Y-axis title size
  residFrameCand2->Draw();
  line->Draw("same");

  TString nameRes = "residuals" + saveName + ".pdf";
  cResiduals->SaveAs(nameRes); // Save to file if needed
  cResiduals->Write();
  cResiduals->Close();

  delete c_ratio;
  delete cResiduals;
  delete histoCand1;
  delete histoCand2;

}

double InvMassFitter2D::calculateWeights(double const &y, double const &pt)
{
  int nbinsPt = _efficiencyMap->GetXaxis()->GetNbins();
  int nbinsY = _efficiencyMap->GetYaxis()->GetNbins();

  int lowerBinPt = _efficiencyMap->GetXaxis()->FindBin(_ptMin);
  int upperBinPt = _efficiencyMap->GetXaxis()->FindBin(_ptMax);

  double weight = 1.;
  const double epsilon = 1e-4; // Small value to avoid division by zero

  for (int ipt = lowerBinPt; ipt <= upperBinPt; ipt++)
  {
    double ptMin = _efficiencyMap->GetXaxis()->GetBinLowEdge(ipt);
    double ptMax = _efficiencyMap->GetXaxis()->GetBinUpEdge(ipt);

    for (int iy = 1; iy <= nbinsY; iy++)
    {
      double yMin = _efficiencyMap->GetYaxis()->GetBinLowEdge(iy);
      double yMax = _efficiencyMap->GetYaxis()->GetBinUpEdge(iy);

      // Include ptMax and yMax for the last bin
      bool isLastPtBin = (ipt == nbinsPt);
      bool isLastYBin = (iy == nbinsY);

      if (pt >= ptMin && (pt < ptMax || (isLastPtBin && pt <= ptMax)) &&
          y >= yMin && (y < yMax || (isLastYBin && y <= yMax)))
      {
        double efficiency = _efficiencyMap->GetBinContent(ipt, iy);
        if (efficiency > epsilon)
        {
          weight = 1. / efficiency;
          return weight;
        }
      }
    }
  }
  return weight;
}

// Create Lorentz vector to calculate the pt of the pair
ROOT::Math::PxPyPzMVector InvMassFitter2D::createLorentzVector(double const &phi, double const &y, double const &pt, double const &m)
{
  // Calculate momentum components
  double const px = pt * cos(phi);
  double const py = pt * sin(phi);
  double const pz = pt * sinh(y);

  // Create the Lorentz vector
  ROOT::Math::PxPyPzMVector lorentzVector(px, py, pz, m);

  return lorentzVector;
}

RooRealVar *InvMassFitter2D::getYieldInRange(RooFitResult *fitResult, RooRealVar *massCand1, RooRealVar *massCand2, RooProdPdf function, RooFormulaVar nCands, TString range)
{
  // Calculate the integrals over the specified range and the full range
  RooAbsReal *intRange = function.createIntegral(RooArgSet(*massCand1, *massCand2), Range(range));
  RooAbsReal *intTotal = function.createIntegral(RooArgSet(*massCand1, *massCand2));

  double yieldRange = nCands.getVal() * intRange->getVal() / intTotal->getVal();

  // Propagate the uncertainties
  double nCandsErr = nCands.getPropagatedError(*fitResult);

  // Calculate the uncertainty using error propagation
  double errorRange = nCandsErr * intRange->getVal() / intTotal->getVal();

  cout << "Integral: " << intRange->getVal() / intTotal->getVal() << endl;

  // Initialize the RooRealVar to hold the yield and error
  RooRealVar *yield = new RooRealVar("yield", "Yield in 3-sigma range", yieldRange);
  yield->setError(errorRange);
  return yield;
}

void InvMassFitter2D::plotFitResults(RooDataSet* dataset, RooRealVar* mass, RooAbsPdf* model, RooFitResult* fitResult,
                                     const char* sgnComponent, const char* bkgComponent, const char* title, const char* canvasName) {
  // Create a RooPlot for the mass variable
  RooPlot* frame = mass->frame(Title(title));

  // Get chi^2/ndf
  // Only works with binned data!!!
  // Convert unbinned dataset into binned dataset
  RooDataHist rooDsBinned("rooDsBinned", "Binned version of dataset", RooArgSet(*mass), *dataset);
  // Compute chi^2
  RooAbsReal* chi2 = model->createChi2(rooDsBinned, RooFit::DataError(RooAbsData::SumW2), RooFit::Extended(true));

  int nPar = fitResult->floatParsFinal().getSize(); // get the number of free parameters
  // Compute number of bins
  int nBins = rooDsBinned.numEntries();
  int ndof = nBins - nPar;
  
  double chi2_ndf = chi2->getVal() / ndof;
  cout << chi2_ndf << endl;

  // Plot the data on the frame
  dataset->plotOn(frame, Name("data"));

  // Plot the components of the model (signal and background)
  model->plotOn(frame, Components(sgnComponent), LineStyle(kDashed), LineColor(kRed), Name("signal"));
  model->plotOn(frame, Components(bkgComponent), LineStyle(kDashed), LineColor(kGreen), Name("background"));

  // Plot the fitted model on the frame
  model->plotOn(frame, Name("model"), LineColor(kBlue));

  // Customize the plot
  frame->SetTitle(title);
  frame->GetXaxis()->SetTitle("Mass");
  frame->GetYaxis()->SetTitle("Events");

  // Create a legend
  TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend->AddEntry(frame->findObject("data"), "Data", "lep");
  legend->AddEntry(frame->findObject("model"), "Fit Model", "l");
  legend->AddEntry(frame->findObject("signal"), "Signal Component", "l");
  legend->AddEntry(frame->findObject("background"), "Background Component", "l");

  // Display chi2/ndf and parameter values
  TPaveText *chi2Text = new TPaveText(0.12, 0.8, 0.4, 0.85, "brNDC");
  chi2Text->SetTextFont(42);
  chi2Text->SetTextSize(0.04);
  chi2Text->SetBorderSize(0);
  chi2Text->SetFillStyle(0);
  chi2Text->SetTextAlign(11);

  // Add text lines
  chi2Text->AddText(Form("#chi^{2}/ndf = %.3f", chi2_ndf));

  // Draw the plot on a TCanvas
  TCanvas* canvas = new TCanvas(canvasName, title, 800, 600);
  frame->Draw();
  chi2Text->Draw("same");
  legend->Draw();
  canvas->Update();

  // Optionally, save the plot to a file
  canvas->SaveAs(Form("%s.png", canvasName));
  delete canvas;
}

void InvMassFitter2D::plotFitResults(RooDataSet* dataset, RooRealVar* mass, RooAbsPdf* model, RooFitResult* fitResult,
                                     const char* sgnComponent, const char* bkgComponent, const char* reflComponent, const char* title, const char* canvasName) {
  // Create a RooPlot for the mass variable
  RooPlot* frame = mass->frame(Title(title));

  // Get chi^2/ndf
  // Only works with binned data!!!
  // Convert unbinned dataset into binned dataset
  RooDataHist rooDsBinned("rooDsBinned", "Binned version of dataset", RooArgSet(*mass), *dataset);
  // Compute chi^2
  RooAbsReal* chi2 = model->createChi2(rooDsBinned, RooFit::DataError(RooAbsData::SumW2), RooFit::Extended(true));

  int nPar = fitResult->floatParsFinal().getSize(); // get the number of free parameters
  // Compute number of bins
  int nBins = rooDsBinned.numEntries();
  int ndof = nBins - nPar;
  
  double chi2_ndf = chi2->getVal() / ndof;
  cout << chi2_ndf << endl;

  // Plot the data on the frame
  dataset->plotOn(frame, Name("data"));

  // Plot the components of the model (signal and background)
  model->plotOn(frame, Components(sgnComponent), LineStyle(kDashed), LineColor(kGreen), Name("signal"));
  model->plotOn(frame, Components(bkgComponent), LineStyle(kDashed), LineColor(kRed), Name("background"));
  model->plotOn(frame, Components(reflComponent), LineStyle(kDashed), LineColor(kPink), Name("reflections"));

  // Plot the fitted model on the frame
  model->plotOn(frame, Name("model"), LineColor(kBlue));

  // Customize the plot
  frame->SetTitle(title);
  frame->GetXaxis()->SetTitle("Mass");
  frame->GetYaxis()->SetTitle("Events");

  // Create a legend
  TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend->AddEntry(frame->findObject("data"), "Data", "lep");
  legend->AddEntry(frame->findObject("model"), "Fit Model", "l");
  legend->AddEntry(frame->findObject("signal"), "Signal Component", "l");
  legend->AddEntry(frame->findObject("background"), "Background Component", "l");
  legend->AddEntry(frame->findObject("reflections"), "Reflected Component", "l");

  // Display chi2/ndf and parameter values
  TPaveText *chi2Text = new TPaveText(0.12, 0.8, 0.4, 0.85, "brNDC");
  chi2Text->SetTextFont(42);
  chi2Text->SetTextSize(0.04);
  chi2Text->SetBorderSize(0);
  chi2Text->SetFillStyle(0);
  chi2Text->SetTextAlign(11);

  // Add text lines
  chi2Text->AddText(Form("#chi^{2}/ndf = %.3f", chi2_ndf));

  // Draw the plot on a TCanvas
  TCanvas* canvas = new TCanvas(canvasName, title, 800, 600);
  frame->Draw();
  chi2Text->Draw("same");
  legend->Draw();
  canvas->Update();

  // Optionally, save the plot to a file
  canvas->SaveAs(Form("%s.png", canvasName));
  delete canvas;
}


void InvMassFitter2D::selectFitFunctions(RooAbsPdf* &sgnPdfCand1,RooAbsPdf* &sgnPdfCand2,RooAbsPdf* &bkgPdfCand1,
                                         RooAbsPdf* &bkgPdfCand2, RooAbsPdf* &reflPdfCand1, RooAbsPdf* &reflPdfCand2, int8_t fitType) {
  // Select signal function
  if (_sgnFuncOption == "gaus")
  {
    cout << "Signal function chosen: GAUSSIAN" << endl;
    if (fitType == Prefit) {
      sgnPdfCand1 = _workspace.pdf("sgnFuncGausCand1Prefit");
      sgnPdfCand2 = _workspace.pdf("sgnFuncGausCand2Prefit");
    } else if (fitType == Rawfit) {
      sgnPdfCand1 = _workspace.pdf("sgnFuncGausCand1");
      sgnPdfCand2 = _workspace.pdf("sgnFuncGausCand2");
    } else if (fitType == Corrfit) {
      sgnPdfCand1 = _workspace.pdf("sgnFuncGausCand1Corr");
      sgnPdfCand2 = _workspace.pdf("sgnFuncGausCand2Corr");
    } else if (fitType == PrefitCorr) {
      sgnPdfCand1 = _workspace.pdf("sgnFuncGausCand1CorrPrefit");
      sgnPdfCand2 = _workspace.pdf("sgnFuncGausCand2CorrPrefit");
    }
  }
  else if (_sgnFuncOption == "CB")
  {
    cout << "Signal function chosen: BIFURCATED CRYSTAL BALL" << endl;
    if (fitType == Prefit) {
      sgnPdfCand1 = _workspace.pdf("sgnFuncCBCand1Prefit");
      sgnPdfCand2 = _workspace.pdf("sgnFuncCBCand2Prefit");
    } else if (fitType == Rawfit) {
      sgnPdfCand1 = _workspace.pdf("sgnFuncCBCand1");
      sgnPdfCand2 = _workspace.pdf("sgnFuncCBCand2");
    } else if (fitType == Corrfit) {
      sgnPdfCand1 = _workspace.pdf("sgnFuncCBCand1Corr");
      sgnPdfCand2 = _workspace.pdf("sgnFuncCBCand2Corr");
    } else if (fitType == PrefitCorr) {
      sgnPdfCand1 = _workspace.pdf("sgnFuncCBCand1CorrPrefit");
      sgnPdfCand2 = _workspace.pdf("sgnFuncCBCand2CorrPrefit");
    }
  } else {
    cerr << "ERROR: signal function not supported! \n Available options: gaus, CB. \n Exit!" << endl;
    return;
  }

  if (!sgnPdfCand1 || !sgnPdfCand2)
  {
    cerr << "ERROR: sgnPdf not found!" << endl;
    return;
  }

  // Select background function
  if (_bkgFuncOption == "expo")
  {
    cout << "Background function chosen: EXPONENTIAL" << endl;
    if (fitType == Prefit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpoCand1Prefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpoCand2Prefit");
    } else if (fitType == Rawfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpoCand1");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpoCand2");
    } else if (fitType == Corrfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpoCand1Corr");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpoCand2Corr");
    } else if (fitType == PrefitCorr) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpoCand1CorrPrefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpoCand2CorrPrefit");
    }
  }
  else if (_bkgFuncOption == "poly0")
  {
    cout << "Background function chosen: POLY 0" << endl;
    if (fitType == Prefit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly0Cand1Prefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly0Cand2Prefit");
    } else if (fitType == Rawfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly0Cand1");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly0Cand2");
    } else if (fitType == Corrfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly0Cand1Corr");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly0Cand2Corr");
    } else if (fitType == PrefitCorr) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly0Cand1CorrPrefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly0Cand2CorrPrefit");
    }
  }
  else if (_bkgFuncOption == "poly1")
  {
    cout << "Background function chosen: POLY 1" << endl;
    if (fitType == Prefit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly1Cand1Prefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly1Cand2Prefit");
    } else if (fitType == Rawfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly1Cand1");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly1Cand2");
    } else if (fitType == Corrfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly1Cand1Corr");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly1Cand2Corr");
    } else if (fitType == PrefitCorr) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly1Cand1CorrPrefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly1Cand2CorrPrefit");
    }
  }
  else if (_bkgFuncOption == "poly2")
  {
    cout << "Background function chosen: POLY 2" << endl;
    if (fitType == Prefit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly2Cand1Prefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly2Cand2Prefit");
    } else if (fitType == Rawfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly2Cand1");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly2Cand2");
    } else if (fitType == Corrfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly2Cand1Corr");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly2Cand2Corr");
    } else if (fitType == PrefitCorr) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncPoly2Cand1CorrPrefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncPoly2Cand2CorrPrefit");
    }
  }
  else if (_bkgFuncOption == "expPoly2") {
    cout << "Background function chosen: EXPO * POLY 2" << endl;
    if (fitType == Prefit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpPoly2Cand1Prefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpPoly2Cand2Prefit");
    } else if (fitType == Rawfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpPoly2Cand1");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpPoly2Cand2");
    } else if (fitType == Corrfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpPoly2Cand1Corr");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpPoly2Cand2Corr");
    } else if (fitType == PrefitCorr) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpPoly2Cand1CorrPrefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpPoly2Cand2CorrPrefit");
    }
  }
  else if (_bkgFuncOption == "expPoly1") {
    cout << "Background function chosen: EXPO * POLY 1" << endl;
    if (fitType == Prefit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpPoly1Cand1Prefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpPoly1Cand2Prefit");
    } else if (fitType == Rawfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpPoly1Cand1");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpPoly1Cand2");
    } else if (fitType == Corrfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpPoly1Cand1Corr");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpPoly1Cand2Corr");
    } else if (fitType == PrefitCorr) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExpPoly1Cand1CorrPrefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExpPoly1Cand2CorrPrefit");
    }
  }
  else if (_bkgFuncOption == "exp2") {
    cout << "Background function chosen: EXPO^(X^2)" << endl;
    if (fitType == Prefit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExp2Cand1Prefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExp2Cand2Prefit");
    } else if (fitType == Rawfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExp2Cand1");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExp2Cand2");
    } else if (fitType == Corrfit) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExp2Cand1Corr");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExp2Cand2Corr");
    } else if (fitType == PrefitCorr) {
      bkgPdfCand1 = _workspace.pdf("bkgFuncExp2Cand1CorrPrefit");
      bkgPdfCand2 = _workspace.pdf("bkgFuncExp2Cand2CorrPrefit");
    }
  }

  else
  {
    cerr << "ERROR: background function not supported! \n Available options: expo, poly2, poly3, expPoly3, exp2. \n Exit!" << endl;
    return;
  }

  if (!bkgPdfCand1 || !bkgPdfCand2)
  {
    cerr << "ERROR: bkgPdf not found!" << endl;
    return;
  }

  // Select reflection function
  if (_reflFuncOption == "gaus")
  {
    cout << "Reflected function chosen: GAUSSIAN" << endl;
    if (fitType == Prefit) {
      reflPdfCand1 = _workspace.pdf("reflFuncGausCand1Prefit");
      reflPdfCand2 = _workspace.pdf("reflFuncGausCand2Prefit");
    } else if (fitType == Rawfit) {
      reflPdfCand1 = _workspace.pdf("reflFuncGausCand1");
      reflPdfCand2 = _workspace.pdf("reflFuncGausCand2");
    } else if (fitType == Corrfit) {
      reflPdfCand1 = _workspace.pdf("reflFuncGausCand1Corr");
      reflPdfCand2 = _workspace.pdf("reflFuncGausCand2Corr");
    } else if (fitType == PrefitCorr) {
      reflPdfCand1 = _workspace.pdf("reflFuncGausCand1CorrPrefit");
      reflPdfCand2 = _workspace.pdf("reflFuncGausCand2CorrPrefit");
    }
  }
  else if (_reflFuncOption == "doubleGaus")
  {
    cout << "Reflected function chosen: DOUBLE GAUSSIAN" << endl;
    if (fitType == Prefit) {
      reflPdfCand1 = _workspace.pdf("reflFuncDoubleGausCand1Prefit");
      reflPdfCand2 = _workspace.pdf("reflFuncDoubleGausCand2Prefit");
    } else if (fitType == Rawfit) {
      reflPdfCand1 = _workspace.pdf("reflFuncDoubleGausCand1");
      reflPdfCand2 = _workspace.pdf("reflFuncDoubleGausCand2");
    } else if (fitType == Corrfit) {
      reflPdfCand1 = _workspace.pdf("reflFuncDoubleGausCand1Corr");
      reflPdfCand2 = _workspace.pdf("reflFuncDoubleGausCand2Corr");
    } else if (fitType == PrefitCorr) {
      reflPdfCand1 = _workspace.pdf("reflFuncDoubleGausCand1CorrPrefit");
      reflPdfCand2 = _workspace.pdf("reflFuncDoubleGausCand2CorrPrefit");
    }
  }
  else
  {
    cerr << "ERROR: refl function not supported! \n Available options: gaus, doubleGaus. \n Exit!" << endl;
    return;
  }

  if (!reflPdfCand1 || !reflPdfCand2)
  {
    cerr << "ERROR: reflPdf function not found!" << endl;
    return;
  }
}

RooFitResult *InvMassFitter2D::fitAndPlot1DCandidate(RooAbsPdf* sgnPdf, RooAbsPdf* bkgPdf, RooRealVar& massVar, RooDataSet* dataset,
                                            const std::string& candidateName, const std::string& plotFilename) {
  // Define the number of signal and background events

  RooRealVar nSgn(Form("nSgn%s", candidateName.c_str()),
  Form("Number of signal events of %s", candidateName.c_str()),
        100000, 0, 1e8);
  RooRealVar nBkg(Form("nBkg%s", candidateName.c_str()),
  Form("Number of background events of %s", candidateName.c_str()),
    500000, 0, 1e8);

  RooRealVar nSgnCorr(Form("nSgnCorr%s", candidateName.c_str()),
    Form("Number of signal events of %s", candidateName.c_str()),
        1e8, 0, 1e12);
  RooRealVar nBkgCorr(Form("nBkgCorr%s", candidateName.c_str()),
    Form("Number of background events of %s", candidateName.c_str()),
        1e8, 0, 1e12);

  // Create the composite model
  RooAddPdf model(Form("model%s", candidateName.c_str()),
  Form("Signal + Background of %s", candidateName.c_str()),
    RooArgList(*sgnPdf, *bkgPdf), RooArgList(nSgn, nBkg));
  
  RooAddPdf modelCorr(Form("modelCorr%s", candidateName.c_str()),
  Form("Signal + Background of %s", candidateName.c_str()),
    RooArgList(*sgnPdf, *bkgPdf), RooArgList(nSgnCorr, nBkgCorr));

  // Fit the model to the data
  RooFitResult* fitResult;
  if (plotFilename.find("Corr") != std::string::npos) {
    // The string contains "Corr"
    fitResult =  modelCorr.fitTo(*dataset, Save(), RooFit::PrintLevel(-1));
  } else {
    fitResult =  model.fitTo(*dataset, Save(), RooFit::PrintLevel(-1));
  }

  fitResult->Print("v");
  if (fitResult->status() != 0) {
    std::cout << "Fit did not converge for " << candidateName << "!" << std::endl;
  }

  // Plot the fit results
  std::string title = "Fit Result for " + candidateName;
  const char* sgnName = sgnPdf->GetName();
  const char* bkgName = bkgPdf->GetName();

  if (plotFilename.find("Corr") != std::string::npos) {
    plotFitResults(dataset, &massVar, &modelCorr, fitResult, sgnName, bkgName, title.c_str(), plotFilename.c_str());
  } else {
    plotFitResults(dataset, &massVar, &model, fitResult, sgnName, bkgName, title.c_str(), plotFilename.c_str());
  }

  return fitResult;
}

RooFitResult *InvMassFitter2D::fitAndPlot1DCandidate(RooAbsPdf* sgnPdf, RooAbsPdf* bkgPdf, RooAbsPdf* reflPdf, RooRealVar& massVar, RooDataSet* dataset,
                                            const std::string& candidateName, const std::string& plotFilename) {
  // Define the number of signal and background events

  RooRealVar nSgn(Form("nSgn%s", candidateName.c_str()),
  Form("Number of signal events of %s", candidateName.c_str()),
        100000, 0, 1e8);
  RooFormulaVar *nRefl = new RooFormulaVar("nRefl", Form("reflected signal yield of %s", candidateName.c_str()), "@0 * @1", RooArgList(_reflOverSgn, nSgn));
  RooRealVar nBkg(Form("nBkg%s", candidateName.c_str()),
  Form("Number of background events of %s", candidateName.c_str()),
    500000, 0, 1e8);

  RooRealVar nSgnCorr(Form("nSgnCorr%s", candidateName.c_str()),
    Form("Number of signal events of %s", candidateName.c_str()),
        1e8, 0, 1e12);
  RooRealVar nBkgCorr(Form("nBkgCorr%s", candidateName.c_str()),
    Form("Number of background events of %s", candidateName.c_str()),
        1e8, 0, 1e12);
  RooFormulaVar *nReflCorr = new RooFormulaVar("nReflCorr", Form("reflected signal yield of %s", candidateName.c_str()), "@0 * @1", RooArgList(_reflOverSgn, nSgnCorr));

  // Create the composite model
  RooAddPdf model(Form("model%s", candidateName.c_str()),
  Form("Signal + Background + Reflections of %s", candidateName.c_str()),
    RooArgList(*sgnPdf, *bkgPdf, *reflPdf), RooArgList(nSgn, nBkg, *nRefl));
  
  RooAddPdf modelCorr(Form("modelCorr%s", candidateName.c_str()),
  Form("Signal + Background + Reflections of %s", candidateName.c_str()),
    RooArgList(*sgnPdf, *bkgPdf, *reflPdf), RooArgList(nSgnCorr, nBkgCorr, *nReflCorr));

  // Fit the model to the data
  RooFitResult* fitResult;
  if (plotFilename.find("Corr") != std::string::npos) {
    // The string contains "Corr"
    fitResult =  modelCorr.fitTo(*dataset, Save(), RooFit::PrintLevel(-1));
  } else {
    fitResult =  model.fitTo(*dataset, Save(), RooFit::PrintLevel(-1));
  }

  fitResult->Print("v");
  if (fitResult->status() != 0) {
    std::cout << "Fit did not converge for " << candidateName << "!" << std::endl;
  }

  // Plot the fit results
  std::string title = "Fit Result for " + candidateName;
  const char* sgnName = sgnPdf->GetName();
  const char* bkgName = bkgPdf->GetName();
  const char* reflName = reflPdf->GetName();

  if (plotFilename.find("Corr") != std::string::npos) {
    plotFitResults(dataset, &massVar, &modelCorr, fitResult, sgnName, bkgName, reflName, title.c_str(), plotFilename.c_str());
  } else {
    plotFitResults(dataset, &massVar, &model, fitResult, sgnName, bkgName, reflName, title.c_str(), plotFilename.c_str());
  }

  return fitResult;
}

void InvMassFitter2D::setPrefitParameters(const RooArgList& prefitParams) {
  RooArgSet workspaceVars = _workspace.allVars();
  // Loop through all parameters and print their values
  for (int i = 0; i < prefitParams.getSize(); i++) {
    RooRealVar* par = dynamic_cast<RooRealVar*>(prefitParams.at(i));
    if (!par) continue;  // Safety check

    std::string parName = par->GetName();
    std::string baseName = parName;

    // Remove "Prefit" suffix if it exists
    std::string suffix = "Prefit";
    if (parName.size() > suffix.size() &&
        parName.compare(parName.size() - suffix.size(), suffix.size(), suffix) == 0) {
      baseName = parName.substr(0, parName.size() - suffix.size());
    }

    std::cout << baseName << " = " << par->getVal() 
              << " ± " << par->getError() << std::endl;

    // Check if the parameter exists in the workspace
    RooRealVar* wsPar = dynamic_cast<RooRealVar*>(workspaceVars.find(baseName.c_str()));
    if (wsPar) {
        // Update the workspace variable with the fitted value
        wsPar->setVal(par->getVal());
        std::cout << "Updated workspace parameter " << parName 
                  << " to " << wsPar->getVal() << std::endl;
        //if (baseName.find("Corr") == std::string::npos) {
            if (baseName != "sigma" && baseName != "mean" && baseName != "sigmaCand2" && baseName != "meanCand2" &&
              baseName != "sigmaCorr" && baseName != "meanCorr" && baseName != "sigmaCand2Corr" && baseName != "meanCand2Corr") {
                wsPar->setConstant(kTRUE);
            }
        //}
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////
///////////////// BACKGROUND SUBTRACTION OF SIGNAL REGION /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

/*   // Create a flag variable
  RooCategory signalFlagCand1("signalFlagCand1", "Signal or Background Flag");
  signalFlagCand1.defineType("Signal", 1);
  signalFlagCand1.defineType("Background", 0);

  RooCategory signalFlagCand2("signalFlagCand2", "Signal or Background Flag");
  signalFlagCand2.defineType("Signal", 1);
  signalFlagCand2.defineType("Background", 0);

  // Create a new dataset with the signal flag
  RooArgSet varsWithFlag(*dataset->get());
  varsWithFlag.add(signalFlagCand1);
  varsWithFlag.add(signalFlagCand2);
  RooDataSet datasetWithFlag("datasetWithFlag", "Dataset with Signal Flags", varsWithFlag);

  double sidebandLeft = mean->getVal() - 4 * sigma->getVal();
  double sidebandRight = mean->getVal() + 4 * sigma->getVal();
  double signalLeft = mean->getVal() - 3 * sigma->getVal();
  double signalRight = mean->getVal() + 3 * sigma->getVal();

  // Iterate over the dataset and set the flag
  for (int i = 0; i < dataset->numEntries(); ++i) {

      const RooArgSet* row = dataset->get(i);
      double massValCand1 = row->getRealValue("fMCand1");
      double massValCand2 = row->getRealValue("fMCand2");

      // Determine if the event is signal or background
      if (massValCand1 > signalLeft && massValCand1 < signalRight) {
          // Calculate likelihoods
          //cout << "mass Cand 1 " << massCand1->getVal() << endl;
          massCand1->setVal(massValCand1);

          double signalLikelihoodCand1 = sgnPdfCand1->getVal();
          double backgroundLikelihoodCand1 = bkgPdfCand1->getVal();
          //double modelLikelihood = modelCand1.getVal();
          //cout << "modelCand1: " << modelCand1.getVal() << " backgroundLikelihoodCand1: " << bkgPdfCand1->getVal() << endl;

          // Calculate signal probability using likelihood ratio
          double signalProbCand1 = signalLikelihoodCand1 / (signalLikelihoodCand1 + backgroundLikelihoodCand1);

          // Classify based on probability (e.g., using a threshold)
          if (signalProbCand1 > 0.5) { // Adjust threshold as needed
              signalFlagCand1.setLabel("Signal");
          } else {
              signalFlagCand1.setLabel("Background");
          }
      } else {
          signalFlagCand1.setLabel("Background");
      }
      // Now for cand 2
      if (massValCand2 > signalLeft && massValCand2 < signalRight) {
        // Calculate likelihoods
        massCand2->setVal(massValCand2);
        double signalLikelihoodCand2 = sgnPdfCand2->getVal(&normSet);
        double backgroundLikelihoodCand2 = bkgPdfCand2->getVal(&normSet);

        // Calculate signal probability using likelihood ratio
        double signalProbCand2 = signalLikelihoodCand2 / (signalLikelihoodCand2 + backgroundLikelihoodCand2);
        //cout << signalLikelihoodCand2 << " bkg likelihood: " << backgroundLikelihoodCand2 << endl;

        // Classify based on probability (e.g., using a threshold)
        if (signalProbCand2 > 0.5) { // Adjust threshold as needed
            signalFlagCand2.setLabel("Signal");
        } else {
            signalFlagCand2.setLabel("Background");
        }
    } else {
        signalFlagCand2.setLabel("Background");
    }

      // Add the event with the flag to the new dataset
      varsWithFlag.setRealValue(massCand1->GetName(), massValCand1);
      varsWithFlag.setRealValue(massCand2->GetName(), massValCand2);
      datasetWithFlag.add(varsWithFlag);
  } 
  // Create datasets for signal and background regions based on the flag
  RooDataSet* signalData = (RooDataSet*)datasetWithFlag.reduce(Cut("signalFlagCand1 == signalFlagCand1::Signal"));*/
