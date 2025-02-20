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
#include "RooCurve.h"
#include "TLegend.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include <cstring> // for strcmp
#include "Math/Vector4D.h"
#include "BifurcatedCB.h"

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

// Constructors
InvMassFitter2D::InvMassFitter2D() : _tree(nullptr), _pairType("OS"), nentries(0),
                                     _massMin(1.70), _massMax(2.05), _ptMin(0.), _ptMax(50.),
                                     rooPtCand1("fPtCand1", "pt of candidate 1", 0., 36.),
                                     rooPtCand2("fPtCand2", "pt of candidate 2", 0., 36.),
                                     rooMCand1("fMCand1", "invariant-mass of the first candidate", 1.70, 2.05),
                                     rooMCand2("fMCand2", "invariant-mass of the second candidate", 1.70, 2.05),
                                     rooYCand1("fYCand1", "y of candidate 1", -1., 1.),
                                     rooYCand2("fYCand2", "y of candidate 2", -1., 1.),
                                     rooPhiCand1("fPhiCand1", "phi of candidate 1", -3.5, 3.5),
                                     rooPhiCand2("fPhiCand2", "phi of candidate 2", -3.5, 3.5),
                                     _efficiencyMap(0x0),
                                     _mean(0x0), _sigma(0x0), _meanRefl(0x0), _sigmaRefl(0x0), _tau(0x0), _fracRefl(0x0),
                                     _meanReflDoubleGaus(0x0), _sigmaReflDoubleGaus(0x0), _rawYield(0x0), _reflOverSgn(0),
                                     _workspace(0x0), _bkgPdfCand1(0x0), _sgnPdfCand1(0x0), _reflPdfCand1(0x0), _totPdfCand1(0x0),
                                     _bkgPdfCand2(0x0), _sgnPdfCand2(0x0), _reflPdfCand2(0x0), _totPdfCand2(0x0) {}
InvMassFitter2D::InvMassFitter2D(TTree *tree, const char *pairType) : _tree(tree), _pairType(pairType), nentries(0),
                                                                      _massMin(1.70), _massMax(2.05), _ptMin(0.), _ptMax(50.),
                                                                      rooPtCand1("fPtCand1", "pt of candidate 1", 0., 36.),
                                                                      rooPtCand2("fPtCand2", "pt of candidate 2", 0., 36.),
                                                                      rooMCand1("fMCand1", "invariant-mass of the first candidate", 1.70, 2.05),
                                                                      rooMCand2("fMCand2", "invariant-mass of the second candidate", 1.70, 2.05),
                                                                      rooYCand1("fYCand1", "y of candidate 1", -1., 1.),
                                                                      rooYCand2("fYCand2", "y of candidate 2", -1., 1.),
                                                                      rooPhiCand1("fPhiCand1", "phi of candidate 1", -3.5, 3.5),
                                                                      rooPhiCand2("fPhiCand2", "phi of candidate 2", -3.5, 3.5),
                                                                      _efficiencyMap(0x0),
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

void InvMassFitter2D::setPtLims(double const &ptMin, double const &ptMax)
{
  _ptMin = ptMin;
  _ptMax = ptMax;
}

void InvMassFitter2D::setLumi(double const &lumi)
{
  _lumi = lumi;
}

void InvMassFitter2D::setMassLims(double const &massMin, double const &massMax)
{
  _massMin = massMin;
  _massMax = massMax;
}

void InvMassFitter2D::setEfficiencyMap(TH2F *h)
{
  _efficiencyMap = (TH2F *)h->Clone();
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
  RooArgSet vars(rooPtCand1, rooPtCand2, rooMCand1, rooMCand2, rooYCand1, rooYCand2, rooPhiCand1, rooPhiCand2); // Cand 1: D, Cand 2: Dbar (OS)

  RooRealVar weightCand1("weightCand1", "weights of cand 1", 1., 0., 100.);
  RooRealVar weightCand2("weightCand2", "weights of cand 2", 1., 0., 100.);
  RooArgSet weightedVars(rooPtCand1, rooPtCand2, rooMCand1, rooMCand2, rooYCand1, rooYCand2, rooPhiCand1, rooPhiCand2, weightCand1, weightCand2);
  RooFormulaVar combinedWeight("combinedWeight", "combined weight", "weightCand1 * weightCand2", RooArgList(weightCand1, weightCand2));

  // Create an empty dataset with the variables and category
  RooDataSet data("data", "data", vars);
  RooDataSet weightedData("weightedData", "weightedData", weightedVars, RooFit::WeightVar(combinedWeight.GetName()));

  fillDataset(data, vars);
  fillDataset(weightedData, weightedVars);
}

void InvMassFitter2D::fillDataset(RooDataSet &data, RooArgSet &vars)
{

  RooRealVar *weightsCand1;
  RooRealVar *weightsCand2;

  double totalWeight = 0.; // Variable to store the total weight for normalization
  if (_efficiencyMap && vars.getSize() != 6)
  {
    weightsCand1 = dynamic_cast<RooRealVar *>(vars.find("weightCand1"));
    weightsCand2 = dynamic_cast<RooRealVar *>(vars.find("weightCand2"));
  }

  // Fill the dataset with info from the tree
  cout << "Number of tree entries " << _tree->GetEntries() << endl;
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

    // Select pT range
    if ((ptCand1 < _ptMin || ptCand2 < _ptMin) || (ptCand1 > _ptMax || ptCand2 > _ptMax))
    {
      continue;
    }

    // Add cuts to avoid ambiguous candidates
    if (TESTBIT(typeCand1, SelectedD) && TESTBIT(typeCand1, SelectedDbar))
    {
      continue;
    }
    if (TESTBIT(typeCand2, SelectedD) && TESTBIT(typeCand2, SelectedDbar))
    {
      continue;
    }

    ROOT::Math::PxPyPzMVector vLorentzCand1 = createLorentzVector(phiCand1, yCand1, ptCand1, mDCand1);
    ROOT::Math::PxPyPzMVector vLorentzCand2 = createLorentzVector(phiCand2, yCand2, ptCand2, mDCand2);
    ROOT::Math::PxPyPzMVector vLorentzPair = vLorentzCand1 + vLorentzCand2;

    TH1F* ptPairHist = new TH1F("ptPairHist", "Pair Pt distribution", 100, 0, 10);
    ptPairHist->Fill(vLorentzPair.Pt());
    ptPairHist->Draw();

    if (TESTBIT(typePair, DD))
    {
      if ((mDCand1 < _massMin || mDCand1 > _massMax) || (mDCand2 < _massMin || mDCand2 > _massMax))
        continue;

      rooMCand1.setVal(mDCand1);
      rooMCand2.setVal(mDCand2);
    }
    if (TESTBIT(typePair, DbarDbar))
    {
      if ((mDbarCand1 < _massMin || mDbarCand1 > _massMax) || (mDbarCand2 < _massMin || mDbarCand2 > _massMax))
        continue;
      rooMCand1.setVal(mDbarCand1);
      rooMCand2.setVal(mDbarCand2);
    }
    if (TESTBIT(typePair, DDbar))
    {
      if ((mDCand1 < _massMin || mDCand1 > _massMax) || (mDbarCand2 < _massMin || mDbarCand2 > _massMax))
        continue;
      rooMCand1.setVal(mDCand1);
      rooMCand2.setVal(mDbarCand2);
    }
    if (TESTBIT(typePair, DbarD))
    {
      if ((mDbarCand1 < _massMin || mDbarCand1 > _massMax) || (mDCand2 < _massMin || mDCand2 > _massMax))
        continue;
      rooMCand1.setVal(mDCand2);
      rooMCand2.setVal(mDbarCand1);
    }

    rooPtCand1.setVal(ptCand1);
    rooPtCand2.setVal(ptCand2);
    rooYCand1.setVal(yCand1);
    rooYCand2.setVal(yCand2);
    rooPhiCand1.setVal(yCand1);
    rooPhiCand2.setVal(yCand2);

    double weightCand1 = 1., weightCand2 = 1.;
    double combinedWeight = 1.;

    totalWeight += combinedWeight;

    if (_efficiencyMap)
    {
      weightCand1 = calculateWeights(yCand1, ptCand1);
      weightCand2 = calculateWeights(yCand2, ptCand2);
      combinedWeight = weightCand1 * weightCand2;
      if (vars.getSize() != 6)
      {
        weightsCand1->setVal(weightCand1);
        weightsCand2->setVal(weightCand2);
      }
    }

    if (strcmp(_pairType, "OS") == 0)
    {
      if (TESTBIT(typePair, DDbar))
      {
        if (vars.getSize() != 6)
          data.add(vars, combinedWeight);
        else
          data.add(vars);
      }
      if (TESTBIT(typePair, DbarD))
      {
        if (vars.getSize() != 6)
          data.add(vars, combinedWeight);
        else
          data.add(vars);
      }
    }
    else if (strcmp(_pairType, "LS") == 0)
    {
      if (TESTBIT(typePair, DD))
      {
        if (vars.getSize() != 6)
          data.add(vars, combinedWeight);
        else
          data.add(vars);
      }
      if (TESTBIT(typePair, DbarDbar))
      {
        if (vars.getSize() != 6)
          data.add(vars, combinedWeight);
        else
          data.add(vars);
      }
    }
    else
    {
      cerr << "ERROR: wrong pairType assigned. Please choose OS or LS" << endl;
      return;
    }
  }
  _workspace.import(data);
  cout << "data loaded" << endl;
}

void InvMassFitter2D::set1DParameters(const RooArgSet *vars1D, double const &reflOverSgn, double const &integratedEfficiency)
{
  _tau = (RooRealVar *)vars1D->find("tau");
  _mean = (RooRealVar *)vars1D->find("mean");
  _sigma = (RooRealVar *)vars1D->find("sigma");
  _meanRefl = (RooRealVar *)vars1D->find("meanRefl");
  _sigmaRefl = (RooRealVar *)vars1D->find("sigmaRefl");
  _fracRefl = (RooRealVar *)vars1D->find("fracRefl");
  _meanReflDoubleGaus = (RooRealVar *)vars1D->find("meanReflDoubleGaus");
  _sigmaReflDoubleGaus = (RooRealVar *)vars1D->find("sigmaReflDoubleGaus");
  _rawYield = (RooRealVar *)vars1D->find("rooRawYield");
  _reflOverSgn = reflOverSgn;
  _integratedEfficiency = integratedEfficiency;
}

/// @brief Fill workspace with relevant functions for signal, bkg and reflections
void InvMassFitter2D::fillWorkspace(RooDataSet *dataset)
{
  const RooArgSet *vars = dataset->get();
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
  RooRealVar tauCand1("tauCand1", "tauCand1", -2, -10., -1.);
  RooRealVar tauCand2("tauCand2", "tauCand2", -2, -10., -1.);

  RooAbsPdf *bkgFuncExpoCand1 = new RooExponential("bkgFuncExpoCand1", "background fit function of candidate 1", *massCand1, tauCand1);
  RooAbsPdf *bkgFuncExpoCand2 = new RooExponential("bkgFuncExpoCand2", "background fit function of candidate 2", *massCand2, tauCand2);
  _workspace.import(*bkgFuncExpoCand1);
  _workspace.import(*bkgFuncExpoCand2);

  // bkg poly1 (useful for MC)
  RooRealVar c1("c1", "Linear coefficient", 1., 0., 5.);      // c1: linear coefficient
  RooRealVar c1Cand2("c1Cand2", "Linear coefficient", 1., 0., 5.);      // c1: linear coefficient

  RooAbsPdf *bkgFuncPoly1Cand1 = new RooPolynomial("bkgFuncPoly1Cand1", "background fit function of candidate 1", *massCand1, RooArgList(c1));
  RooAbsPdf *bkgFuncPoly1Cand2 = new RooPolynomial("bkgFuncPoly1Cand2", "background fit function of candidate 2", *massCand2, RooArgList(c1Cand2));
  _workspace.import(*bkgFuncPoly1Cand1);
  _workspace.import(*bkgFuncPoly1Cand2);

  // bkg poly2
  RooRealVar c2("c2", "Quadratic coefficient", -1., -5., 0.); // c2: quadratic coefficient
  RooRealVar c2Cand2("c2Cand2", "Quadratic coefficient", -1., -5., 0.); // c2: quadratic coefficient

  RooAbsPdf *bkgFuncPoly2Cand1 = new RooPolynomial("bkgFuncPoly2Cand1", "background fit function of candidate 1", *massCand1, RooArgList(c1, c2));
  RooAbsPdf *bkgFuncPoly2Cand2 = new RooPolynomial("bkgFuncPoly2Cand2", "background fit function of candidate 2", *massCand2, RooArgList(c1Cand2, c2Cand2));
  _workspace.import(*bkgFuncPoly2Cand1);
  _workspace.import(*bkgFuncPoly2Cand2);

  // bkg poly3
  RooRealVar c3("c3", "Cubic coefficient", 1., -3., 3.);           // c3: cubic coefficient
  RooRealVar c3Cand2("c3Cand2", "Cubic coefficient", 1., -3., 3.); // c3: cubic coefficient

  RooAbsPdf *bkgFuncPoly3Cand1 = new RooPolynomial("bkgFuncPoly3Cand1", "background fit function of candidate 1", *massCand1, RooArgList(c1, c2, c3));
  RooAbsPdf *bkgFuncPoly3Cand2 = new RooPolynomial("bkgFuncPoly3Cand2", "background fit function of candidate 2", *massCand2, RooArgList(c1Cand2, c2Cand2, c3Cand2));
  _workspace.import(*bkgFuncPoly3Cand1);
  _workspace.import(*bkgFuncPoly3Cand2);

  /// | ------------------------------------------------------------------ |
  /// | ----------------------- SIGNAL FUNCTIONS ------------------------- |
  /// | ------------------------------------------------------------------ |
  // signal pdf
  RooRealVar mean("mean", "mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigma("sigma", "sigma for signal", 0.02, 0.01, 0.07);

  RooRealVar meanCand2("meanCand2", "mean for signal fit", 1.85, 1.83, 1.9);
  RooRealVar sigmaCand2("sigmaCand2", "sigma for signal", 0.02, 0.01, 0.07);

  RooAbsPdf *sgnFuncGausCand1 = new RooGaussian("sgnFuncGausCand1", "signal pdf of candidate 1", *massCand1, mean, sigma);
  RooAbsPdf *sgnFuncGausCand2 = new RooGaussian("sgnFuncGausCand2", "signal pdf of candidate 2", *massCand2, meanCand2, sigmaCand2);
  _workspace.import(*sgnFuncGausCand1);
  _workspace.import(*sgnFuncGausCand2);

  RooRealVar alpha1("alpha1", "Alpha (transition)", 1.5, 0.5, 2.0);
  RooRealVar n1("n1", "n (steepness of tail)", 5.0, 0.5, 3.0);
  RooRealVar alpha2("alpha2", "Alpha (transition)", 1.5, 0.5, 2.0);
  RooRealVar n2("n2", "n (steepness of tail)", 5.0, 0.5, 20.0);

  /* RooAbsPdf *sgnFuncGausCand1 = new RooGaussian("sgnFuncGausCand1", "signal pdf of candidate 1", *massCand1, mean, sigma);
  RooAbsPdf *sgnFuncGausCand2 = new RooGaussian("sgnFuncGausCand2", "signal pdf of candidate 2", *massCand2, mean, sigma); */
  RooAbsPdf *sgnFuncCBCand1 = new BifurcatedCB("sgnFuncCBCand1", "signal pdf of candidate 1", *massCand1, mean, sigma, alpha1, n1, alpha2, n2);
  RooAbsPdf *sgnFuncCBCand2 = new BifurcatedCB("sgnFuncCBCand2", "signal pdf of candidate 2", *massCand2, meanCand2, sigmaCand2, alpha1, n1, alpha2, n2);
  _workspace.import(*sgnFuncCBCand1);
  _workspace.import(*sgnFuncCBCand2);

  /// | ------------------------------------------------------------------ |
  /// | ------------------- REFLECTION FUNCTIONS ------------------------- |
  /// | ------------------------------------------------------------------ |
  // reflection Gaussian
  RooRealVar meanRefl("meanRefl", "mean for reflections", 1.85, 0.0, 2.15);
  RooRealVar sigmaRefl("sigmaRefl", "sigma for reflection", 0.012, 0, 0.3);

  RooAbsPdf *reflFuncGausCand1 = new RooGaussian("reflFuncGausCand1", "reflection pdf of candidate 1", *massCand1, meanRefl, sigmaRefl);
  RooAbsPdf *reflFuncGausCand2 = new RooGaussian("reflFuncGausCand2", "reflection pdf of candidate 2", *massCand2, meanRefl, sigmaRefl);
  _workspace.import(*reflFuncGausCand1);
  _workspace.import(*reflFuncGausCand2);

  // reflection double gaussian
  RooRealVar meanReflDoubleGaus("meanReflDoubleGaus", "mean for reflection double gaussian", 1.85, 0.0, 1.90);
  RooRealVar sigmaReflDoubleGaus("sigmaReflDoubleGaus", "sigmaReflDoubleGaus", 0.012, 0.0, 0.2);
  RooRealVar fracRefl("fracRefl", "frac of the two reflected gaussians", 0.5, 0, 1.);
  //    Second gaussian definition
  RooAbsPdf *relfFuncSecondGausCand1 = new RooGaussian("relfFuncSecondGausCand1", "relfFuncSecondGausCand1", *massCand1, meanReflDoubleGaus, sigmaReflDoubleGaus);
  RooAbsPdf *relfFuncSecondGausCand2 = new RooGaussian("relfFuncSecondGausCand2", "relfFuncSecondGausCand2", *massCand2, meanReflDoubleGaus, sigmaReflDoubleGaus);
  //    Compose the final double gaussian
  RooAbsPdf *reflFuncDoubleGausCand1 = new RooAddPdf("reflFuncDoubleGausCand1", "reflection pdf of candidate 1", RooArgList(*reflFuncGausCand1, *relfFuncSecondGausCand1), fracRefl);
  RooAbsPdf *reflFuncDoubleGausCand2 = new RooAddPdf("reflFuncDoubleGausCand2", "reflection pdf of candidate 2", RooArgList(*reflFuncGausCand2, *relfFuncSecondGausCand2), fracRefl);
  _workspace.import(*reflFuncDoubleGausCand1);
  _workspace.import(*reflFuncDoubleGausCand2);

  cout << "Workspace filled with functions" << endl;
}

void InvMassFitter2D::do2DFit(Bool_t draw, Bool_t doReflections, Bool_t isMc, TFile *fout)
{
  createDataset();
  // Declare observable variable
  _workspace.Print("v");
  RooDataSet *dataset = dynamic_cast<RooDataSet *>(_workspace.data("data"));
  RooDataSet *weightedDataset = dynamic_cast<RooDataSet *>(_workspace.data("weightedData"));

  if (!dataset)
  {
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

  RooAbsPdf *bkgPdfCand1;
  RooAbsPdf *bkgPdfCand2;
  RooAbsPdf *sgnPdfCand1;
  RooAbsPdf *sgnPdfCand2;
  RooAbsPdf *reflPdfCand1;
  RooAbsPdf *reflPdfCand2;

  // Select signal function
  if (_sgnFuncOption == "gaus")
  {
    cout << "Signal function chosen: GAUSSIAN" << endl;
    sgnPdfCand1 = _workspace.pdf("sgnFuncGausCand1");
    sgnPdfCand2 = _workspace.pdf("sgnFuncGausCand2");
  }
  else if (_sgnFuncOption == "CB")
  {
    cout << "Signal function chosen: BIFURCATED CRYSTAL BALL" << endl;
    sgnPdfCand1 = _workspace.pdf("sgnFuncCBCand1");
    sgnPdfCand2 = _workspace.pdf("sgnFuncCBCand2");
  }
  else
  {
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
    bkgPdfCand1 = _workspace.pdf("bkgFuncExpoCand1");
    bkgPdfCand2 = _workspace.pdf("bkgFuncExpoCand2");
  }
  else if (_bkgFuncOption == "poly1")
  {
    cout << "Background function chosen: POLY 1" << endl;
    bkgPdfCand1 = _workspace.pdf("bkgFuncPoly1Cand1");
    bkgPdfCand2 = _workspace.pdf("bkgFuncPoly1Cand2");
  }
  else if (_bkgFuncOption == "poly2")
  {
    cout << "Background function chosen: POLY 2" << endl;
    bkgPdfCand1 = _workspace.pdf("bkgFuncPoly2Cand1");
    bkgPdfCand2 = _workspace.pdf("bkgFuncPoly2Cand2");
  }
  else if (_bkgFuncOption == "poly3")
  {
    cout << "Background function chosen: POLY 3" << endl;
    bkgPdfCand1 = _workspace.pdf("bkgFuncPoly3Cand1");
    bkgPdfCand2 = _workspace.pdf("bkgFuncPoly3Cand2");
  }
  else
  {
    cerr << "ERROR: background function not supported! \n Available options: expo, poly1, poly2, poly3. \n Exit!" << endl;
    return;
  }

  if (!bkgPdfCand1 || !bkgPdfCand2)
  {
    cerr << "ERROR: bkgPdf function not found!" << endl;
    return;
  }

  // Select reflection function
  if (_reflFuncOption == "gaus")
  {
    cout << "Reflected function chosen: GAUSSIAN" << endl;
    reflPdfCand1 = _workspace.pdf("reflFuncGausCand1");
    reflPdfCand2 = _workspace.pdf("reflFuncGausCand2");
  }
  else if (_reflFuncOption == "doubleGaus")
  {
    cout << "Reflected function chosen: DOUBLE GAUSSIAN" << endl;
    reflPdfCand1 = _workspace.pdf("reflFuncDoubleGausCand1");
    reflPdfCand2 = _workspace.pdf("reflFuncDoubleGausCand2");
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

  /// ----------------------------------------------------------
  /// ------------ LOAD AND FIX PARAMETERS ---------------------
  /// ----------------------------------------------------------
  //  Most parameters are just set as initial values

  RooRealVar *mean = _workspace.var("mean");
  RooRealVar *sigma = _workspace.var("sigma");
  RooRealVar *tauCand1 = _workspace.var("tauCand1");
  RooRealVar *tauCand2 = _workspace.var("tauCand2");

  mean->setVal(_mean->getVal());
  mean->setError(_mean->getError());
  sigma->setVal(_sigma->getVal());
  sigma->setError(_sigma->getError());
  tauCand1->setVal(_tau->getVal());
  tauCand1->setError(_tau->getError());
  tauCand2->setVal(_tau->getVal());
  tauCand2->setError(_tau->getError());

  // mean->setConstant(kTRUE);
  // sigma->setConstant(kTRUE);

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

  if (doReflections)
  {
    meanRefl->setVal(_meanRefl->getVal());
    meanRefl->setError(_meanRefl->getError());
    sigmaRefl->setVal(_sigmaRefl->getVal());
    sigmaRefl->setError(_sigmaRefl->getError());

    cout << "sigmaRefl: " << sigmaRefl->getVal() << endl;
    cout << "_sigmaRefl: " << _sigmaRefl->getVal() << endl;

    meanReflDoubleGaus->setVal(_meanReflDoubleGaus->getVal());
    meanReflDoubleGaus->setError(_meanReflDoubleGaus->getError());
    sigmaReflDoubleGaus->setVal(_sigmaReflDoubleGaus->getVal());
    sigmaReflDoubleGaus->setError(_sigmaReflDoubleGaus->getError());
    fracRefl->setVal(_fracRefl->getVal());

    /* meanRefl->setConstant(kTRUE);
    sigmaRefl->setConstant(kTRUE);
    fracRefl->setConstant(kTRUE);
    meanReflDoubleGaus->setConstant(kTRUE);
    sigmaReflDoubleGaus->setConstant(kTRUE); */
  }

  /// -----------------------------------------------------------------
  /// --------------- DEFINE YIELD PARAMETERS FOR FIT -----------------
  /// -----------------------------------------------------------------

  RooArgSet normSet(*massCand1, *massCand2);
  RooArgSet weightedNormSet(*massCand1, *massCand2);

  // RAW PARAMETERS
  RooRealVar *nBkg1 = new RooRealVar("nBkg1", "background yield of cand 1", 100.0, 0.0, 500.0);
  RooRealVar *nSgn1 = new RooRealVar("nSgn1", "signal yield of cand1", 30.0, 0.0, 300.0);
  RooFormulaVar *nRefl1 = new RooFormulaVar("nRefl1", "reflected signal yield of cand1", "@0 * @1", RooArgList(_reflOverSgn, *nSgn1));

  RooRealVar *nBkg2 = new RooRealVar("nBkg2", "background yield of cand2", 100.0, 0.0, 500.0);
  RooRealVar *nSgn2 = new RooRealVar("nSgn2", "signal yield of cand2", 30.0, 0.0, 300.0);
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
  RooRealVar *nBkgCorr1 = new RooRealVar("nBkgCorr1", "Corrected background yield of cand 1", 9000, 1000, 20000);
  RooRealVar *nSgnCorr1 = new RooRealVar("nSgnCorr1", "Corrected signal yield of cand1", 1000, 100.0, 5000.0);
  RooFormulaVar *nReflCorr1 = new RooFormulaVar("nReflCorr1", "Corrected reflected signal yield of cand1", "@0 * @1", RooArgList(_reflOverSgn, *nSgnCorr1));

  RooRealVar *nBkgCorr2 = new RooRealVar("nBkgCorr2", "Corrected background yield of cand2", 9000, 1000, 20000);
  RooRealVar *nSgnCorr2 = new RooRealVar("nSgnCorr2", "Corrected signal yield of cand2", 1000, 100.0, 5000.0);
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

  // There isn't a direct equivalent of RooPlot for >1 dimensions
  TH2D *hMassCorrelations = new TH2D("hMassCorrelations", "data of mass correlations", 45, _massMin, _massMax, 45, _massMin, _massMax);
  dataset->fillHistogram(hMassCorrelations, RooArgList(*massCand1, *massCand2));

  TH2D *hMassWeightedCorrelations = new TH2D("hMassWeightedCorrelations", "data of weighted mass correlations", 45, _massMin, _massMax, 45, _massMin, _massMax);
  weightedDataset->fillHistogram(hMassWeightedCorrelations, RooArgList(*massCand1, *massCand2));

  // Define the difference between nSgn1 and nSgn2
  RooFormulaVar diffSgn("diffSgn", "difference between nSgn1 and nSgn2", "@0 - @1", RooArgList(*nSgn1, *nSgn2));
  RooFormulaVar diffBkg("diffBkg", "difference between nBkg1 and nBkg2", "@0 - @1", RooArgList(*nBkg1, *nBkg2));
  RooFormulaVar diffSgnCorr("diffSgnCorr", "difference between nSgnCorr1 and nSgnCorr2", "@0 - @1", RooArgList(*nSgnCorr1, *nSgnCorr2));
  RooFormulaVar diffBkgCorr("diffBkgCorr", "difference between nBkgCorr1 and nBkgCorr2", "@0 - @1", RooArgList(*nBkgCorr1, *nBkgCorr2));

  // Define the Gaussian constraint to keep nSgn1 and nSgn2 symmetric
  RooRealVar sigmaConstraint("sigmaConstraint", "width of Gaussian constraint", 50.0, 1.0, 100.0); // adjust this as needed
  sigmaConstraint.setConstant(kTRUE);
  RooRealVar sigmaConstraintCorr("sigmaConstraintCorr", "width of Gaussian constraint corr", 7.0, 1.0, 100.0); // adjust this as needed
  sigmaConstraintCorr.setConstant(kTRUE);
  RooGaussian *constraintSgn = new RooGaussian("constraintSgn", "Gaussian constraint on nSgn1 - nSgn2", diffBkg, RooConst(0), sigmaConstraint);
  RooGaussian *constraintSgnCorr = new RooGaussian("constraintSgnCorr", "Gaussian constraint on nSgnCorr1 - nSgnCorr2", diffSgnCorr, RooConst(0), sigmaConstraintCorr);
  RooGaussian *constraintBkgCorr = new RooGaussian("constraintBkgCorr", "Gaussian constraint on nBkgCorr1 - nBkgCorr2", diffBkgCorr, RooConst(0), sigmaConstraintCorr);
  // RooArgSet constraintsCorr(*constraintSgnCorr, *constraintBkgCorr);

  // Final fit
  if (isMc)
  {
    if (doReflections)
    {
      _totPdf2D = new RooAddPdf("_totPdf2D", "signal + reflection pdf 2D",
        RooArgList(sgnSgnFunc2D, sgnReflFunc2D, reflSgnFunc2D, sgnBkgFunc2D,
                   bkgSgnFunc2D, bkgBkgFunc2D, bkgReflFunc2D, reflBkgFunc2D, reflReflFunc2D),
        RooArgList(nSgnSgn, nSgnRefl, nReflSgn, nSgnBkg, nBkgSgn, nBkgBkg, nBkgRefl, nReflBkg, nReflRefl));

      _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "Corrected signal + reflection pdf 2D",
                RooArgList(sgnSgnFunc2D, sgnReflFunc2D, reflSgnFunc2D, reflReflFunc2D),
                RooArgList(nSgnSgnCorr, nSgnReflCorr, nReflSgnCorr, nReflReflCorr));
    }
    else
    {
      _totPdf2D = new RooAddPdf("_totPdf2D", "signal pdf 2D",
                                RooArgList(sgnSgnFunc2D),
                                RooArgList(nSgnSgn));

      _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "Corrected signal pdf 2D",
                                        RooArgList(sgnSgnFunc2D),
                                        RooArgList(nSgnSgnCorr));
    }
  }
  else
  {
    if (doReflections)
    {
      _totPdf2D = new RooAddPdf("_totPdf2D", "background + signal + reflection pdf 2D",
                                RooArgList(sgnSgnFunc2D, sgnReflFunc2D, reflSgnFunc2D, sgnBkgFunc2D,
                                           bkgSgnFunc2D, bkgBkgFunc2D, bkgReflFunc2D, reflBkgFunc2D, reflReflFunc2D),
                                RooArgList(nSgnSgn, nSgnRefl, nReflSgn, nSgnBkg, nBkgSgn, nBkgBkg, nBkgRefl, nReflBkg, nReflRefl));

      _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "Corrected background + signal + reflection pdf 2D",
                                        RooArgList(sgnSgnFunc2D, sgnReflFunc2D, reflSgnFunc2D, sgnBkgFunc2D, bkgSgnFunc2D,
                                                   bkgBkgFunc2D, bkgReflFunc2D, reflBkgFunc2D, reflReflFunc2D),
                                        RooArgList(nSgnSgnCorr, nSgnReflCorr, nReflSgnCorr, nSgnBkgCorr, nBkgSgnCorr, nBkgBkgCorr,
                                                   nBkgReflCorr, nReflBkgCorr, nReflReflCorr));
    }
    else
    {
      _totPdf2D = new RooAddPdf("_totPdf2D", "background + signal pdf 2D",
                                RooArgList(bkgBkgFunc2D, sgnSgnFunc2D, sgnBkgFunc2D, bkgSgnFunc2D),
                                RooArgList(nBkgBkg, nSgnSgn, nSgnBkg, nBkgSgn));

      _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "background + signal pdf 2D",
                                        RooArgList(bkgBkgFunc2D, sgnSgnFunc2D, sgnBkgFunc2D, bkgSgnFunc2D),
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
                                                               RooFit::SumW2Error(kTRUE),
                                                               RooFit::Minimizer("Minuit2", "Migrad"),
                                                               RooFit::Strategy(0),
                                                               RooFit::MaxCalls(10000), // Adjust max calls as needed
                                                               RooFit::PrintLevel(-1)   // Suppress verbose output
  );
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

  // sigma->setConstant(kTRUE);
  // mean->setConstant(kTRUE);

  if (doReflections)
  {
    fracRefl->setConstant(kTRUE);
    meanRefl->setConstant(kTRUE);
    sigmaRefl->setConstant(kTRUE);
    meanReflDoubleGaus->setConstant(kTRUE);
    sigmaReflDoubleGaus->setConstant(kTRUE);
  }

  RooFitResult *resultConstrainedCorr = modelWithConstraintCorr->fitTo(*weightedDataset, RooFit::Save(),
                                                                       RooFit::SumW2Error(kTRUE),
                                                                       RooFit::Minimizer("Minuit2", "Migrad"),
                                                                       RooFit::Strategy(0),
                                                                       RooFit::MaxCalls(10000), // Adjust max calls as needed
                                                                       RooFit::PrintLevel(-1)   // Suppress verbose output
  );

  if (resultConstrainedCorr->status() == 0)
  {
    std::cout << "Fit Converged!" << std::endl;
  }
  else
  {
    std::cout << "Fit did not converge." << std::endl;
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

  RooRealVar *yieldSgnRefl;
  RooRealVar *yieldReflSgn;
  RooRealVar *yieldBkgRefl;
  RooRealVar *yieldReflBkg;
  RooRealVar *yieldReflRefl;
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

  // ------------------------------------------------------
  // -- Get efficiency corrected yields per contribution --
  // ------------------------------------------------------

  modelWithConstraintCorr->getParameters(*massCand1)->assignValueOnly(resultConstrainedCorr->floatParsFinal());
  modelWithConstraintCorr->getParameters(*massCand2)->assignValueOnly(resultConstrainedCorr->floatParsFinal());

  cout << "\n\n ||||| EFFICIENCY CORRECTED INTEGRAL VALUES |||||" << endl;
  cout << "Signal x Signal ";
  RooRealVar *yieldSgnSgnCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, sgnSgnFunc2D, nSgnSgnCorr, "3sigmaRange");

  cout << "Signal x Background ";
  RooRealVar *yieldSgnBkgCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, sgnBkgFunc2D, nSgnBkgCorr, "3sigmaRange");

  cout << "Background x Signal ";
  RooRealVar *yieldBkgSgnCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, bkgSgnFunc2D, nBkgSgnCorr, "3sigmaRange");

  cout << "Background x Background ";
  RooRealVar *yieldBkgBkgCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, bkgBkgFunc2D, nBkgBkgCorr, "3sigmaRange");

  RooRealVar *yieldSgnReflCorr;
  RooRealVar *yieldReflSgnCorr;
  RooRealVar *yieldBkgReflCorr;
  RooRealVar *yieldReflBkgCorr;
  RooRealVar *yieldReflReflCorr;
  if (doReflections)
  {

    cout << "Signal x Reflected ";
    yieldSgnReflCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, sgnReflFunc2D, nSgnReflCorr, "3sigmaRange");

    cout << "Reflected x Signal ";
    yieldReflSgnCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, reflSgnFunc2D, nReflSgnCorr, "3sigmaRange");

    cout << "Background x Reflected ";
    yieldBkgReflCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, bkgReflFunc2D, nBkgReflCorr, "3sigmaRange");

    cout << "Reflected x Background ";
    yieldReflBkgCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, reflBkgFunc2D, nReflBkgCorr, "3sigmaRange");

    cout << "Reflected x Reflected ";
    yieldReflReflCorr = getYieldInRange(resultConstrainedCorr, massCand1, massCand2, reflReflFunc2D, nReflReflCorr, "3sigmaRange");
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

  float const BR = 0.03950;
  float const BRErr = 0.0003;

  float xSec = yieldSgnSgnCorr->getVal() / (BR * BR * _lumi);
  cout << "CROSS SECTION:  " << xSec << endl;

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
    plotProjectionsAfterFit(modelWithConstraint, dataset, "cRawMassProjections", fout, doReflections);
    plotProjectionsAfterFit(modelWithConstraintCorr, weightedDataset, "cWeightedMassProjections", fout, doReflections);
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

  analyseKinematicDistributions(fout);

  // Create a frame for ptCand1 with a specified range and number of bins
  RooPlot *framePtCand1 = ptCand1->frame(RooFit::Bins(50), RooFit::Range(ptCand1->getMin(), ptCand1->getMax()));
  RooPlot *framePtCand2 = ptCand2->frame(RooFit::Bins(50), RooFit::Range(ptCand2->getMin(), ptCand2->getMax()));

  RooPlot *frameWeightedPtCand1 = ptCand1->frame(RooFit::Bins(50), RooFit::Range(ptCand1->getMin(), ptCand1->getMax()));
  RooPlot *frameWeightedPtCand2 = ptCand2->frame(RooFit::Bins(50), RooFit::Range(ptCand2->getMin(), ptCand2->getMax()));

  dataset->plotOn(framePtCand1);
  dataset->plotOn(framePtCand2);
  weightedDataset->plotOn(frameWeightedPtCand1);
  weightedDataset->plotOn(frameWeightedPtCand2);

  // Draw the frame on a canvas
  TCanvas *cPtCand1 = new TCanvas("cPtCand1", "pt of Cand1 vs Number of Candidates", 800, 600);
  cPtCand1->SetLogy();
  framePtCand1->Draw();
  framePtCand1->SetTitle("pT of cand1 BEFORE correcting by efficiencies");
  cPtCand1->Write();

  TCanvas *cPtCand2 = new TCanvas("cPtCand2", "pt of Cand2 vs Number of Candidates", 800, 600);
  cPtCand2->SetLogy();
  framePtCand2->Draw();
  framePtCand2->SetTitle("pT of cand2 BEFORE correcting by efficiencies");
  cPtCand2->Write();

  TCanvas *cWeightedPtCand1 = new TCanvas("cWeightedPtCand1", "Weighted pt of Cand1 vs Number of Candidates", 800, 600);
  cWeightedPtCand1->SetLogy();
  frameWeightedPtCand1->Draw();
  frameWeightedPtCand1->SetTitle("pT of cand1 AFTER correcting by efficiencies");
  cWeightedPtCand1->Write();

  TCanvas *cWeightedPtCand2 = new TCanvas("cWeightedPtCand2", "Weighted pt of Cand2 vs Number of Candidates", 800, 600);
  cWeightedPtCand2->SetLogy();
  frameWeightedPtCand2->Draw();
  frameWeightedPtCand2->SetTitle("pT of cand2 AFTER correcting by efficiencies");
  cWeightedPtCand2->Write();

  cPtCand1->Close();
  cPtCand2->Close();
  cWeightedPtCand1->Close();
  cWeightedPtCand2->Close();
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
    histFit->SetLineColor(kRed);
    histFit->SetLineWidth(2);
    histFit->SetTitle("");
    /* histFit->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
    histFit->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
    histFit->GetZaxis()->SetTitle("Counts per (3.5 MeV/#it{c}^{2})^{2}");
    histFit->GetXaxis()->CenterTitle();
    histFit->GetYaxis()->CenterTitle(); */
    histFit->GetXaxis()->SetDecimals(2);
    histFit->GetYaxis()->SetDecimals(2);

    hMassCorrelations->GetXaxis()->SetDecimals(2);
    hMassCorrelations->GetYaxis()->SetDecimals(2);
    hMassCorrelations->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
    hMassCorrelations->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
    hMassCorrelations->GetZaxis()->SetTitle("Counts per (3.5 MeV/#it{c}^{2})^{2}");
    hMassCorrelations->GetXaxis()->CenterTitle();
    hMassCorrelations->GetYaxis()->CenterTitle();
    hMassCorrelations->GetXaxis()->SetDecimals(2);
    hMassCorrelations->GetYaxis()->SetDecimals(2);
    hMassCorrelations->GetXaxis()->SetTitleOffset(2);    // Adjust X-axis title offset
    hMassCorrelations->GetYaxis()->SetTitleOffset(2);    // Adjust Y-axis title offset
    hMassCorrelations->GetZaxis()->SetTitleOffset(1.55); // Adjust Z-axis title offset
    hMassCorrelations->SetTitle("");

    TPaveText *ptTitle = new TPaveText(0.12, 0.85, 0.4, 0.94, "brNDC");
    ptTitle->SetTextFont(42);
    ptTitle->SetTextSize(0.04);
    ptTitle->SetBorderSize(0);
    ptTitle->SetFillStyle(0);
    ptTitle->SetTextAlign(11);

    // Add text lines
    ptTitle->AddText("#scale[1.35]{ALICE Work in progress}");
    ptTitle->AddText("#scale[1.2]{pp, #sqrt{#it{s}} = 13.6 TeV}");

    TPaveText *pt = new TPaveText(0.12, 0.69, 0.4, 0.85, "brNDC");
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
    pt->AddText(Form("%.1f < #it{p}_{T} < %.1f GeV/#it{c}", _ptMin, _ptMax));

    gStyle->SetOptStat(0);
    hMassCorrelations->Draw("LEGO2");
    histFit->Draw("sameLEGO0");

    ptTitle->Draw("same");
    pt->Draw("same");
    fout->cd();
    c2d->Write();

    TString c1dName = cName + "_colz";
    TCanvas *c1d = new TCanvas(c1dName, c1dName, 1000, 700);
    hMassCorrelations->Draw("colz");
    c1d->Write();
    c1d->Close();
  }
}

void InvMassFitter2D::plotProjectionsAfterFit(RooProdPdf *model, RooDataSet *dataset, TString saveName, TFile *fout, bool doReflections)
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

  dataset->plotOn(frameCand1, RooFit::Name("data_cand1")); // plot data histogram on the frame
  dataset->plotOn(frameCand2, RooFit::Name("data_cand2")); // plot data histogram on the frame

  TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  TLegend *legend_2 = new TLegend(0.6, 0.7, 0.9, 0.9);

  RooCurve *sgnReflFunc2DCurve;
  RooCurve *bkgReflFunc2DCurve;
  RooCurve *reflBkgFunc2DCurve;
  RooCurve *reflSgnFunc2DCurve;
  RooCurve *reflReflFunc2DCurve;
  if (doReflections)
  {
    model->plotOn(frameCand1, Components("bkgReflFunc2D"), LineWidth(2), LineStyle(2), LineColor(kViolet), Name("bkgReflFunc2D_curve"));
    bkgReflFunc2DCurve = (RooCurve *)frameCand1->findObject("bkgReflFunc2D_curve", RooCurve::Class());
    model->plotOn(frameCand1, Components("sgnReflFunc2D"), LineWidth(2), LineStyle(2), LineColor(kMagenta - 7), Name("sgnReflFunc2D_curve"));
    sgnReflFunc2DCurve = (RooCurve *)frameCand1->findObject("sgnReflFunc2D_curve", RooCurve::Class());
    model->plotOn(frameCand1, Components("reflBkgFunc2D"), LineWidth(2), LineStyle(2), LineColor(kCyan - 3), Name("reflBkgFunc2D_curve"));
    reflBkgFunc2DCurve = (RooCurve *)frameCand1->findObject("reflBkgFunc2D_curve", RooCurve::Class());
    model->plotOn(frameCand1, Components("reflSgnFunc2D"), LineWidth(2), LineStyle(2), LineColor(kMagenta - 5), Name("reflSgnFunc2D_curve"));
    reflSgnFunc2DCurve = (RooCurve *)frameCand1->findObject("reflSgnFunc2D_curve", RooCurve::Class());
    model->plotOn(frameCand1, Components("reflReflFunc2D"), LineWidth(2), LineStyle(2), LineColor(kAzure + 3), Name("reflReflFunc2D_curve"));
    reflReflFunc2DCurve = (RooCurve *)frameCand1->findObject("reflReflFunc2D_curve", RooCurve::Class());

    model->plotOn(frameCand2, Components("bkgReflFunc2D"), LineWidth(2), LineStyle(2), LineColor(kViolet), Name("bkgReflFunc2D_curve"));
    model->plotOn(frameCand2, Components("sgnReflFunc2D"), LineWidth(2), LineStyle(2), LineColor(kMagenta - 7), Name("sgnReflFunc2D_curve"));
    model->plotOn(frameCand2, Components("reflBkgFunc2D"), LineWidth(2), LineStyle(2), LineColor(kCyan - 3), Name("reflBkgFunc2D_curve"));
    model->plotOn(frameCand2, Components("reflSgnFunc2D"), LineWidth(2), LineStyle(2), LineColor(kMagenta - 5), Name("reflSgnFunc2D_curve"));
    model->plotOn(frameCand2, Components("reflReflFunc2D"), LineWidth(2), LineStyle(2), LineColor(kAzure + 3), Name("reflReflFunc2D_curve"));
  }

  // Plot the components of the model and keep references to the returned RooCurve objects
  model->plotOn(frameCand1, Components("bkgBkgFunc2D"), LineWidth(2), LineStyle(kDashed), LineColor(kRed), Name("bkgBkgFunc2D_curve"));
  RooCurve *bkgBkgFunc2DCurve = (RooCurve *)frameCand1->findObject("bkgBkgFunc2D_curve", RooCurve::Class());
  model->plotOn(frameCand1, Components("sgnBkgFunc2D"), LineStyle(2), LineWidth(2), LineColor(kOrange + 6), Name("sgnBkgFunc2D_curve"));
  RooCurve *sgnBkgFunc2DCurve = (RooCurve *)frameCand1->findObject("sgnBkgFunc2D_curve", RooCurve::Class());
  model->plotOn(frameCand1, Components("bkgSgnFunc2D"), LineStyle(2), LineWidth(2), LineColor(kOrange), Name("bkgSgnFunc2D_curve"));
  RooCurve *bkgSgnFunc2DCurve = (RooCurve *)frameCand1->findObject("bkgSgnFunc2D_curve", RooCurve::Class());
  model->plotOn(frameCand1, Components("sgnSgnFunc2D"), LineWidth(2), LineStyle(1), LineColor(kGreen), Name("sgnSgnFunc2D_curve"));
  RooCurve *sgnSgnFunc2DCurve = (RooCurve *)frameCand1->findObject("sgnSgnFunc2D_curve", RooCurve::Class());

  // Now cand 2
  model->plotOn(frameCand2, Components("bkgBkgFunc2D"), LineWidth(2), LineStyle(kDashed), LineColor(kRed), Name("bkgBkgFunc2D_curve"));
  model->plotOn(frameCand2, Components("sgnBkgFunc2D"), LineStyle(2), LineWidth(2), LineColor(kOrange + 6), Name("sgnBkgFunc2D_curve"));
  model->plotOn(frameCand2, Components("bkgSgnFunc2D"), LineStyle(2), LineWidth(2), LineColor(kOrange), Name("bkgSgnFunc2D_curve"));
  model->plotOn(frameCand2, Components("sgnSgnFunc2D"), LineWidth(2), LineStyle(1), LineColor(kGreen), Name("sgnSgnFunc2D_curve"));

  model->plotOn(frameCand1, Name("model_curve")); // Plot the full model
  RooCurve *modelMCCurve = (RooCurve *)frameCand1->findObject("model_curve", RooCurve::Class());
  model->plotOn(frameCand2, Name("model_curve")); // Plot the full model

  // Add entries to the legend
  legend->AddEntry("data_cand1", "Data", "ep");              // Entry for data
  legend->AddEntry(sgnSgnFunc2DCurve, "Sig 1 - Sig 2", "l"); // Entry for signal
  legend->AddEntry(bkgBkgFunc2DCurve, "Bkg 1 - Bkg 2", "l"); // Entry for background
  legend->AddEntry(sgnBkgFunc2DCurve, "Sig 1 - Bkg 2", "l"); // Entry for cross term sgn x bkg
  legend->AddEntry(bkgSgnFunc2DCurve, "Bkg 1 - Sgn 2", "l"); // Entry for cross term bkg x sgn

  if (doReflections)
  {
    legend_2->AddEntry(sgnReflFunc2DCurve, "Sig 1 - Refl 2", "l");   // Entry for cross term sgn x refl
    legend_2->AddEntry(bkgReflFunc2DCurve, "Bkg 1 - Refl 2", "l");   // Entry for cross term bkg x refl
    legend_2->AddEntry(reflSgnFunc2DCurve, "Refl 1 - Sig 2", "l");   // Entry for cross term refl x sgn
    legend_2->AddEntry(reflBkgFunc2DCurve, "Refl 1 - Bkg 2", "l");   // Entry for cross term refl x bkg
    legend_2->AddEntry(reflReflFunc2DCurve, "Refl 1 - Refl 2", "l"); // Entry for cross term refl x refl
  }

  // Create a canvas to draw the projections
  TCanvas *c_proj = new TCanvas(saveName, saveName, 1250, 600);
  c_proj->Divide(2, 1);     // Divide canvas into two pads
  gPad->SetLeftMargin(0.2); // Increase the left margin (default is ~0.1)

  // Draw the X projection
  c_proj->cd(1); // Go to the first pad
  frameCand1->SetXTitle("M(K#pi) (GeV/#it{c}^{2})");
  frameCand1->SetYTitle("Counts per 5 MeV/#it{c}^{2}");

  frameCand1->GetXaxis()->SetTitleSize(0.05); // Set X-axis title size
  frameCand1->GetXaxis()->SetTitleSize(0.05); // Set X-axis title size

  frameCand1->GetYaxis()->SetTitleSize(0.05); // Set Y-axis title size
  frameCand1->GetYaxis()->SetTitleSize(0.05); // Set Y-axis title size
  // frameCand1->GetYaxis()->SetLabelSize(0.05); // Set Y-axis label size
  // frameCand1->GetYaxis()->SetLabelSize(0.05); // Set Y-axis label size
  frameCand1->Draw();
  // Draw the legend on the canvas
  legend->Draw();

  // Draw the Y projection
  c_proj->cd(2); // Go to the second pad
  frameCand2->SetYTitle("Counts per 5 MeV/#it{c}^{2}");
  frameCand2->SetXTitle("M(K#pi) (GeV/#it{c}^{2})");

  frameCand2->GetXaxis()->SetTitleSize(0.05); // Set X-axis title size
  frameCand2->GetXaxis()->SetTitleSize(0.05); // Set X-axis title size

  frameCand2->GetYaxis()->SetTitleSize(0.05); // Set Y-axis title size
  frameCand2->GetYaxis()->SetTitleSize(0.05); // Set Y-axis title size
  // frameCand2->GetYaxis()->SetLabelSize(0.05); // Set Y-axis label size
  // frameCand2->GetYaxis()->SetLabelSize(0.05); // Set Y-axis label size
  frameCand2->Draw();
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
    double y1 = frameCand1->getHist("data_cand1")->Eval(x);
    double y2 = frameCand2->getHist("data_cand2")->Eval(x);

    histoCand1->SetBinContent(i, y1);
    histoCand2->SetBinContent(i, y2);

    // Optionally set errors, assuming you have uncertainties
    histoCand1->SetBinError(i, std::sqrt(y1)); // for Poisson stats
    histoCand2->SetBinError(i, std::sqrt(y2));
  }

  histoCand1->Divide(histoCand2);
  TLine *line = new TLine(1.7, 1, 2.05, 1);

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

// Function to analyse and compare the pT and y distributions of my candidates
// in a sideband region and the signal region after subtracting the background
void InvMassFitter2D::analyseKinematicDistributions(TFile *fout)
{
  RooDataSet *dataset = dynamic_cast<RooDataSet *>(_workspace.data("data"));
  if (!dataset)
    {
        cerr << "ERROR: dataset not found!" << endl;
        return;
    }
  // Get inv. mass, pt and y of the candidates
  RooRealVar *massCand1 = dynamic_cast<RooRealVar *>(dataset->get()->find("fMCand1"));
  RooRealVar *massCand2 = dynamic_cast<RooRealVar *>(dataset->get()->find("fMCand2"));
  RooRealVar *ptCand1 = dynamic_cast<RooRealVar *>(dataset->get()->find("fPtCand1"));
  RooRealVar *ptCand2 = dynamic_cast<RooRealVar *>(dataset->get()->find("fPtCand2"));
  RooRealVar *yCand1 = dynamic_cast<RooRealVar *>(dataset->get()->find("fYCand1"));
  RooRealVar *yCand2 = dynamic_cast<RooRealVar *>(dataset->get()->find("fYCand2"));
  // TODO: add phi (only available when I rerun the dataset)

  // Create formula variables for deltaPt and deltaY
  RooRealVar deltaPt("deltaPt", "ptCand1 - ptCand2", 0., -24.0, 24.0);
  RooRealVar deltaY("deltaY", "yCand1 - yCand2", 0., -2.0, 2.0);

  // Create a set of the variables for the new dataset
  RooArgSet vars(*ptCand1, *ptCand2, *yCand1, *yCand2, deltaPt, deltaY);

  RooRealVar *mean = _workspace.var("mean");
  RooRealVar *sigma = _workspace.var("sigma");
  RooRealVar *meanCand2 = _workspace.var("meanCand2");
  RooRealVar *sigmaCand2 = _workspace.var("sigmaCand2");

  // Create dataset in sideband region
  RooDataSet *datasetSidebandRegion = new RooDataSet("datasetSidebandRegion", "datasetSidebandRegion", vars);

  cout << "bkgLeftCand1: [" << _massMin << ", " << mean->getVal() - 4 * sigma->getVal() << "]" << endl;
  cout << "bkgRightCand1: [" << mean->getVal() + 4 * sigma->getVal() << ", " << _massMax << "]" << endl;

  for (int i = 0; i < dataset->numEntries(); ++i)
  {
    const RooArgSet *row = dataset->get(i);
    // Select sideband region
    if ((massCand1->getVal() > (mean->getVal() - 4 * sigma->getVal()) && massCand1->getVal() < (mean->getVal() + 4 * sigma->getVal())) ||
        (massCand2->getVal() > (meanCand2->getVal() - 4 * sigmaCand2->getVal()) && massCand2->getVal() < (meanCand2->getVal() + 4 * sigmaCand2->getVal()))) {
      continue;
    }

    ptCand1->setVal(row->getRealValue("fPtCand1"));
    ptCand2->setVal(row->getRealValue("fPtCand2"));
    yCand1->setVal(row->getRealValue("fYCand1"));
    yCand2->setVal(row->getRealValue("fYCand2"));

    // Calculate deltaPt and deltaY
    double deltaPtValue = ptCand1->getVal() - ptCand2->getVal();
    double deltaYValue = yCand1->getVal() - yCand2->getVal();

    // Set the delta values manually in the new dataset
    deltaPt.setVal(deltaPtValue);
    deltaY.setVal(deltaYValue);

    datasetSidebandRegion->add(vars);
  }

  // Create dataset in signal region
  // Manually fill the dataset with deltaPt and deltaY values
  // TODO: remove background candidates from signal region

  RooAbsPdf *bkgPdfCand1;
  RooAbsPdf *bkgPdfCand2;
  RooAbsPdf *sgnPdfCand1;
  RooAbsPdf *sgnPdfCand2;

  // Select signal function
  if (_sgnFuncOption == "gaus")
  {
    cout << "Signal function chosen: GAUSSIAN" << endl;
    sgnPdfCand1 = _workspace.pdf("sgnFuncGausCand1");
    sgnPdfCand2 = _workspace.pdf("sgnFuncGausCand2");
  }
  else if (_sgnFuncOption == "CB")
  {
    cout << "Signal function chosen: BIFURCATED CRYSTAL BALL" << endl;
    sgnPdfCand1 = _workspace.pdf("sgnFuncCBCand1");
    sgnPdfCand2 = _workspace.pdf("sgnFuncCBCand2");
  }
  else
  {
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
    bkgPdfCand1 = _workspace.pdf("bkgFuncExpoCand1");
    bkgPdfCand2 = _workspace.pdf("bkgFuncExpoCand2");
  }
  else if (_bkgFuncOption == "poly1")
  {
    cout << "Background function chosen: POLY 1" << endl;
    bkgPdfCand1 = _workspace.pdf("bkgFuncPoly1Cand1");
    bkgPdfCand2 = _workspace.pdf("bkgFuncPoly1Cand2");
  }
  else if (_bkgFuncOption == "poly2")
  {
    cout << "Background function chosen: POLY 2" << endl;
    bkgPdfCand1 = _workspace.pdf("bkgFuncPoly2Cand1");
    bkgPdfCand2 = _workspace.pdf("bkgFuncPoly2Cand2");
  }
  else if (_bkgFuncOption == "poly3")
  {
    cout << "Background function chosen: POLY 3" << endl;
    bkgPdfCand1 = _workspace.pdf("bkgFuncPoly3Cand1");
    bkgPdfCand2 = _workspace.pdf("bkgFuncPoly3Cand2");
  }
  else
  {
    cerr << "ERROR: background function not supported! \n Available options: expo, poly2, poly3. \n Exit!" << endl;
    return;
  }

  if (!bkgPdfCand1 || !bkgPdfCand2)
  {
    cerr << "ERROR: bkgPdf function not found!" << endl;
    return;
  }

  RooArgSet normSet(*massCand1, *massCand2);

  // Define the number of signal and background events for candidate 1
  RooRealVar nSgnCand1("nSgnCand1", "Number of signal events of candidate 1", 10000, 0, 100000);
  RooRealVar nBkgCand1("nBkgCand1", "Number of background events of candidate 1", 50000, 0, 200000);

  // Create the composite model for candidate 1
  RooAddPdf modelCand1("modelCand1", "Signal + Background of candidate 1", RooArgList(*sgnPdfCand1, *bkgPdfCand1), RooArgList(nSgnCand1, nBkgCand1));

  // Fit the model to the data for candidate 1
  RooFitResult* fitResultCand1 = modelCand1.fitTo(*dataset, Save());
  fitResultCand1->Print("v");
  if (fitResultCand1->status() != 0) {
    cout << "Fit did not converge for candidate 1!" << endl;
  }

  // Plot the fit results for candidate 1
  plotFitResults(dataset, massCand1, &modelCand1, fitResultCand1, "Fit Result for Candidate 1", "fit_result_cand1");

  // Define the number of signal and background events for candidate 2
  RooRealVar nSgnCand2("nSgnCand2", "Number of signal events of candidate 2", 10000, 0, 100000);
  RooRealVar nBkgCand2("nBkgCand2", "Number of background events of candidate 2", 50000, 0, 200000);

  // Create the composite model for candidate 2
  RooAddPdf modelCand2("modelCand2", "Signal + Background of candidate 2", RooArgList(*sgnPdfCand2, *bkgPdfCand2), RooArgList(nSgnCand2, nBkgCand2));

  // Fit the model to the data for candidate 2
  RooFitResult* fitResultCand2 = modelCand2.fitTo(*dataset, Save());
  fitResultCand2->Print("v");

  // Plot the fit results for candidate 2
  plotFitResults(dataset, massCand2, &modelCand2, fitResultCand2, "Fit Result for Candidate 2", "fit_result_cand2");

  // Create a flag variable
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
        cout << signalLikelihoodCand2 << " bkg likelihood: " << backgroundLikelihoodCand2 << endl;

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
  RooDataSet* signalData = (RooDataSet*)datasetWithFlag.reduce(Cut("signalFlagCand1 == signalFlagCand1::Signal"));

  RooDataSet *datasetSignalRegion = new RooDataSet("datasetSignalRegion", "datasetSignalRegion", vars);
  for (int i = 0; i < dataset->numEntries(); ++i)
  {
    const RooArgSet *row = dataset->get(i);
    // Select signal region
    if (massCand1->getVal() < (mean->getVal() - 3 * sigma->getVal()) || massCand1->getVal() > (mean->getVal() + 3 * sigma->getVal())) {
      continue;
    }
    if (massCand2->getVal() < (meanCand2->getVal() - 3 * sigmaCand2->getVal()) || massCand2->getVal() > (meanCand2->getVal() + 3 * sigmaCand2->getVal())) {
      continue;
    }

    ptCand1->setVal(row->getRealValue("fPtCand1"));
    ptCand2->setVal(row->getRealValue("fPtCand2"));
    yCand1->setVal(row->getRealValue("fYCand1"));
    yCand2->setVal(row->getRealValue("fYCand2"));

    // Calculate deltaPt and deltaY
    double deltaPtValue = ptCand1->getVal() - ptCand2->getVal();
    double deltaYValue = yCand1->getVal() - yCand2->getVal();

    // Set the delta values manually in the new dataset
    deltaPt.setVal(deltaPtValue);
    deltaY.setVal(deltaYValue);

    datasetSignalRegion->add(vars);
  }

  double binWidth = 1.; // Desired bin width
  int nBins = (int)((_ptMax - _ptMin) / binWidth);

  /// --------------------------------------------------------
  /// ------------------ PT DISTRIBUTIONS --------------------
  /// --------------------------------------------------------

  gROOT->SetBatch(kTRUE);

  double nSideband = datasetSidebandRegion->sumEntries();
  double nSignal = datasetSignalRegion->sumEntries();

  // Create histos with pt distributions
  TH1F* histSidebandCand1 = (TH1F*)datasetSidebandRegion->createHistogram("histSidebandCand1", *ptCand1, Binning(nBins, _ptMin, _ptMax));
  TH1F* histSignalCand1 = (TH1F*)datasetSignalRegion->createHistogram("histSignalCand1", *ptCand1, Binning(nBins, _ptMin, _ptMax));

  TH1F* histSidebandCand2 = (TH1F*)datasetSidebandRegion->createHistogram("histSidebandCand2", *ptCand2, Binning(nBins, _ptMin, _ptMax));
  TH1F* histSignalCand2 = (TH1F*)datasetSignalRegion->createHistogram("histSignalCand2", *ptCand2, Binning(nBins, _ptMin, _ptMax));
  // Normalise
  histSidebandCand1->Scale(1.0 / nSideband);
  histSignalCand1->Scale(1.0 / nSignal);
  histSidebandCand2->Scale(1.0 / nSideband);
  histSignalCand2->Scale(1.0 / nSignal);

  setHistoSignalSidebandStyle(histSidebandCand1, histSignalCand1, 1, "#it{p}_{T} (GeV/#it{c})");
  setHistoSignalSidebandStyle(histSidebandCand2, histSignalCand2, 2, "#it{p}_{T} (GeV/#it{c})");

  TLegend *legendPt = new TLegend(0.6, 0.7, 0.9, 0.9);
  legendPt->AddEntry(histSidebandCand1, "Sideband Region", "pl");
  legendPt->AddEntry(histSignalCand1, "Signal Region (not subtracted)", "pl");

  TCanvas* canvasPt = new TCanvas("canvasPt", "canvasPt", 1400, 650);
  canvasPt->Divide(2, 1); // Divide the canvas into 2 columns and 1 row
  canvasPt->cd(1); // Select the first pad
  histSidebandCand1->Draw("PE");
  histSignalCand1->Draw("samePE");
  legendPt->Draw();
  canvasPt->cd(2); // Select the second pad
  histSidebandCand1->Draw("PE");
  histSignalCand1->Draw("samePE");

  canvasPt->SaveAs("pt_sidebandRegion.png");
  fout->cd();
  canvasPt->Write();
  canvasPt->Close();

  /// --------------------------------------------------------
  /// ------------------- Y DISTRIBUTIONS --------------------
  /// --------------------------------------------------------
  gROOT->SetBatch(kTRUE);

  // Create histos with pt distributions
  TH1F* histSidebandYCand1 = (TH1F*)datasetSidebandRegion->createHistogram("histSidebandYCand1", *yCand1, Binning(20, -0.6, 0.6));
  TH1F* histSignalYCand1 = (TH1F*)datasetSignalRegion->createHistogram("histSignalYCand1", *yCand1, Binning(20, -0.6, 0.6));

  TH1F* histSidebandYCand2 = (TH1F*)datasetSidebandRegion->createHistogram("histSidebandYCand2", *yCand2, Binning(20, -0.6, 0.6));
  TH1F* histSignalYCand2 = (TH1F*)datasetSignalRegion->createHistogram("histSignalYCand2", *yCand2, Binning(20, -0.6, 0.6));
  // Normalise
  histSidebandYCand1->Scale(1.0 / nSideband);
  histSignalYCand1->Scale(1.0 / nSignal);
  histSidebandYCand2->Scale(1.0 / nSideband);
  histSignalYCand2->Scale(1.0 / nSignal);

  setHistoSignalSidebandStyle(histSidebandYCand1, histSignalYCand1, 1, "#it{y}");
  setHistoSignalSidebandStyle(histSidebandYCand2, histSignalYCand2, 2, "#it{y}");

  TLegend *legendY = new TLegend(0.3, 0.3, 0.7, 0.5);
  legendY->AddEntry(histSidebandYCand1, "Sideband Region", "pl");
  legendY->AddEntry(histSignalYCand1, "Signal Region (not subtracted)", "pl");

  TCanvas* canvasY = new TCanvas("canvasY", "canvasY", 1400, 650);
  canvasY->Divide(2, 1); // Divide the canvas into 2 columns and 1 row
  canvasY->cd(1); // Select the first pad
  histSidebandYCand1->Draw("PE");
  histSignalYCand1->Draw("samePE");
  legendY->Draw();
  canvasY->cd(2); // Select the second pad
  histSidebandYCand2->Draw("PE");
  histSignalYCand2->Draw("samePE");

  canvasY->SaveAs("y_sidebandRegion.png");
  fout->cd();
  canvasY->Write();
  canvasY->Close();

  /// --------------------------------------------------------------
  /// ------------------ DELTA Y, PT DISTRIBUTIONS --------------------
  /// --------------------------------------------------------------

  gROOT->SetBatch(kTRUE);

  // Create histos with pt distributions
  TH1F* histSidebandDeltaPt = (TH1F*)datasetSidebandRegion->createHistogram("histSidebandDeltaPt", deltaPt, Binning(40, -_ptMax, _ptMax));
  TH1F* histSidebandDeltaY = (TH1F*)datasetSidebandRegion->createHistogram("histSidebandDeltaY", deltaY, Binning(40, -2.0, 2.0));

  TH1F* histSignalDeltaPt = (TH1F*)datasetSignalRegion->createHistogram("histSignalDeltaPt", deltaPt, Binning(40, -_ptMax, _ptMax));
  TH1F* histSignalDeltaY = (TH1F*)datasetSignalRegion->createHistogram("histSignalDeltaY", deltaY, Binning(40, -2.0, 2.0));
  // Normalise
  histSidebandDeltaPt->Scale(1.0 / nSideband);
  histSignalDeltaPt->Scale(1.0 / nSignal);
  histSidebandDeltaY->Scale(1.0 / nSideband);
  histSignalDeltaY->Scale(1.0 / nSignal);

  setHistoSignalSidebandStyle(histSidebandDeltaPt, histSignalDeltaPt, "#Delta#it{p}_{T} (GeV/#it{c})");
  setHistoSignalSidebandStyle(histSidebandDeltaY, histSignalDeltaY, "#Delta#it{y}");

  TLegend *legendDelta = new TLegend(0.3, 0.6, 0.7, 0.8);
  legendDelta->AddEntry(histSidebandDeltaPt, "Sideband Region", "pl");
  legendDelta->AddEntry(histSignalDeltaPt, "Signal Region (not subtracted)", "pl");

  TCanvas* canvasDelta = new TCanvas("canvasDelta", "canvasDelta", 1400, 650);
  canvasDelta->Divide(2, 1); // Divide the canvas into 2 columns and 1 row
  canvasDelta->cd(1); // Select the first pad
  histSidebandDeltaPt->Draw("PE");
  histSignalDeltaPt->Draw("samePE");
  legendDelta->Draw();
  canvasDelta->cd(2); // Select the second pad
  histSidebandDeltaY->Draw("PE");
  histSignalDeltaY->Draw("samePE");

  canvasDelta->SaveAs("delta_sidebandRegion.png");
  fout->cd();
  canvasDelta->Write();
  canvasDelta->Close();

}

void InvMassFitter2D::plotFitResults(RooDataSet* dataset, RooRealVar* mass, RooAbsPdf* model, RooFitResult* fitResult, const char* title, const char* canvasName) {
  // Create a RooPlot for the mass variable
  RooPlot* frame = mass->frame(Title(title));

  // Plot the data on the frame
  dataset->plotOn(frame, Name("data"));

  // Plot the fitted model on the frame
  model->plotOn(frame, Name("model"), LineColor(kBlue));

  // Plot the components of the model (signal and background)
  model->plotOn(frame, Components("sgnPdfCand1"), LineStyle(kDashed), LineColor(kRed), Name("signal"));
  model->plotOn(frame, Components("bkgPdfCand1"), LineStyle(kDashed), LineColor(kGreen), Name("background"));

  // Customize the plot
  frame->SetTitle(title);
  frame->GetXaxis()->SetTitle("Mass");
  frame->GetYaxis()->SetTitle("Events");

  // Create a legend
  TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
  legend->AddEntry("data", "Data", "lep");
  legend->AddEntry("model", "Fit Model", "l");
  legend->AddEntry("signal", "Signal Component", "l");
  legend->AddEntry("background", "Background Component", "l");

  // Draw the plot on a TCanvas
  TCanvas* canvas = new TCanvas(canvasName, title, 800, 600);
  frame->Draw();
  legend->Draw();
  canvas->Update();

  // Optionally, save the plot to a file
  canvas->SaveAs(Form("%s.png", canvasName));
}

void InvMassFitter2D::setHistoSignalSidebandStyle(TH1F *hSideband, TH1F *hSignal, int const& candNum, TString physVar) {
  hSideband->SetLineColor(kRed);
  hSideband->SetMarkerColor(kRed);
  hSideband->SetMarkerStyle(21);

  hSignal->SetLineColor(kBlue);
  hSignal->SetMarkerColor(kBlue);
  hSignal->SetMarkerStyle(22);

  hSideband->GetYaxis()->SetRangeUser(0, 0.5);
  hSignal->GetYaxis()->SetRangeUser(0, 0.5);

  TString title = Form("%s of candidate %d", physVar.Data(), candNum);

  hSideband->SetTitle(title); // Set the plot title
  hSideband->GetXaxis()->SetTitle(physVar);   // Set the x-axis label
  hSideband->GetYaxis()->SetTitle("Normalised counts");   // Set the y-axis label

  hSignal->SetTitle(title); // Set the plot title
  hSignal->GetXaxis()->SetTitle(physVar);   // Set the x-axis label
  hSignal->GetYaxis()->SetTitle("Normalised counts");   // Set the y-axis label
}

void InvMassFitter2D::setHistoSignalSidebandStyle(TH1F *hSideband, TH1F *hSignal, TString physVar) {
  hSideband->SetLineColor(kRed);
  hSideband->SetMarkerColor(kRed);
  hSideband->SetMarkerStyle(21);

  hSignal->SetLineColor(kBlue);
  hSignal->SetMarkerColor(kBlue);
  hSignal->SetMarkerStyle(22);

  hSideband->GetYaxis()->SetRangeUser(0, 0.2);
  hSignal->GetYaxis()->SetRangeUser(0, 0.2);

  TString title = Form("%s", physVar.Data());

  hSideband->SetTitle(title); // Set the plot title
  hSideband->GetXaxis()->SetTitle(physVar);   // Set the x-axis label
  hSideband->GetYaxis()->SetTitle("Normalised counts");   // Set the y-axis label

  hSignal->SetTitle(title); // Set the plot title
  hSignal->GetXaxis()->SetTitle(physVar);   // Set the x-axis label
  hSignal->GetYaxis()->SetTitle("Normalised counts");   // Set the y-axis label
}