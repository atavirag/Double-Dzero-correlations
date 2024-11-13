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
    RooArgSet vars(rooPtCand1, rooPtCand2, rooMCand1, rooMCand2, rooYCand1, rooYCand2); // Cand 1: D, Cand 2: Dbar (OS)

    RooRealVar weightCand1("weightCand1", "weights of cand 1", 1., 0., 100.);
    RooRealVar weightCand2("weightCand2", "weights of cand 2", 1., 0., 100.);
    RooArgSet weightedVars(rooPtCand1, rooPtCand2, rooMCand1, rooMCand2, rooYCand1, rooYCand2, weightCand1, weightCand2);
    RooFormulaVar combinedWeight("combinedWeight", "combined weight", "weightCand1 * weightCand2", RooArgList(weightCand1, weightCand2));

    // Create an empty dataset with the variables and category
    RooDataSet data("data", "data", vars);
    RooDataSet weightedData("weightedData", "weightedData", weightedVars, RooFit::WeightVar(combinedWeight.GetName()));

    fillDataset(data, vars);
    fillDataset(weightedData, weightedVars);

    checkCorrelations(data);
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
    /// | ------------------------------------------------------------------ |
    /// | ------------------- BACKGROUND FUNCTIONS ------------------------- |
    /// | ------------------------------------------------------------------ |
    cout << "loading functions" << endl;
    // bkg expo
    RooRealVar tauCand1("tauCand1", "tauCand1", -1, -10., 6.);
    RooRealVar tauCand2("tauCand2", "tauCand2", -1, -10., 6.);

    RooAbsPdf *bkgFuncExpoCand1 = new RooExponential("bkgFuncExpoCand1", "background fit function of candidate 1", *massCand1, tauCand1);
    RooAbsPdf *bkgFuncExpoCand2 = new RooExponential("bkgFuncExpoCand2", "background fit function of candidate 2", *massCand2, tauCand2);
    _workspace.import(*bkgFuncExpoCand1);
    _workspace.import(*bkgFuncExpoCand2);

    /// | ------------------------------------------------------------------ |
    /// | ----------------------- SIGNAL FUNCTIONS ------------------------- |
    /// | ------------------------------------------------------------------ |
    // signal pdf
    RooRealVar mean("mean", "mean for signal fit", 1.85, 1.83, 1.9);
    RooRealVar sigma("sigma", "sigma for signal", 0.02, 0.01, 0.07);

    RooAbsPdf *sgnFuncGausCand1 = new RooGaussian("sgnFuncGausCand1", "signal pdf of candidate 1", *massCand1, mean, sigma);
    RooAbsPdf *sgnFuncGausCand2 = new RooGaussian("sgnFuncGausCand2", "signal pdf of candidate 2", *massCand2, mean, sigma);
    _workspace.import(*sgnFuncGausCand1);
    _workspace.import(*sgnFuncGausCand2);

    /// | ------------------------------------------------------------------ |
    /// | ------------------- REFLECTION FUNCTIONS ------------------------- |
    /// | ------------------------------------------------------------------ |
    // reflection Gaussian
    RooRealVar meanRefl("meanRefl", "mean for reflections", 1.85, 0.0, 2.15);
    RooRealVar sigmaRefl("sigmaRefl", "sigma for reflection", 0.012, 0, 0.3);
    RooRealVar meanReflDoubleGaus("meanReflDoubleGaus", "mean for reflection double gaussian", 1.85, 0.0, 1.90);
    RooRealVar sigmaReflDoubleGaus("sigmaReflDoubleGaus", "sigmaReflDoubleGaus", 0.012, 0.0, 0.2);

    RooGaussian gausRefl1Cand1("gausRefl1Cand1", "gausRefl1Cand1", *massCand1, meanRefl, sigmaRefl);
    RooGaussian gausRefl2Cand1("gausRefl2Cand1", "gausRefl2Cand1", *massCand1, meanReflDoubleGaus, sigmaReflDoubleGaus);
    RooGaussian gausRefl1Cand2("gausRefl1Cand2", "gausRefl1Cand2", *massCand2, meanRefl, sigmaRefl);
    RooGaussian gausRefl2Cand2("gausRefl2Cand2", "gausRefl2Cand2", *massCand2, meanReflDoubleGaus, sigmaReflDoubleGaus);

    RooRealVar fracRefl("fracRefl", "frac of two gauss of candidate 1", 0.5, 0, 1.);
    RooAbsPdf *reflFuncDoubleGausCand1 = new RooAddPdf("reflFuncDoubleGausCand1", "reflection pdf of candidate 1", RooArgList(gausRefl1Cand1, gausRefl2Cand1), fracRefl);
    RooAbsPdf *reflFuncDoubleGausCand2 = new RooAddPdf("reflFuncDoubleGausCand2", "reflection pdf of candidate 2", RooArgList(gausRefl1Cand2, gausRefl2Cand2), fracRefl);
    _workspace.import(*reflFuncDoubleGausCand1);
    _workspace.import(*reflFuncDoubleGausCand2);

    cout << "Workspace filled with functions" << endl;
}

void InvMassFitter2D::do2DFit(Bool_t draw, Bool_t doReflections, TFile *fout)
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

    cout << "\n\n                  Integrated efficiency: " << _integratedEfficiency << endl;

    RooAbsPdf *bkgPdfCand1 = _workspace.pdf("bkgFuncExpoCand1");
    RooAbsPdf *sgnPdfCand1 = _workspace.pdf("sgnFuncGausCand1");
    RooAbsPdf *bkgPdfCand2 = _workspace.pdf("bkgFuncExpoCand2");
    RooAbsPdf *sgnPdfCand2 = _workspace.pdf("sgnFuncGausCand2");

    if (!bkgPdfCand1 || !sgnPdfCand1 || !bkgPdfCand2 || !sgnPdfCand2)
    {
        cerr << "ERROR: attempted to access an empty pointer!" << endl;
        return;
    }

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

    mean->setConstant(kTRUE);
    sigma->setConstant(kTRUE);

    cout << "Mean, sigma and tau values set" << endl;

    // Por ahora, asumo que tengo una doble gaussiana
    RooAbsPdf *reflPdfCand1 = _workspace.pdf("reflFuncDoubleGausCand1");
    RooAbsPdf *reflPdfCand2 = _workspace.pdf("reflFuncDoubleGausCand2");

    if (!reflPdfCand1 || !reflPdfCand2)
    {
        cerr << "ERROR: attempted to access an empty pointer (reflPdf functions)!" << endl;
        return;
    }

    if (doReflections)
    {
        RooRealVar *sigmaRefl = _workspace.var("sigmaRefl");
        RooRealVar *meanRefl = _workspace.var("meanRefl");

        RooRealVar *sigmaReflDoubleGaus = _workspace.var("sigmaReflDoubleGaus");
        RooRealVar *meanReflDoubleGaus = _workspace.var("meanReflDoubleGaus");

        RooRealVar *fracRefl = _workspace.var("fracRefl");

        meanRefl->setVal(_meanRefl->getVal());
        meanRefl->setError(_meanRefl->getError());
        sigmaRefl->setVal(_sigmaRefl->getVal());
        sigmaRefl->setError(_sigmaRefl->getError());

        cout << "sigmaRefl: " << sigmaRefl->getVal() << endl;
        cout << "_sigmaRefl: " << _sigmaRefl->getVal() << endl;

        meanRefl->setConstant(kTRUE);
        sigmaRefl->setConstant(kTRUE);

        meanReflDoubleGaus->setVal(_meanReflDoubleGaus->getVal());
        meanReflDoubleGaus->setError(_meanReflDoubleGaus->getError());
        sigmaReflDoubleGaus->setVal(_sigmaReflDoubleGaus->getVal());
        sigmaReflDoubleGaus->setError(_sigmaReflDoubleGaus->getError());

        fracRefl->setVal(_fracRefl->getVal());
        fracRefl->setConstant(kTRUE);

        meanReflDoubleGaus->setConstant(kTRUE);
        sigmaReflDoubleGaus->setConstant(kTRUE);
    }

    RooArgSet normSet(*massCand1, *massCand2);

    // RAW PARAMETERS
    RooRealVar *nBkg1 = new RooRealVar("nBkg1", "background yield of cand 1", 100.0, 0.0, 200.0);
    // RooRealVar *nSgn1 = new RooRealVar("nSgn1", "signal yield of cand1", 0.3 * _rawYield->getVal(), 0.0, 1.2 * numEntries);
    RooRealVar *nSgn1 = new RooRealVar("nSgn1", "signal yield of cand1", 30.0, 0.0, 300.0);
    RooFormulaVar *nRefl1 = new RooFormulaVar("nRefl1", "reflected signal yield of cand1", "@0 * @1", RooArgList(_reflOverSgn, *nSgn1));

    RooRealVar *nBkg2 = new RooRealVar("nBkg2", "background yield of cand2", 100.0, 0.0, 200.0);
    // RooRealVar *nSgn2 = new RooRealVar("nSgn2", "signal yield of cand2", 0.3 * _rawYield->getVal(), 0.0, 1.2 * numEntries);
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
    RooRealVar *nBkgCorr1 = new RooRealVar("nBkgCorr1", "Corrected background yield of cand 1", 900, 0.0, 2000.0);
    RooRealVar *nSgnCorr1 = new RooRealVar("nSgnCorr1", "Corrected signal yield of cand1", 100, 0.0, 1000.0);
    RooFormulaVar *nReflCorr1 = new RooFormulaVar("nReflCorr1", "Corrected reflected signal yield of cand1", "@0 * @1", RooArgList(_reflOverSgn, *nSgnCorr1));

    RooRealVar *nBkgCorr2 = new RooRealVar("nBkgCorr2", "Corrected background yield of cand2", 900, 0.0, 2000.0);
    RooRealVar *nSgnCorr2 = new RooRealVar("nSgnCorr2", "Corrected signal yield of cand2", 100, 0.0, 1000.0);
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
    RooFormulaVar diffSgnCorr("diffSgnCorr", "difference between nSgnCorr1 and nSgnCorr2", "@0 - @1", RooArgList(*nSgnCorr1, *nSgnCorr2));
    RooFormulaVar diffBkgCorr("diffBkgCorr", "difference between nBkgCorr1 and nBkgCorr2", "@0 - @1", RooArgList(*nBkgCorr1, *nBkgCorr2));

    // Define the Gaussian constraint to keep nSgn1 and nSgn2 symmetric
    RooRealVar sigmaConstraint("sigmaConstraint", "width of Gaussian constraint", 10.0);              // adjust this as needed
    RooRealVar sigmaConstraintCorr("sigmaConstraintCorr", "width of Gaussian constraint corr", 5.0); // adjust this as needed
    RooGaussian *constraintSgn = new RooGaussian("constraintSgn", "Gaussian constraint on nSgn1 - nSgn2", diffSgn, RooConst(0), sigmaConstraint);
    RooGaussian *constraintSgnCorr = new RooGaussian("constraintSgnCorr", "Gaussian constraint on nSgnCorr1 - nSgnCorr2", diffSgnCorr, RooConst(0), sigmaConstraintCorr);
    RooGaussian *constraintBkgCorr = new RooGaussian("constraintBkgCorr", "Gaussian constraint on nBkgCorr1 - nBkgCorr2", diffBkgCorr, RooConst(0), sigmaConstraintCorr);
    // RooArgSet constraintsCorr(*constraintSgnCorr, *constraintBkgCorr);

    // Final fit
    if (doReflections)
    {
        _totPdf2D = new RooAddPdf("_totPdf2D", "background + signal + reflection pdf 2D", RooArgList(sgnSgnFunc2D, sgnReflFunc2D, reflSgnFunc2D, sgnBkgFunc2D, bkgSgnFunc2D, bkgBkgFunc2D, bkgReflFunc2D, reflBkgFunc2D, reflReflFunc2D),
                                  RooArgList(nSgnSgn, nSgnRefl, nReflSgn, nSgnBkg, nBkgSgn, nBkgBkg, nBkgRefl, nReflBkg, nReflRefl));
        _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "Corrected background + signal + reflection pdf 2D", RooArgList(sgnSgnFunc2D, sgnReflFunc2D, reflSgnFunc2D, sgnBkgFunc2D, bkgSgnFunc2D, bkgBkgFunc2D, bkgReflFunc2D, reflBkgFunc2D, reflReflFunc2D),
                                          RooArgList(nSgnSgnCorr, nSgnReflCorr, nReflSgnCorr, nSgnBkgCorr, nBkgSgnCorr, nBkgBkgCorr, nBkgReflCorr, nReflBkgCorr, nReflReflCorr));
    }
    else
    {
        _totPdf2D = new RooAddPdf("_totPdf2D", "background + signal pdf 2D", RooArgList(bkgBkgFunc2D, sgnSgnFunc2D, sgnBkgFunc2D, bkgSgnFunc2D), RooArgList(nBkgBkg, nSgnSgn, nSgnBkg, nBkgSgn));
        _weightedTotPdf2D = new RooAddPdf("_weightedTotPdf2D", "background + signal pdf 2D", RooArgList(bkgBkgFunc2D, sgnSgnFunc2D, sgnBkgFunc2D, bkgSgnFunc2D), RooArgList(nBkgBkgCorr, nSgnSgnCorr, nSgnBkgCorr, nBkgSgnCorr));
    }

    _totPdf2D->fixCoefNormalization(normSet);
    _weightedTotPdf2D->fixCoefNormalization(normSet);

    // Include the constraint in the total model
    RooProdPdf *modelWithConstraint = new RooProdPdf("modelWithConstraint", "model with signal symmetry constraint", RooArgList(*_totPdf2D, *constraintSgn));
    RooProdPdf *modelWithConstraintCorr = new RooProdPdf("modelWithConstraintCorr", "model with signal symmetry constraint Corr", RooArgList(*_weightedTotPdf2D, *constraintSgnCorr));

    // Perform the fit
    RooFitResult *resultConstrained = modelWithConstraint->fitTo(*dataset, Save(), NormSet(normSet));

    // tauCand1->setConstant(kTRUE);
    // tauCand2->setConstant(kTRUE);

    RooFitResult *resultConstrainedCorr = modelWithConstraintCorr->fitTo(*weightedDataset, RooFit::SumW2Error(kTRUE), Save(), NormSet(normSet));

    if (resultConstrainedCorr->status() == 0)
    {
        std::cout << "Fit Converged!" << std::endl;
    }
    else
    {
        std::cout << "Fit did not converge." << std::endl;
    }

    resultConstrained->Print("v");
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
    TH2D *histFit = new TH2D("histFit", "Fit Results;m(D_{1}^{0}) (GeV/#it{c}^2);#it{M}(D_{2}^{0}) (GeV/#it{c}^2)",
                             45, _massMin, _massMax, 45, _massMin, _massMax);

    TH2D *histFitWeighted = new TH2D("histFitWeighted", "Efficiency corrected fit Results;m(D_{1}^{0}) (GeV/#it{c}^2);#it{M}(D_{2}^{0}) (GeV/#it{c}^2)",
                                     45, _massMin, _massMax, 45, _massMin, _massMax);

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

    TCanvas *cRatioFunctions = new TCanvas("cRatioFunctions", "cRatioFunctions", 1000, 700);
    ratioHist2D->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
    ratioHist2D->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
    ratioHist2D->SetTitle("Normalised Efficiency-corrected / non-corrected fit functions");
    ratioHist2D->Draw("colz");
    cRatioFunctions->Write();

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

    TCanvas *cRatioMassIntEff = new TCanvas("cRatioMassIntEff", "cRatioMassIntEff", 1000, 700);
    ratioHistMass2DIntEff->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
    ratioHistMass2DIntEff->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
    ratioHistMass2DIntEff->SetTitle("(raw yield/(integrated efficiency)^{2}) / efficiency-corrected distr.");
    ratioHistMass2DIntEff->Draw("colz");
    cRatioMassIntEff->Write();

    // -----------------------------------------------------------------------
    // -- Plot pT distributions before and after correcting by efficiencies --
    // -----------------------------------------------------------------------

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
        histFit->GetXaxis()->SetTitle("#it{M}(D^{0}_{1}) (GeV/#it{c}^{2})");
        histFit->GetYaxis()->SetTitle("#it{M}(D^{0}_{2}) (GeV/#it{c}^{2})");
        histFit->GetZaxis()->SetTitle("Counts per (3.5 MeV/#it{c}^{2})^{2}");
        histFit->GetXaxis()->CenterTitle();
        histFit->GetXaxis()->CenterTitle();
        histFit->GetXaxis()->SetDecimals(2);
        histFit->GetYaxis()->SetDecimals(2);
        hMassCorrelations->GetXaxis()->SetDecimals(2);
        hMassCorrelations->GetYaxis()->SetDecimals(2);

        hMassCorrelations->SetTitle("");

        TPaveText *ptTitle = new TPaveText(0.12, 0.85, 0.4, 0.94, "brNDC");
        ptTitle->SetTextFont(42);
        ptTitle->SetTextSize(0.04);
        ptTitle->SetBorderSize(0);
        ptTitle->SetFillStyle(0);
        ptTitle->SetTextAlign(11);

        // Add text lines
        ptTitle->AddText("#scale[1.35]{ALICE Performance}");
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

        hMassCorrelations->Draw("LEGO2");
        histFit->Draw("sameLEGO0");
        gStyle->SetOptStat(0);
        ptTitle->Draw("same");
        pt->Draw("same");
        fout->cd();
        c2d->Write();

        TString c1dName = cName + "_colz";
        TCanvas *c1d = new TCanvas(c1dName, c1dName, 1000, 700);
        hMassCorrelations->Draw("colz");
        c1d->Write();
    }
}

void InvMassFitter2D::plotProjectionsAfterFit(RooProdPdf *model, RooDataSet *dataset, TString saveName, TFile *fout, bool doReflections)
{

    RooRealVar *massCand1 = _workspace.var("fMCand1");
    RooRealVar *massCand2 = _workspace.var("fMCand2");

    massCand1->setBins(50);
    massCand2->setBins(50);

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
    TCanvas *c_proj = new TCanvas(saveName, saveName, 1200, 500);
    c_proj->Divide(2, 1); // Divide canvas into two pads

    // Draw the X projection
    c_proj->cd(1); // Go to the first pad
    frameCand1->Draw();
    // Draw the legend on the canvas
    legend->Draw();

    // Draw the Y projection
    c_proj->cd(2); // Go to the second pad
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
    c_proj->Close();

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
    double weight = 1.;
    const double epsilon = 1e-4; // Small value to avoid division by zero

    for (int ipt = 1; ipt <= nbinsPt; ipt++)
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

void InvMassFitter2D::checkCorrelations(RooDataSet &data) {
    RooRealVar *ptCand1 = _workspace.var("fPtCand1");
    RooRealVar *ptCand2 = _workspace.var("fPtCand2");
    RooRealVar *yCand1 = _workspace.var("fYCand1");
    RooRealVar *yCand2 = _workspace.var("fYCand2");

    // Create formula variables for deltaPt and deltaY
    RooRealVar deltaPt("deltaPt", "ptCand1 - ptCand2", 0., -24.0, 24.0);
    RooRealVar deltaY("deltaY", "yCand1 - yCand2", 0., -2.0, 2.0);

    // Create a set of the variables for the new dataset
    RooArgSet vars(*ptCand1, *ptCand2, *yCand1, *yCand2, deltaPt, deltaY);

    // Manually fill the dataset with deltaPt and deltaY values
    RooDataSet *newData = new RooDataSet("newData", "newData", vars);
    for (int i = 0; i < data.numEntries(); ++i) {
        const RooArgSet *row = data.get(i);
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

        newData->add(vars);
    }

    // Create frames for plotting the delta values
    TCanvas *cDeltaPt = new TCanvas("cDeltaPt", "Delta Pt", 800, 600);
    RooPlot *framePt = deltaPt.frame();
    newData->plotOn(framePt);  // This will plot deltaPt
    framePt->Draw();
    cDeltaPt->Write();

    TCanvas *cDeltaY = new TCanvas("cDeltaY", "Delta y", 800, 600);
    RooPlot *frameY = deltaY.frame();
    newData->plotOn(frameY);  // This will plot deltaY
    frameY->Draw();
    cDeltaY->Write();
}