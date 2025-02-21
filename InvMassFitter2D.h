// Implementation of InvMassFitter2D class
// Author: Andrea Tavira García, IJCLab (Orsay, France)

#ifndef INVMASSFITTER2D_H
#define INVMASSFITTER2D_H

#include "TTree.h"
#include "Math/Vector4D.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "TCanvas.h"
#include "TFile.h"

class InvMassFitter2D {
    public:
        // Constructor
        InvMassFitter2D();
        InvMassFitter2D(TTree* tree, const char *pairType); // pairType can be "OS" or "LS"

        // Destructor
        //~InvMassFitter2D();

        // Member functions
        void createDataset();
        void fillDataset(RooDataSet &data, RooArgSet &vars);
        void fillWorkspace(RooDataSet *dataset);
        // Functions to set parameters manually
        void set1DParameters(const RooArgSet *vars1D, double const &reflOverSgn, double const &integratedEfficiency);
        void setLumi(double const& lumi);
        void setPtLims(double const& ptMin, double const& ptMax);
        void setMassLims(double const& massMin, double const& massMax);
        void setSgnFunc(TString const& sgnFunc);
        void setBkgFunc(TString  const& bkgFunc);
        void setReflFunc(TString  const& reflFunc);
        void setEfficiencyMap(TH2F *h);
        // Functions for checks and calculations
        double calculateWeights(double const& y, double const& pt);
        void analyseKinematicDistributions(TFile *fout);
        //double calculateIntegratedEfficiency();
        ROOT::Math::PxPyPzMVector createLorentzVector(double const& phi, double const& y, double const& pt, double const& m);
        RooRealVar *getYieldInRange(RooFitResult *fitResult, RooRealVar *massCand1, RooRealVar *massCand2, RooProdPdf function, RooFormulaVar nCands, TString range);
        // Fitting and plotting functions
        void do2DFit(Bool_t draw, Bool_t doReflections, Bool_t isMc, TFile *fout);
        void plotProjectionsAfterFit(RooProdPdf *model, RooDataSet *dataset, TString saveName, TFile *fout, bool doReflections);
        void plot2DFit(TH2D *hMassCorrelations, TH2D* histFit, RooProdPdf *model, Bool_t draw, TFile *fout, TString const& cName);
        void plotFitResults(RooDataSet* dataset, RooRealVar* mass, RooAbsPdf* model, RooFitResult* fitResult, const char* title, const char* canvasName);
        void setHistoSignalSidebandStyle(TH1F *hSideband, TH1F *hSignal, int const& candNum, TString physVar);
        void setHistoSignalSidebandStyle(TH1F *hSideband, TH1F *hSignal, TString physVar);

    private:
        // Member variables
        TTree* _tree;
        const char *_pairType;

        int nentries;
        double _lumi = 0.;
        float _massMin = 1.70;
        float _massMax = 2.05;
        float _ptMin = 0.;
        float _ptMax = 50.;

        float ptCand1 = 0., ptCand2 = 0., yCand1 = 0., yCand2 = 0., phiCand1 = 0., phiCand2 = 0., mDCand1 = 0., mDCand2 = 0., mDbarCand1 = 0., mDbarCand2 = 0.;
        uint8_t typeCand1 = 0, typeCand2 = 0, typePair = 0;

        int _typeOfBkgPdf;
        int _typeOfSgnPdf;
        int _typeOfReflPdf;
        TString _sgnFuncOption;
        TString _bkgFuncOption;
        TString _reflFuncOption;

        RooRealVar* _tau;
        RooRealVar* _mean;
        RooRealVar* _sigma;
        RooRealVar* _meanRefl;
        RooRealVar* _sigmaRefl;
        RooRealVar* _fracRefl;
        RooRealVar* _meanReflDoubleGaus;
        RooRealVar* _sigmaReflDoubleGaus;
        RooRealVar* _rawYield;
        double _reflOverSgn = 0.;
        double _integratedEfficiency = 0.;
        TH2F *_efficiencyMap;

        RooRealVar rooPtCand1;
        RooRealVar rooPtCand2;
        RooRealVar rooMCand1;
        RooRealVar rooMCand2;
        RooRealVar rooYCand1;
        RooRealVar rooYCand2;
        RooRealVar rooPhiCand1;
        RooRealVar rooPhiCand2;

        RooWorkspace _workspace;
        void _loadTreeInfo();

        RooAddPdf* _bkgPdfCand1;
        RooAddPdf* _sgnPdfCand1;
        RooAddPdf* _reflPdfCand1;
        RooAddPdf* _totPdfCand1;

        RooAddPdf* _bkgPdfCand2;
        RooAddPdf* _sgnPdfCand2;
        RooAddPdf* _reflPdfCand2;
        RooAddPdf* _totPdfCand2;

        RooAddPdf* _totPdf2D;
        RooAddPdf* _weightedTotPdf2D;

};

#endif // INVMASSFITTER2D_H
