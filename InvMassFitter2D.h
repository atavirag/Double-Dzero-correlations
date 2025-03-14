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
#include "TH2D.h"
#include "RooBernstein.h"
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
        void removeAmbiguous(bool remove);
        void setPtLims(double const& ptMin, double const& ptMax);
        void setPtPairLims(double const& ptMinPair, double const& ptMaxPair);
        void setMassLims(double const& massMin, double const& massMax);
        void setSgnFunc(TString const& sgnFunc);
        void setBkgFunc(TString  const& bkgFunc);
        void setReflFunc(TString  const& reflFunc);
        void setEfficiencyMap(TH2F *h);
        // Functions for checks and calculations
        double calculateWeights(double const& y, double const& pt);
        void analyseKinematicDistributions(TFile *fout, RooDataSet *dataset, const char *suffix);
        void selectFitFunctions(RooAbsPdf* &sgnPdfCand1,RooAbsPdf* &sgnPdfCand2,RooAbsPdf* &bkgPdfCand1,
                                RooAbsPdf* &bkgPdfCand2, RooAbsPdf* &reflPdfCand1, RooAbsPdf* &reflPdfCand2, int8_t fitType);
        RooFitResult *fitAndPlot1DCandidate(RooAbsPdf* sgnPdf, RooAbsPdf* bkgPdf, RooRealVar& massVar, RooDataSet* dataset,
                                   const std::string& candidateName, const std::string& plotFilename);
        //double calculateIntegratedEfficiency();
        void setPrefitParameters(const RooArgList& prefitParams);
        void setCorrectedParameters(const RooArgList& corrfitParams);
        ROOT::Math::PxPyPzMVector createLorentzVector(double const& phi, double const& y, double const& pt, double const& m);
        RooRealVar *getYieldInRange(RooFitResult *fitResult, RooRealVar *massCand1, RooRealVar *massCand2, RooProdPdf function, RooFormulaVar nCands, TString range);
        // Fitting and plotting functions
        void do2DFit(Bool_t draw, Bool_t doReflections, Bool_t isMc, TFile *fout);
        void plotProjectionsAfterFit(RooFitResult *fitResult, RooProdPdf *model, RooDataSet *dataset, TString saveName, TFile *fout, bool doReflections, const char* suffix);
        void plot2DFit(TH2D *hMassCorrelations, TH2D* histFit, RooProdPdf *model, Bool_t draw, TFile *fout, TString const& cName);
        void plotFitResults(RooDataSet* dataset, RooRealVar* mass, RooAbsPdf* model, RooFitResult* fitResult, const char* sgnComponent, const char* bkgComponent, const char* title, const char* canvasName);
        void setHistoSignalSidebandStyle(TH1F *hSideband, TH1F *hSignal, int const& candNum, TString physVar);
        void setHistoSignalSidebandStyle(TH1F *hSideband, TH1F *hSignal, TString physVar);
        void plotKinematicDistributions(TH1F* histSidebandCand1, TH1F* histSignalCand1, TH1F* histSidebandCand2, TH1F* histSignalCand2,
                                        float const nSideband, float const nSignal, TString const varName, TString canvasName, TFile* fout);
        void plotDeltaKinematicDistributions(TH1F* histSidebandDeltaPt, TH1F* histSignalDeltaPt, TH1F* histSidebandDeltaY,
                                            TH1F* histSignalDeltaY, TH1F* histSignalDeltaPhi, TH1F* histSidebandDeltaPhi,
                                            float const nSideband, float const nSignal, TString canvasName, TFile *fout);

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
        float _ptMinPair = 0.;
        float _ptMaxPair = 50.;

        float ptCand1 = 0., ptCand2 = 0., yCand1 = 0., yCand2 = 0., phiCand1 = 0., phiCand2 = 0., mDCand1 = 0., mDCand2 = 0., mDbarCand1 = 0., mDbarCand2 = 0.;
        uint8_t typeCand1 = 0, typeCand2 = 0, typePair = 0;

        int _typeOfBkgPdf;
        int _typeOfSgnPdf;
        int _typeOfReflPdf;
        TString _sgnFuncOption;
        TString _bkgFuncOption;
        TString _reflFuncOption;
        bool _removeAmbiguous;

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
        RooRealVar rooPtPair;

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
