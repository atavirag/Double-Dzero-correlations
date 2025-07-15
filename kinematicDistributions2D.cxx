// Implementation of InvMassFitter2D class
// Author: Andrea Tavira García, IJCLab (Orsay, France)

#include "InvMassFitter2D.h"

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "RooPlot.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TBox.h"
#include "RooCurve.h"
#include "TLegend.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include <cstring> // for strcmp
#include "Math/Vector4D.h"
#include "TH2D.h"
#include "TLine.h"
#include "RooHist.h"

using namespace RooFit;

// Function to analyse and compare the pT and y distributions of my candidates
// in a sideband region and the signal region after subtracting the background
void InvMassFitter2D::analyseKinematicDistributions(TFile *fout, bool isWeighted, const char *suffix)
{

    RooDataSet *dataset = dynamic_cast<RooDataSet *>(_workspace.data("weightedData"));

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
    RooRealVar *phiCand1 = dynamic_cast<RooRealVar *>(dataset->get()->find("fPhiCand1"));
    RooRealVar *phiCand2 = dynamic_cast<RooRealVar *>(dataset->get()->find("fPhiCand2"));
    RooRealVar *ptPair = dynamic_cast<RooRealVar *>(dataset->get()->find("fPtPair"));
    RooRealVar *weightCand1 = dynamic_cast<RooRealVar *>(dataset->get()->find("weightCand1"));
    RooRealVar *weightCand2 = dynamic_cast<RooRealVar *>(dataset->get()->find("weightCand2"));

    massCand1->setRange(_massMin, _massMax);
    massCand2->setRange(_massMin, _massMax);

    // Create formula variables for deltaPt and deltaY
    RooRealVar deltaPt("deltaPt", "ptCand1 - ptCand2", 0., -24.0, 24.0);
    RooRealVar deltaY("deltaY", "yCand1 - yCand2", 0., -2.0, 2.0);
    RooRealVar deltaPhi("deltaPhi", "phiCand1 - phiCand2", 0., -10.0, 10.0);

    // Create a set of the variables for the new dataset
    RooArgSet vars(*massCand1, *massCand2, *ptCand1, *ptCand2, *yCand1, *yCand2, *phiCand1, *phiCand2, *ptPair, deltaPt, deltaY, deltaPhi);
    RooArgSet weightedVars(*massCand1, *massCand2, *ptCand1, *ptCand2, *yCand1, *yCand2, *phiCand1, *phiCand2, *ptPair, deltaPt, deltaY, deltaPhi, *weightCand1, *weightCand2);

    RooRealVar *mean = _workspace.var("mean");
    RooRealVar *sigma = _workspace.var("sigma");
    // RooRealVar *meanCand2 = _workspace.var("meanCand2");
    // RooRealVar *sigmaCand2 = _workspace.var("sigmaCand2");

    TH1F *hPtPair = new TH1F("hPtPair", Form("Pair Pt distribution_%s", suffix), 100, 0, 15);
    TH2F *hPtPairVsDeltaPt = new TH2F(Form("ptPairVsDeltaPtHist_%s", suffix), "Pair Pt vs. delta Pt distribution", 96, 0, 24, 30, -1.5, 1.5);
    TH2F *hPtPairVsDeltaY = new TH2F(Form("ptPairVsDeltaYHist_%s", suffix), "Pair Pt vs. delta Y distribution", 96, 0, 24, 30, -1., 1.);
    TH2F *hPtPairVsDeltaPhi = new TH2F(Form("ptPairVsDeltaPhiHist_%s", suffix), "Pair Pt vs. delta Phi distribution", 96, 0, 24, 28, -6.28, 6.28);

    TH2F *hPtPairVsDeltaPt_02 = new TH2F(Form("ptPairVsDeltaPtHist_0-2_%s", suffix), "Pair Pt vs. delta Pt distribution", 96, 0, 24, 75, -1.5, 1.5);
    TH2F *hPtPairVsDeltaY_02 = new TH2F(Form("ptPairVsDeltaYHist_0-2_%s", suffix), "Pair Pt vs. delta Y distribution", 96, 0, 24, 30, -1.0, 1.0);
    TH2F *hPtPairVsDeltaPhi_02 = new TH2F(Form("ptPairVsDeltaPhiHist_0-2_%s", suffix), "Pair Pt vs. delta Phi distribution", 96, 0, 24, 28, -6.28, 6.28);

    TH2F *hPtPairVsAbsDeltaPt = new TH2F(Form("ptPairVsAbsDeltaPtHist_%s", suffix), "Pair Pt vs. |delta Pt| distribution", 96, 0, 24, 45, 0, 1.5);
    TH2F *hPtPairVsAbsDeltaY = new TH2F(Form("ptPairVsAbsDeltaYHist_%s", suffix), "Pair Pt vs. |delta Y| distribution", 96, 0, 24, 60, 0., 1.5);
    TH2F *hPtPairVsAbsDeltaPhi = new TH2F(Form("ptPairVsAbsDeltaPhiHist_%s", suffix), "Pair Pt vs. |delta Phi| distribution", 72, 0, 36, 28, 0., 3.14);

    TH2F *hPtPairVsAbsDeltaPt_02 = new TH2F(Form("ptPairVsAbsDeltaPtHist_0-2_%s", suffix), "Pair Pt vs. |delta Pt| distribution", 96, 0, 24, 45, 0, 1.5);
    TH2F *hPtPairVsAbsDeltaY_02 = new TH2F(Form("ptPairVsAbsDeltaYHist_0-2_%s", suffix), "Pair Pt vs. |delta Y| distribution", 96, 0, 24, 60, 0., 1.5);
    TH2F *hPtPairVsAbsDeltaPhi_02 = new TH2F(Form("ptPairVsAbsDeltaPhiHist_0-2_%s", suffix), "Pair Pt vs. |delta Phi| distribution", 72, 0, 36, 28, 0., 3.14);

    // Create dataset in signal region
    // Manually fill the dataset with deltaPt and deltaY values
    // TODO: remove background candidates from signal region

    RooDataSet *datasetSignalRegion = new RooDataSet("datasetSignalRegion", "datasetSignalRegion", vars);

    for (int i = 0; i < dataset->numEntries(); ++i)
    {
        const RooArgSet *row = dataset->get(i);
        // Select signal region
        if (massCand1->getVal() < (mean->getVal() - 3 * sigma->getVal()) || massCand1->getVal() > (mean->getVal() + 3 * sigma->getVal()))
        {
            continue;
        }
        if (massCand2->getVal() < (mean->getVal() - 3 * sigma->getVal()) || massCand2->getVal() > (mean->getVal() + 3 * sigma->getVal()))
        {
            continue;
        }

        massCand1->setVal(row->getRealValue("fMCand1"));
        massCand2->setVal(row->getRealValue("fMCand2"));
        ptCand1->setVal(row->getRealValue("fPtCand1"));
        ptCand2->setVal(row->getRealValue("fPtCand2"));
        yCand1->setVal(row->getRealValue("fYCand1"));
        yCand2->setVal(row->getRealValue("fYCand2"));
        phiCand1->setVal(row->getRealValue("fPhiCand1"));
        phiCand2->setVal(row->getRealValue("fPhiCand2"));
        ptPair->setVal(row->getRealValue("fPtPair"));
        weightCand1->setVal(row->getRealValue("weightCand1"));
        weightCand2->setVal(row->getRealValue("weightCand2"));

        // Calculate deltaPt, deltaPhi, and deltaY
        double deltaPtValue = (ptCand1->getVal() - ptCand2->getVal()) / ptPair->getVal();
        double deltaYValue = yCand1->getVal() - yCand2->getVal();
        double deltaPhiValue = phiCand1->getVal() - phiCand2->getVal();
        double combinedWeight = weightCand1->getVal() * weightCand2->getVal();

        double wrappedDphi = std::fabs(std::atan2(std::sin(deltaPhiValue), std::cos(deltaPhiValue)));

        // Set the delta values manually in the new dataset
        deltaPt.setVal(deltaPtValue);
        deltaY.setVal(deltaYValue);
        deltaPhi.setVal(wrappedDphi);

        datasetSignalRegion->add(vars);

        if (isWeighted)
        {
            hPtPair->Fill(ptPair->getVal(), combinedWeight); // Fill histogram for each pair
            hPtPairVsDeltaPt->Fill(ptPair->getVal(), deltaPtValue, combinedWeight);
            hPtPairVsDeltaY->Fill(ptPair->getVal(), deltaYValue, combinedWeight);
            hPtPairVsDeltaPhi->Fill(ptPair->getVal(), deltaPhiValue, combinedWeight);

            hPtPairVsAbsDeltaPt->Fill(ptPair->getVal(), abs(deltaPtValue), combinedWeight);
            hPtPairVsAbsDeltaY->Fill(ptPair->getVal(), abs(deltaYValue), combinedWeight);
            hPtPairVsAbsDeltaPhi->Fill(ptPair->getVal(), wrappedDphi, combinedWeight);

            if (ptPair->getVal() >= 4)
            {
                hPtPairVsDeltaPt_02->Fill(ptPair->getVal(), deltaPtValue, combinedWeight);
                hPtPairVsDeltaY_02->Fill(ptPair->getVal(), deltaYValue, combinedWeight);
                hPtPairVsDeltaPhi_02->Fill(ptPair->getVal(), deltaPhiValue, combinedWeight);

                hPtPairVsAbsDeltaPt_02->Fill(ptPair->getVal(), abs(deltaPtValue), combinedWeight);
                hPtPairVsAbsDeltaY_02->Fill(ptPair->getVal(), abs(deltaYValue), combinedWeight);
                hPtPairVsAbsDeltaPhi_02->Fill(ptPair->getVal(), wrappedDphi, combinedWeight);
            }
        }
        else
        {
            hPtPair->Fill(ptPair->getVal()); // Fill histogram for each pair
            hPtPairVsDeltaPt->Fill(ptPair->getVal(), deltaPtValue);
            hPtPairVsDeltaY->Fill(ptPair->getVal(), deltaYValue);
            hPtPairVsDeltaPhi->Fill(ptPair->getVal(), deltaPhiValue);

            hPtPairVsAbsDeltaPt->Fill(ptPair->getVal(), abs(deltaPtValue));
            hPtPairVsAbsDeltaY->Fill(ptPair->getVal(), abs(deltaYValue));
            hPtPairVsAbsDeltaPhi->Fill(ptPair->getVal(), wrappedDphi);

            if (ptPair->getVal() >= 4)
            {
                hPtPairVsDeltaPt_02->Fill(ptPair->getVal(), deltaPtValue);
                hPtPairVsDeltaY_02->Fill(ptPair->getVal(), deltaYValue);
                hPtPairVsDeltaPhi_02->Fill(ptPair->getVal(), deltaPhiValue);

                hPtPairVsAbsDeltaPt_02->Fill(ptPair->getVal(), abs(deltaPtValue));
                hPtPairVsAbsDeltaY_02->Fill(ptPair->getVal(), abs(deltaYValue));
                hPtPairVsAbsDeltaPhi_02->Fill(ptPair->getVal(), wrappedDphi);
            }
        }
    }

    double binWidth = 1.0;
    int nBins = (int)((_ptMax - _ptMin) / binWidth);

    // Create candidate histograms
    TH1D *histPtCand1 = (TH1D *)dataset->createHistogram(Form("histPtCand1_%s", suffix),
                                                         *ptCand1, Binning(48, _ptMin, _ptMax));
    TH1D *histPtCand2 = (TH1D *)dataset->createHistogram(Form("histPtCand2_%s", suffix),
                                                         *ptCand2, Binning(48, _ptMin, _ptMax));

    // Configure and plot candidate histograms
    TCanvas *cPtCands = new TCanvas(Form("cPtCands_%s", suffix), Form("cPtCands_%s", suffix), 800, 600);
    cPtCands->SetLogy();

    configureHistogram(histPtCand1, 20, kBlue);
    configureHistogram(histPtCand2, 21, kRed);

    histPtCand1->SetTitle("#it{p}_{T} of the D-meson candidate");
    histPtCand1->GetXaxis()->SetTitle("#it{p}_{T}^{D} (GeV/#it{c})");
    histPtCand1->GetYaxis()->SetTitle("Entries");
    histPtCand1->Draw("PE");
    histPtCand2->Draw("samePE");

    TLegend *legend_ptCand = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_ptCand->AddEntry(histPtCand1, "Candidate 1", "ep");
    legend_ptCand->AddEntry(histPtCand2, "Candidate 2", "ep");
    legend_ptCand->Draw();

    cPtCands->SaveAs(Form("ptCands_%s.png", suffix));

    // Create pair histogram
    TCanvas *cPtPair = new TCanvas(Form("cPtPair_%s", suffix), Form("cPtPair_%s", suffix), 800, 600);
    hPtPair->SetTitle("#it{p}_{T} of the D-meson pair");
    hPtPair->GetXaxis()->SetTitle("#it{p}_{T}^{DD} (GeV/#it{c})");
    hPtPair->GetYaxis()->SetTitle("Entries");
    hPtPair->Draw();
    cPtPair->SaveAs(Form("ptPairHist_%s.png", suffix));

    // Create 2D histogram canvases using helper function
    TCanvas *cPtPairVsDeltaPt = create2DCanvas(hPtPairVsDeltaPt, "cPtPairVsDeltaPt", suffix,
                                               "#it{p}_{T}^{DD}(GeV/#it{c})", "#Delta#it{p}_{T} (GeV/#it{c})", false, false);

    TCanvas *cPtPairVsDeltaPhi = create2DCanvas(hPtPairVsDeltaPhi, "cPtPairVsDeltaPhi", suffix,
                                                "#it{p}_{T}^{DD}(GeV/#it{c})", "#Delta#varphi", false, false);

    TCanvas *cPtPairVsDeltaY = create2DCanvas(hPtPairVsDeltaY, "cPtPairVsDeltaY", suffix,
                                              "#it{p}_{T}^{DD}(GeV/#it{c})", "#Delta#it{y}", false, true);

    TCanvas *cPtPairVsAbsDeltaPt = create2DCanvas(hPtPairVsAbsDeltaPt, "cPtPairVsAbsDeltaPt", suffix,
                                                  "#it{p}_{T}^{DD}(GeV/#it{c})", "|#Delta#it{p}_{T}|/#it{p}_{T}^{DD}", true, false);
    cPtPairVsAbsDeltaPt->SaveAs(Form("ptPairVsAbsDeltaPtHist_%s.png", suffix));

    TCanvas *cPtPairVsAbsDeltaPhi = create2DCanvas(hPtPairVsAbsDeltaPhi, "cPtPairVsAbsDeltaPhi", suffix,
                                                   "#it{p}_{T}^{DD}(GeV/#it{c})", "|#Delta#it{#varphi}| (rad)", false, false);

    // Add white box to delta phi plot
    TBox *box = new TBox(4, 4, 6, 6);
    box->SetFillColor(kWhite);
    box->SetLineColor(kWhite);
    box->Draw("same");
    cPtPairVsAbsDeltaPhi->SaveAs(Form("ptPairVsAbsDeltaPhiHist_%s.png", suffix));
    cPtPairVsAbsDeltaPhi->SetLogz();

    TCanvas *cPtPairVsAbsDeltaY = create2DCanvas(hPtPairVsAbsDeltaY, "cPtPairVsAbsDeltaY", suffix,
                                                 "#it{p}_{T}^{DD}(GeV/#it{c})", "|#Delta#it{y}|", true, false);
    cPtPairVsAbsDeltaY->SaveAs(Form("ptPairVsAbsDeltaYHist_%s.png", suffix));

    // Create similar canvases for _02 histograms
    TCanvas *cPtPairVsDeltaPt_02 = create2DCanvas(hPtPairVsDeltaPt_02, "cPtPairVsDeltaPt_02", suffix,
                                                  "#it{p}_{T}^{Pair}(GeV/#it{c})", "#Delta#it{p}_{T} (GeV/#it{c})", false, false);

    TCanvas *cPtPairVsDeltaPhi_02 = create2DCanvas(hPtPairVsDeltaPhi_02, "cPtPairVsDeltaPhi_02", suffix,
                                                   "#it{p}_{T}^{Pair}(GeV/#it{c})", "#Delta#varphi", false, false);

    TCanvas *cPtPairVsDeltaY_02 = create2DCanvas(hPtPairVsDeltaY_02, "cPtPairVsDeltaY_02", suffix,
                                                 "#it{p}_{T}^{D}(GeV/#it{c})", "#Delta#it{y}", false, false);

    TCanvas *cPtPairVsAbsDeltaPt_02 = create2DCanvas(hPtPairVsAbsDeltaPt_02, "cPtPairVsAbsDeltaPt_02", suffix,
                                                     "#it{p}_{T}^{D}(GeV/#it{c})", "#Delta#it{p}_{T} (GeV/#it{c})", false, false);

    TCanvas *cPtPairVsAbsDeltaPhi_02 = create2DCanvas(hPtPairVsAbsDeltaPhi_02, "cPtPairVsAbsDeltaPhi_02", suffix,
                                                      "#it{p}_{T}^{D}(GeV/#it{c})", "|#Delta#it{#varphi}| (rad)", false, false);

    TCanvas *cPtPairVsAbsDeltaY_02 = create2DCanvas(hPtPairVsAbsDeltaY_02, "cPtPairVsAbsDeltaY_02", suffix,
                                                    "#it{p}_{T}^{D}(GeV/#it{c})", "|#Delta#it{y}|", false, false);

    // Create projections using helper function
    TH1D *hDeltaPt = createProjectionWithLabels(hPtPairVsAbsDeltaPt, Form("hDeltaAbsPt_%s", suffix),
                                                "#Delta p_{T} distribution for %.1f < p_{T}^{DD} < %.1f GeV/#it{c}",
                                                _ptMin, _ptMax, "|#Delta#it{p}_{T}|/#it{p}_{T}^{DD}");

    TH1D *hDeltaPhi = createProjectionWithLabels(hPtPairVsAbsDeltaPhi, Form("hDeltaAbsPhi_%s", suffix),
                                                 "#Delta #varphi distribution for %.1f < p_{T}^{Pair} < %.1f GeV/#it{c}",
                                                 _ptMin, _ptMax, "|#Delta#varphi| (rad)");

    TH1D *hDeltaY = createProjectionWithLabels(hPtPairVsAbsDeltaY, Form("hDeltaAbsY_%s", suffix),
                                               "#Delta #it{y} distribution for %.1f < p_{T}^{DD} < %.1f GeV/#it{c}",
                                               _ptMin, _ptMax, "|#Delta#it{y}|");

    // Similar projections for _02 histograms
    TH1D *hDeltaPt_02 = createProjectionWithLabels(hPtPairVsAbsDeltaPt_02, Form("hDeltaAbsPt_02_%s", suffix),
                                                   "#Delta p_{T} distribution for %.1f < p_{T}^{D} < %.1f GeV/#it{c}",
                                                   4.0, _ptMax, "|#Delta#it{p}_{T}| (GeV/#it{c})");

    TH1D *hDeltaPhi_02 = createProjectionWithLabels(hPtPairVsAbsDeltaPhi_02, Form("hDeltaAbsPhi_02_%s", suffix),
                                                    "#Delta #varphi distribution for %.1f < p_{T}^{D} < %.1f GeV/#it{c}",
                                                    4.0, _ptMax, "#|Delta#varphi|");

    TH1D *hDeltaY_02 = createProjectionWithLabels(hPtPairVsAbsDeltaY_02, Form("hDeltaAbsY_02_%s", suffix),
                                                  "#Delta #it{y} distribution for %.1f < p_{T}^{D} < %.1f GeV/#it{c}",
                                                  4.0, _ptMax, "|#Delta#it{y}|");

    // Plot Delta pT comparison
    TCanvas *cDeltaPt = new TCanvas(Form("cDeltaPt_%s", suffix), "Delta pT", 800, 600);
    cDeltaPt->SetLogy();
    cDeltaPt->SetTicks(1, 1);

    configureHistogram(hDeltaPt, 21, kBlue);
    configureHistogram(hDeltaPt_02, 22, kRed);
    hDeltaPt->SetMinimum(0.5);

    hDeltaPt->Draw("EP");
    hDeltaPt_02->Draw("EPsame");
    addText(false);

    TLegend *legend_deltaPt = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend_deltaPt->SetBorderSize(0);
    legend_deltaPt->SetFillStyle(0);
    legend_deltaPt->SetTextSize(0.04);
    legend_deltaPt->AddEntry(hDeltaPt, "p_{T}^{DD} min. = 0", "lp");
    legend_deltaPt->AddEntry(hDeltaPt_02, "p_{T}^{DD} min. = 4 GeV/#it{c}", "lp");
    legend_deltaPt->Draw();

    cDeltaPt->SaveAs(Form("deltaPtDistribution_%s.png", suffix));

    // Plot Delta phi comparison
    TCanvas *cDeltaPhi = new TCanvas(Form("cDeltaPhi_%s", suffix), "Delta phi", 800, 600);
    cDeltaPhi->SetLeftMargin(0.12);
    cDeltaPhi->SetRightMargin(0.02);
    cDeltaPhi->SetBottomMargin(0.12);
    cDeltaPhi->SetTicks(1, 1);

    configureHistogram(hDeltaPhi, 21, kBlue);
    configureHistogram(hDeltaPhi_02, 22, kRed);

    hDeltaPhi->Draw("EP");
    hDeltaPhi_02->Draw("EPsame");
    addText(true);

    TLegend *legend_deltaPhi = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend_deltaPhi->SetBorderSize(0);
    legend_deltaPhi->SetFillStyle(0);
    legend_deltaPhi->SetTextSize(0.04);
    legend_deltaPhi->AddEntry(hDeltaPhi, "#it{p}_{T}^{Pair} min. = 0", "lp");
    legend_deltaPhi->AddEntry(hDeltaPhi_02, "#it{p}_{T}^{Pair} min. = 4 GeV/#it{c}", "lp");
    legend_deltaPhi->Draw();

    cDeltaPhi->SaveAs(Form("deltaPhiDistribution_%s.png", suffix));

    // Plot Delta Y comparison
    TCanvas *cDeltaY = new TCanvas(Form("cDeltaY_%s", suffix), "Delta Y", 800, 600);
    cDeltaY->SetTicks(1, 1);

    configureHistogram(hDeltaY, 21, kBlue);
    configureHistogram(hDeltaY_02, 22, kRed);

    hDeltaY->Draw("EP");
    hDeltaY_02->Draw("EPsame");
    addText(false);

    TLegend *legend_deltaY = new TLegend(0.65, 0.75, 0.85, 0.85);
    legend_deltaY->SetBorderSize(0);
    legend_deltaY->SetFillStyle(0);
    legend_deltaY->SetTextSize(0.04);
    legend_deltaY->AddEntry(hDeltaY, "p_{T}^{DD} min. = 0", "lp");
    legend_deltaY->AddEntry(hDeltaY_02, "p_{T}^{DD} min. = 4 GeV/#it{c}", "lp");
    legend_deltaY->Draw();

    cDeltaY->SaveAs(Form("deltaYDistribution_%s.png", suffix));

    // Write to file
    fout->cd();

    // Write canvases
    std::vector<TCanvas *> canvasesToWrite = {
        cDeltaPt, cDeltaPhi, cDeltaY, cPtPair, cPtPairVsDeltaPt, cPtPairVsDeltaY,
        cPtPairVsDeltaPhi, cPtPairVsAbsDeltaPt, cPtPairVsAbsDeltaY, cPtPairVsAbsDeltaPhi, cPtCands};

    for (auto *canvas : canvasesToWrite)
    {
        canvas->Write();
    }

    // Write histograms
    std::vector<TH1 *> histsToWrite = {
        hPtPairVsDeltaPt, hPtPairVsDeltaY, hPtPairVsDeltaPhi,
        hPtPairVsAbsDeltaPt, hPtPairVsAbsDeltaY, hPtPairVsAbsDeltaPhi,
        hPtPairVsDeltaPt_02, hPtPairVsDeltaY_02, hPtPairVsDeltaPhi_02,
        hPtPairVsAbsDeltaPt_02, hPtPairVsAbsDeltaY_02, hPtPairVsAbsDeltaPhi_02,
        histPtCand1, histPtCand2};

    for (auto *hist : histsToWrite)
    {
        hist->Write();
    }

    // Clean up memory
    for (auto *canvas : canvasesToWrite)
    {
        delete canvas;
    }

    // Delete additional canvases
    std::vector<TCanvas *> additionalCanvases = {
        cPtPairVsDeltaPt_02, cPtPairVsDeltaPhi_02, cPtPairVsDeltaY_02,
        cPtPairVsAbsDeltaPt_02, cPtPairVsAbsDeltaPhi_02, cPtPairVsAbsDeltaY_02};

    for (auto *canvas : additionalCanvases)
    {
        delete canvas;
    }

    for (auto *hist : histsToWrite)
    {
        delete hist;
    }
}

void InvMassFitter2D::addText(bool drawStats)
{
    // Helper lambda to create a styled TPaveText
    auto createTPave = [](double x1, double y1, double x2, double y2, int align = 11)
    {
        TPaveText *pave = new TPaveText(x1, y1, x2, y2, "brNDC");
        pave->SetTextFont(42);
        pave->SetTextSize(0.04);
        pave->SetBorderSize(0);
        pave->SetFillStyle(0);
        pave->SetTextAlign(align);
        return pave;
    };

    // ALICE Preliminary (Title)
    TPaveText *ptTitle = createTPave(0.2, 0.89, 0.85, 0.94);
    ptTitle->AddText("#scale[1.35]{ALICE Preliminary}");

    // Collision System
    TPaveText *ptCol = createTPave(0.55, 0.89, 0.85, 0.94);
    ptCol->AddText("#scale[1.2]{pp, #sqrt{#it{s}} = 13.6 TeV}");

    // Decay channel
    TPaveText *ptDecay = createTPave(0.2, 0.75, 0.87, 0.82);
    ptDecay->AddText("D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.");

    // Pair type
    TPaveText *ptType = createTPave(0.2, 0.69, 0.87, 0.74);
    if (strcmp(_pairType, "OS") == 0)
    {
        ptType->AddText("D^{0}#bar{D}^{0} + #bar{D}^{0}D^{0} pairs");
    }
    else if (strcmp(_pairType, "LS") == 0)
    {
        ptType->AddText("D^{0}D^{0} + #bar{D}^{0}#bar{D}^{0} pairs");
    }

    // pT range
    TPaveText *ptPt = createTPave(0.2, 0.63, 0.87, 0.68);
    ptPt->AddText(Form("%.1f < #it{p}^{D}_{T} < %.1f GeV/#it{c}", _ptMin, _ptMax));

    // Statistical uncertainty note
    TPaveText *ptStats = createTPave(0.2, 0.72, 0.4, 0.95);
    ptStats->AddText("Statistical uncertainties only");

    // Mass window cut
    TPaveText *ptMass = createTPave(0.2, 0.56, 0.87, 0.61);
    ptMass->AddText("|#it{m}_{K#pi} #minus #it{m}_{D^{0}}| < 3#it{#sigma}");

    ptTitle->Draw("SAME");
    ptType->Draw("same");
    ptPt->Draw("same");
    ptCol->Draw("same");
    ptMass->Draw("same");

    if (drawStats)
    {
        ptStats->Draw("same");
    }
}

TCanvas *InvMassFitter2D::create2DCanvas(TH2F *hist, const std::string &canvasName, const char *suffix,
                                         const std::string &xTitle, const std::string &yTitle,
                                         bool setLogz = false, bool setLogy = false)
{
    TCanvas *canvas = new TCanvas(Form("%s_%s", canvasName.c_str(), suffix),
                                  Form("%s_%s", canvasName.c_str(), suffix), 800, 600);

    if (setLogz)
        canvas->SetLogz();
    if (setLogy)
        canvas->SetLogy();

    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetXaxis()->SetTitleSize(0.05);
    hist->GetYaxis()->SetTitleSize(0.05);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->SetTitle("");
    hist->SetStats(0);
    hist->Draw("colz");

    canvas->SetTicks(1, 1);
    addText(false);

    return canvas;
}

void InvMassFitter2D::configureHistogram(TH1D *hist, int markerStyle, int color)
{
    hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetTitle("");
    hist->SetStats(0);
}

TH1D *InvMassFitter2D::createProjectionWithLabels(TH2F *hist2D, const char *name,
                                                  const char *titleFormat, double ptMin, double ptMax,
                                                  const char *axisLabel)
{
    TH1D *proj = hist2D->ProjectionY(name);
    proj->SetTitle(Form(titleFormat, ptMin, ptMax));
    proj->GetXaxis()->SetTitle(axisLabel);
    proj->GetYaxis()->SetTitle("Counts");
    return proj;
}

void InvMassFitter2D::saveCanvasIfNeeded(TCanvas *canvas, const std::string &filename, const char *suffix)
{
    if (!filename.empty())
    {
        canvas->SaveAs(Form("%s_%s.png", filename.c_str(), suffix));
    }
}