#include "TStyle.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "Riostream.h"
#include "TMath.h"
#include "TSystem.h"

TGraphErrors *getGraphFromHisto(TH1F *h, bool is2024) {
    auto g = new TGraphErrors();
    for (int i = 1; i <= h->GetNbinsX(); i++ ) {
        // Get bin content and error
        double y = h->GetBinContent(i);
        double y_err = h->GetBinError(i);

        // Calculate x-value (bin center)
        double x = h->GetBinCenter(i);
        double x_err = h->GetBinWidth(i)/2;
        //g->SetPoint(i-1, hRun3->GetBinCenter(i), hRun3->GetBinContent(i));
        // Set point in the TGraphErrors object
        g->SetPoint(i-1, x, y);
        g->SetPointError(i-1, x_err, y_err);
    }
    return g;
}
bool verbose = true;
void compareRun2() {
    TString const fileNameRun2 = "/home/andrea/MyMacros/Cut_variation/Run2Results.root";
    TString const fileNamePrelim = "/home/andrea/MyMacros/Cut_variation/CutVarD0_pp136TeV_final.root";
    TString const fileName2024 = "/home/andrea/MyMacros/Double-Dzero-correlations/cut_variation/CutVarDplus_2024_triggered.root";
    TString const outName = "comparisonWPreliminary.root";

    // Cargamos los datos
    TFile *fRun2 = TFile::Open(fileNameRun2);
    TFile *fPrelim = TFile::Open(fileNamePrelim);
    TFile *f2024 = TFile::Open(fileName2024);

    if (!fRun2 || !fPrelim || !f2024) {
        cerr << "ERROR: files not found\n";
        return;
    }

    // Archivo que me junta todos los plots que necesito
    TFile* outputFile = new TFile(outName, "RECREATE");
    outputFile->cd();
    TString dirnameRun2 = "Table 1";
    TDirectory *dirRun2 = (TDirectory *)fRun2->Get(dirnameRun2);

    // Cogemos los histogramas del MC
    TGraph *gRun2 = (TGraph *)dirRun2->Get("Graph1D_y5");
    TH1F *hPrelim = (TH1F *)fPrelim->Get("hCorrFracNonPrompt");
    TH1F *h2024 = (TH1F *)f2024->Get("hCorrFracNonPrompt");

    // Transformamos el histograma en un TGraph
    auto gPrelim = getGraphFromHisto(hPrelim, false);
    auto g2024 = getGraphFromHisto(h2024, true);

    if (verbose) {
        cout << "            Data loaded \n";
    }
    TCanvas *canvas = new TCanvas("canvas", "TGraph Example", 800, 600);
    canvas->SetGrid();

    gPrelim->SetMarkerColor(kBlue-3);
    gPrelim->SetMarkerStyle(23);
    gPrelim->SetLineColor(kBlue-3);
    gPrelim->SetLineWidth(2);

    g2024->SetMarkerColor(kGreen+2);
    g2024->SetMarkerStyle(21);
    g2024->SetLineColor(kGreen+2);
    g2024->SetLineWidth(2);

    gRun2->SetMarkerStyle(22);
    gRun2->SetMarkerColor(kRed);
    gRun2->SetLineColor(kRed);
    gRun2->SetLineWidth(2);

    gRun2->GetYaxis()->SetRangeUser(0.0, 0.15);
    gRun2->Draw("PA");
    gPrelim->Draw("Psame");
    g2024->Draw("Psame");

    gRun2->SetTitle("Non-prompt fraction");
    gRun2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    gRun2->GetYaxis()->SetTitle("Non-prompt fraction");

    TLegend* legend = new TLegend(0.5, 0.2, 0.85, 0.5); // Adjust coordinates as needed
    legend->AddEntry(gPrelim, "Run 3 Preliminary", "lep");
    legend->AddEntry(g2024, "Run 3 result 2024 sample", "lep");
    legend->AddEntry(gRun2, "Run 2", "lep");
    legend->Draw();

    canvas->Draw();
    canvas->SaveAs("Run2Comparison.png");

    /* // calculate the ratio wrt the updated results
    TH1F *hRatioPtSmearing = (TH1F *)hRun3PtSmearing->Clone();
    hRatioPtSmearing->Divide(hRun3New);
    TH1F *hRatioWorseReso = (TH1F *)hRun3WorseReso->Clone();
    hRatioWorseReso->Divide(hRun3New);

    TLine *line = new TLine(0, 1, 24, 1);

    TCanvas *cRatio = new TCanvas("canvas", "TGraph Example", 800, 600);
    cRatio->SetGrid();

    line->SetLineColor(kBlack);
    line->SetLineWidth(2);

    hRatioPtSmearing->SetMarkerColor(kMagenta-4);
    hRatioPtSmearing->SetMarkerStyle(20);
    hRatioPtSmearing->SetLineColor(kMagenta-4);
    hRatioPtSmearing->SetLineWidth(2);

    hRatioWorseReso->SetMarkerColor(kBlue-3);
    hRatioWorseReso->SetMarkerStyle(23);
    hRatioWorseReso->SetLineColor(kBlue-3);
    hRatioWorseReso->SetLineWidth(2);


    hRatioPtSmearing->GetYaxis()->SetRangeUser(0.6, 1.4);
    hRatioPtSmearing->Draw("PE");
    hRatioWorseReso->Draw("Psame");
    line->Draw("Lsame");
    //g2024->Draw("Psame");
    //gMingyu->Draw("Psame");
    hRatioPtSmearing->SetTitle("Ratio of non-prompt fractions wrt Run 3 result");
    hRatioPtSmearing->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hRatioPtSmearing->GetYaxis()->SetTitle("Non-prompt fraction ratio");

    TLegend* lRatio = new TLegend(0.35, 0.1, 0.65, 0.3); // Adjust coordinates as needed
    lRatio->AddEntry(hRatioPtSmearing, "Run 3 w/pt smearing", "lep");
    //legend->AddEntry(gMingyu, "Updated Run 3 result Mingyu", "lep");
    lRatio->AddEntry(hRatioWorseReso, "Run 3 w/worse reso", "lep");
    //legend->AddEntry(gRun2, "Run 2", "lep");
    lRatio->Draw();

    cRatio->Draw();
    cRatio->SaveAs("comparisonRatio.png"); */

}
