#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <string>

void computeTrackTunerSys() {
    TFile* fmain = TFile::Open("~/MyMacros/Common/AccEffPreselD0ToKPi_bothBDTs_no_ambiguous.root", "READ");
    if (!fmain || fmain->IsZombie()) {
        std::cerr << "Error: Could not open fmain" << std::endl;
        return;
    }

    TFile* fmain_full = TFile::Open("~/MyMacros/Common/AccEffPreselD0ToKPi_bothBDTs.root", "READ");
    if (!fmain_full || fmain_full->IsZombie()) {
        std::cerr << "Error: Could not open fmain" << std::endl;
        return;
    }
    TFile* fworse1p = TFile::Open("~/MyMacros/Common/AccEffPreselD0ToKPi_bothBDTs_no_ambiguous_reso_1p.root", "READ");
    if (!fworse1p || fworse1p->IsZombie()) {
        std::cerr << "Error: Could not open fworse1p" << std::endl;
        return;
    }
    TFile* fworse1p5 = TFile::Open("~/MyMacros/Common/AccEffPreselD0ToKPi_bothBDTs_no_ambiguous_worse_1p5.root", "READ");
    if (!fworse1p5 || fworse1p5->IsZombie()) {
        std::cerr << "Error: Could not open fworse1p5" << std::endl;
        return;
    }

    TFile* fworse1p_full = TFile::Open("~/MyMacros/Common/AccEffPreselD0ToKPi_bothBDTs_worse_reso.root", "READ");
    if (!fworse1p || fworse1p->IsZombie()) {
        std::cerr << "Error: Could not open fworse1p_full" << std::endl;
        return;
    }
    TFile* fworse1p5_full = TFile::Open("~/MyMacros/Common/AccEffPreselD0ToKPi_bothBDTs_worse_reso_1p5.root", "READ");
    if (!fworse1p5 || fworse1p5->IsZombie()) {
        std::cerr << "Error: Could not open fworse1p5_full" << std::endl;
        return;
    }

    TH1F* hEff_main = dynamic_cast<TH1F*>(fmain->Get("hAccEffPreselD0ToKPi_bothBDTs_no_ambiguousAll"));
    if (!hEff_main) {
        std::cerr << "Error: Histogram 'hAccEffPreselD0ToKPi_bothBDTs_no_ambiguousAll' not found in fmain" << std::endl;
        fmain->Close();
        return;
    }
    TH1F* hEff_main_full = dynamic_cast<TH1F*>(fmain_full->Get("hAccEffPreselD0ToKPi_bothBDTsAll"));
    if (!hEff_main_full) {
        std::cerr << "Error: Histogram 'hAccEffPreselD0ToKPi_bothBDTsAll' not found in fmain" << std::endl;
        fmain_full->Close();
        return;
    }
    TH1F* hEff_1p = dynamic_cast<TH1F*>(fworse1p->Get("hAccEffPreselD0ToKPi_bothBDTs_no_ambiguous_reso_1pAll"));
    if (!hEff_1p) {
        std::cerr << "Error: Histogram 'hAccEffPreselD0ToKPi_bothBDTs_no_ambiguous_worse_1pAll' not found in fmain" << std::endl;
        fmain->Close();
        return;
    }
    TH1F* hEff_1p5 = dynamic_cast<TH1F*>(fworse1p5->Get("hAccEffPreselD0ToKPi_bothBDTs_no_ambiguous_worse_1p5All"));
    if (!hEff_1p5) {
        std::cerr << "Error: Histogram 'hAccEffPreselD0ToKPi_bothBDTs_no_ambiguous_worse1p5All' not found in fmain" << std::endl;
        fmain->Close();
        return;
    }
    TH1F* hEff_1p_full = dynamic_cast<TH1F*>(fworse1p_full->Get("hAccEffPreselD0ToKPi_bothBDTs_worse_resoAll"));
    if (!hEff_1p_full) {
        std::cerr << "Error: Histogram 'hAccEffPreselD0ToKPi_bothBDTs_worse_resoAll' not found in fmain" << std::endl;
        fmain->Close();
        return;
    }
    TH1F* hEff_1p5_full = dynamic_cast<TH1F*>(fworse1p5_full->Get("hAccEffPreselD0ToKPi_bothBDTs_worse_reso_1p5All"));
    if (!hEff_1p5_full) {
        std::cerr << "Error: Histogram 'hAccEffPreselD0ToKPi_bothBDTs_worse_reso_1p5All' not found in fmain" << std::endl;
        fmain->Close();
        return;
    }

    // Create ratio histograms
    TH1F *hRatio_1p = (TH1F*)hEff_1p->Clone("hRatio_1p");
    hRatio_1p->Divide(hEff_main);

    TH1F *hRatio_1p5 = (TH1F*)hEff_1p5->Clone("hRatio_1p5");
    hRatio_1p5->Divide(hEff_main);

    TH1F *hRatio_1p_full = (TH1F*)hEff_1p_full->Clone("hRatio_1p_full");
    hRatio_1p_full->Divide(hEff_main_full);

    TH1F *hRatio_1p5_full = (TH1F*)hEff_1p5_full->Clone("hRatio_1p5_full");
    hRatio_1p5_full->Divide(hEff_main_full);

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Efficiency Ratios", 800, 600);
    c1->SetGrid();

    // Customize the ratio histograms
    hRatio_1p->SetLineColor(kRed);
    hRatio_1p->SetMarkerColor(kRed);
    hRatio_1p->SetMarkerStyle(20);
    hRatio_1p->SetTitle("Efficiency Ratios; p_{T} (GeV/c); Ratio");
    hRatio_1p->GetYaxis()->SetRangeUser(0.5, 1.5);  // Adjust range if needed

    hRatio_1p5->SetLineColor(kBlue);
    hRatio_1p5->SetMarkerColor(kBlue);
    hRatio_1p5->SetMarkerStyle(21);

    hRatio_1p_full->SetLineColor(kViolet);
    hRatio_1p_full->SetMarkerColor(kViolet);
    hRatio_1p_full->SetMarkerStyle(22);

    hRatio_1p5_full->SetLineColor(kCyan+2);
    hRatio_1p5_full->SetMarkerColor(kCyan+2);
    hRatio_1p5_full->SetMarkerStyle(23);

    // Draw histograms
    hRatio_1p->Draw("E P");
    hRatio_1p5->Draw("E P SAME");
    hRatio_1p_full->Draw("E P SAME");
    hRatio_1p5_full->Draw("E P SAME");

    // Add legend
    TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg->AddEntry(hRatio_1p, "10%% smearing, no ambiguous", "lp");
    leg->AddEntry(hRatio_1p5, "15%% smearing, no ambiguous", "lp");
    leg->AddEntry(hRatio_1p_full, "10%% smearing, all candidates", "lp");
    leg->AddEntry(hRatio_1p5_full, "15%% smearing, all candidates", "lp");
    leg->Draw();

    // Save the plot
    c1->SaveAs("efficiency_ratios.png");

    // Create a canvas for efficiency/pT
    TCanvas *c2 = new TCanvas("c2", "Efficiency / pT", 800, 600);
    c2->SetGrid();

    // Create efficiency/pT histograms
    TH1F *hEff_main_pT = (TH1F*)hEff_main->Clone("hEff_main_pT");
    TH1F *hEff_main_pT_full = (TH1F*)hEff_main_full->Clone("hEff_main_pT_full");
    TH1F *hEff_1p_pT = (TH1F*)hEff_1p->Clone("hEff_1p_pT");
    TH1F *hEff_1p5_pT = (TH1F*)hEff_1p5->Clone("hEff_1p5_pT");

    hEff_main_pT->Divide(hEff_main_pT_full);

    // Customize efficiency/pT histograms
    hEff_main_pT->SetLineColor(kBlack);
    hEff_main_pT->SetMarkerColor(kBlack);
    hEff_main_pT->SetMarkerStyle(20);
    hEff_main_pT->SetTitle("Efficiency ratio; p_{T} (GeV/c); Efficiency ratio");

    hEff_1p_pT->SetLineColor(kRed);
    hEff_1p_pT->SetMarkerColor(kRed);
    hEff_1p_pT->SetMarkerStyle(21);

    hEff_1p5_pT->SetLineColor(kBlue);
    hEff_1p5_pT->SetMarkerColor(kBlue);
    hEff_1p5_pT->SetMarkerStyle(22);

    // Draw histograms
    hEff_main_pT->Draw("E P");
    //hEff_1p_pT->Draw("E P SAME");
    //hEff_1p5_pT->Draw("E P SAME");

    // Add legend
    TLegend *leg2 = new TLegend(0.4, 0.8, 0.9, 0.9);
    leg2->AddEntry(hEff_main_pT, "Efficiency without ambiguous / efficiency with ambiguous", "lp");
    //leg2->AddEntry(hEff_1p_pT, "hEff_1p / pT", "lp");
    //leg2->AddEntry(hEff_1p5_pT, "hEff_1p5 / pT", "lp");
    leg2->Draw();

    // Save the plot
    c2->SaveAs("efficiency_over_pT.png");

    c1->Close();
    c2->Close();

    // Close files
    fmain->Close();
    fworse1p->Close();
    fworse1p5->Close();
}

