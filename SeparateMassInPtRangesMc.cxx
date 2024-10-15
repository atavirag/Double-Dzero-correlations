// Task para ajustar el pico de masa invariante del archivo de datos del Run3

#include "TStyle.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TMath.h"
#include "THnSparse.h"
#include "THStack.h"
#include "TF1.h"
#include "TSystem.h"

#include "nlohmann/json.hpp"

using namespace std;
using json = nlohmann::json;

void DrawProjection(TH1F* hD, int const& ptMin, int const& ptMax, TString filename);
void DrawProjectionSum(TH1F* hD, TH1F* hDbar, int const& ptMin, int const& ptMax, TString filename);

void SeparateMassInPtRangesMc()
{
    // Cargamos la configuración escrita en el JSON
    ifstream jsonFile("config-SeparateMassInPtRanges.json");
    if (!jsonFile.is_open()) {
        cerr << "Error: Could not open JSON file." << endl;
        return;
    }

    // Read the JSON file into a string.
    std::string jsonString;
    jsonFile.seekg(0, std::ios::end);
    jsonString.reserve(jsonFile.tellg());
    jsonFile.seekg(0, std::ios::beg);
    jsonString.assign((std::istreambuf_iterator<char>(jsonFile)), std::istreambuf_iterator<char>());

    // Parse the JSON data.
    json jsonData = json::parse(jsonString);
    bool const verbose = jsonData["verbose"];

    // Cargamos los datos
    TString const fileName = jsonData["filename"];
    TString const dirName = jsonData["dirname"];
    TString const histName = jsonData["histname"];
    TString const outputName = jsonData["outputfile"];
    TString const outputDir = jsonData["outputdir"];
    TFile *f = TFile::Open(fileName);
    TDirectory *dir = (TDirectory *)f->Get(dirName);

    if (!jsonData["configMC"]["isMC"]) {
        cerr << "Error: Enable MC." << endl;
        return;
    }
    TString const histNamePrompt = jsonData["configMC"]["histnameMcPrompt"];
    TString const histNameNonPrompt = jsonData["configMC"]["histnameMcNonPrompt"];
    TString const histNameReflections = jsonData["configMC"]["histnameReflections"];

    // Cargamos los datos

    TH2F *hMass = (TH2F *)dir->Get(histName);
    TH2F *hMassPrompt = (TH2F *)dir->Get(histNamePrompt);
    TH2F *hMassNonPrompt = (TH2F *)dir->Get(histNameNonPrompt);
    TH2F *hMassRefl = (TH2F *)dir->Get(histNameReflections);

    if (verbose) {std::cout << "Datos cargados \n";}

    // Proyectamos hMass en función del pT
    int const nbinsPt = hMassNonPrompt->GetYaxis()->GetNbins();
    if (verbose) {
        std::cout << " Bines pT: "<< nbinsPt << std::endl;
        std::cout << "   TAMAÑO BINES   \n";
        for (int i=0; i<nbinsPt; i++) {
            std::cout << i << ".-   Low: " << hMassPrompt->GetYaxis()->GetBinLowEdge(i) << "     Up: " << hMassPrompt->GetYaxis()->GetBinUpEdge(i) << std::endl;
        }
    }

    TFile *outputFile = TFile::Open(outputDir+outputName, "RECREATE");

    // Rangos de pT
    int const nPtBins = jsonData["nPtBins"]; 
    vector<double> ptBinsEdges = jsonData["ptBinsEdges"];
    TH1F *hMassPtPrompt = (TH1F *)hMassPrompt->ProjectionX("hMassPtPrompt", 1, nbinsPt);
    TH1F *hMassPtNonPrompt = (TH1F *)hMassNonPrompt->ProjectionX("hMassPtNonPrompt", 1, nbinsPt);
    std::vector<TH1F *> vMassPtPrompt = {};
    std::vector<TH1F *> vMassPtNonPrompt = {};
    std::vector<TH1F *> vMassPtRefl = {};

    for (int i=0; i < nPtBins; i++) {
        std::cout << "Patata /n" << std::endl;
        auto nameMassPrompt = "pT"+std::to_string(i)+"MassPrompt";
        auto nameMassNonPrompt = "pT"+std::to_string(i)+"MassNonPrompt";
        auto nameReflMass = "pT"+std::to_string(i)+"ReflMass";

        vMassPtPrompt.push_back((TH1F *)hMassPrompt->ProjectionX(nameMassPrompt.c_str(), hMassPrompt->GetYaxis()->FindBin(ptBinsEdges[i]), hMassPrompt->GetYaxis()->FindBin(ptBinsEdges[i+1])-1));
        vMassPtNonPrompt.push_back((TH1F *)hMassNonPrompt->ProjectionX(nameMassNonPrompt.c_str(), hMassNonPrompt->GetYaxis()->FindBin(ptBinsEdges[i]), hMassNonPrompt->GetYaxis()->FindBin(ptBinsEdges[i+1])-1));
        vMassPtRefl.push_back((TH1F *)hMassRefl->ProjectionX(nameReflMass.c_str(), hMassRefl->GetYaxis()->FindBin(ptBinsEdges[i]), hMassRefl->GetYaxis()->FindBin(ptBinsEdges[i+1])-1));

        if (verbose) {
            std::cout << "                                   " << hMassNonPrompt->GetYaxis()->FindBin(ptBinsEdges[i]) << " " << hMassNonPrompt->GetYaxis()->FindBin(ptBinsEdges[i+1])-1 << std::endl;
        }
    }

    // Pintamos las proyecciones
    for (int i=0; i < nPtBins; i++) {
        auto projName = "InvMassPt"+std::to_string(i);
        DrawProjection(vMassPtPrompt[i], ptBinsEdges[i], ptBinsEdges[i+1], projName.c_str());
        DrawProjection(vMassPtNonPrompt[i], ptBinsEdges[i], ptBinsEdges[i+1], projName.c_str());
        DrawProjection(vMassPtRefl[i], ptBinsEdges[i], ptBinsEdges[i+1], projName.c_str());
    }
}

void DrawProjection(TH1F* hD, int const& ptMin, int const& ptMax, TString filename) {
    TCanvas *cD = new TCanvas();
    cD->SetTickx();
    cD->SetTicky();
    cD->SetGrid();

    auto titleD = TString::Format("Invariant Mass, pT =(%d, %d) D^{0}; Mass (GeV); Entries", ptMin, ptMax);
    hD->SetMarkerStyle(7);
    hD->SetMarkerColor(kBlue);
    hD->SetTitle(titleD);
    hD->GetXaxis()->SetRangeUser(1.65, 2.10);
    hD->Draw("PE");

    hD->Write("" ,TObject::kOverwrite);

    if (cD)
    {
        cD->Close();
        gSystem->ProcessEvents();
        cD->Clear();
    }
}

void DrawProjectionSum(TH1F* hD, TH1F* hDbar, int const& ptMin, int const& ptMax, TString filename) {
    TCanvas *cD = new TCanvas();
    cD->SetTickx();
    cD->SetTicky();
    cD->SetGrid();

    auto titleD = TString::Format("Invariant Mass, pT =(%d, %d) D^{0} + D^{0}bar; Mass (GeV); Entries", ptMin, ptMax);
    hD->Add(hDbar);
    hD->SetMarkerStyle(7);
    hD->SetMarkerColor(kBlue);
    hD->SetTitle(titleD);
    hD->GetXaxis()->SetRangeUser(1.65, 2.10);
    hD->Draw("PE");

    hD->Write("" ,TObject::kOverwrite);

    if (cD)
    {
        cD->Close();
        gSystem->ProcessEvents();
        cD->Clear();
    }
}

/*     auto title = TString::Format("Invariant Mass, pT =(%d, %d) D^{0} + D^{0}bar; Mass (GeV); Entries", ptMin, ptMax);
    hD->Add(hDbar);
    hD->SetMarkerStyle(7);
    hD->SetMarkerColor(kBlue);
    hD->SetTitle(title);
    hD->GetXaxis()->SetRangeUser(1.65, 2.10);
    hD->Draw("PE");

    hD->Write("",TObject::kOverwrite);
    auto savePdf = filename+".pdf";
    //c->SaveAs(savePdf);
    if (c)
    {
        c->Close();
        gSystem->ProcessEvents();
        c->Clear();
    } */
