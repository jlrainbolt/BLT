#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TRandom3.h"

using namespace std;

void makeSmearMuonHZZ()
{
    // Number of histograms to create

    const int N = 100;


    // Get histograms
 
    TString inPath = "/uscms/home/jrainbol/nobackup/CMSSW_9_4_9_cand2/src/BLT/BLTAnalysis/data/muon_id/";
    TFile *inFile = TFile::Open(inPath + "ScaleFactors_mu_Moriond2018_final.root");
    TH2 *h_err;
    inFile->GetObject("ERROR", h_err);

    h_err->SetStats(0);
    h_err->Draw("COLZ");



    // Retrieve pt, eta binning

    Int_t x_bins, y_bins;
    x_bins = h_err->GetNbinsX();
    y_bins = h_err->GetNbinsY();

    Double_t x_edge[20], y_edge[20];
    h_err->GetXaxis()->GetLowEdge(x_edge);
    h_err->GetYaxis()->GetLowEdge(y_edge);

    // Have to do the overflow bin manually because ROOT sucks
    x_edge[x_bins] = h_err->GetXaxis()->GetBinLowEdge(x_bins+1);
    y_edge[y_bins] = h_err->GetYaxis()->GetBinLowEdge(y_bins+1);


/*
    // Print binning

    cout << "eta: ";
    for (unsigned i = 0; i < x_bins+1; i++)
        cout << x_edge[i] << ", ";
    cout << endl;

    cout << "pt: ";
    for (unsigned j = 0; j < y_bins+1; j++)
        cout << y_edge[j] << ", ";
    cout << endl;
*/


    // Seed RNG (don't change!!)

    TRandom3 rng(13);



    // Create and fill histogram

    TString histName;
    TH2D *h_smr[N];

    for (unsigned n = 0; n < N; n++)
    {
        histName.Form("SMEAR%i", n);
        h_smr[n] = new TH2D(histName, histName, x_bins, x_edge, y_bins, y_edge);
        h_smr[n]->SetDirectory(0);

        for (unsigned i = 1; i <= x_bins; i++)
        {
            for (unsigned j = 1; j <= y_bins; j++)
            {
                Double_t err = h_err->GetBinContent(i, j);
                Double_t smr = rng.Gaus(0, err);

                h_smr[n]->SetBinContent(i, j, smr);
            }
        }
    }
/*
    h_smr[0]->GetXaxis()->SetTitle(h_err->GetXaxis()->GetTitle());
    h_smr[0]->GetYaxis()->SetTitle(h_err->GetYaxis()->GetTitle());
    h_smr[0]->SetStats(0);
    h_smr[0]->Draw("COLZ");
*/


    // Write to file

    inFile->Close();

    TFile *outFile = new TFile("hzz_muon_id_smear.root", "RECREATE");
    for (unsigned n = 0; n < N; n++)
        h_smr[n]->Write();
    outFile->Close();
}
