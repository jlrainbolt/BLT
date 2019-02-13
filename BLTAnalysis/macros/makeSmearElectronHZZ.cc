#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TRandom3.h"

using namespace std;

void makeSmearElectronHZZ()
{
    // HZZ ID SF FOR NON-GAP ELECTRONS



    // Number of histograms to create

    const int N = 100;


    // Get histograms
 
    TString inPath = "/uscms/home/jrainbol/nobackup/CMSSW_9_4_9_cand2/src/BLT/BLTAnalysis/data/electron_id/";

    TFile *inFile = TFile::Open(inPath + "egammaEffi.txt_EGM2D_Moriond2018v1.root");
    TH2 *h_sf;
    inFile->GetObject("EGamma_SF2D", h_sf);

//  h_sf->Draw("COL");



    // Retrieve pt, eta binning

    Int_t x_bins, y_bins;
    x_bins = h_sf->GetNbinsX();
    y_bins = h_sf->GetNbinsY();

    Double_t x_edge[20], y_edge[20];
    h_sf->GetXaxis()->GetLowEdge(x_edge);
    h_sf->GetYaxis()->GetLowEdge(y_edge);

    // Have to do the overflow bin manually because ROOT sucks
    x_edge[x_bins] = h_sf->GetXaxis()->GetBinLowEdge(x_bins+1);
    y_edge[y_bins] = h_sf->GetYaxis()->GetBinLowEdge(y_bins+1);



    // Print binning

    cout << "eta: ";
    for (unsigned i = 0; i < x_bins+1; i++)
        cout << x_edge[i] << ", ";
    cout << endl;

    cout << "pt: ";
    for (unsigned j = 0; j < y_bins+1; j++)
        cout << y_edge[j] << ", ";
    cout << endl;



    // Seed RNG (don't change!!)

    TRandom3 rng(11);



    // Create and fill histogram

    TString histName;
    TH2D *h_smr[N];

    TH2D *h_err = new TH2D("ERROR", "ERROR", x_bins, x_edge, y_bins, y_edge);
    h_err->GetXaxis()->SetTitle(h_sf->GetXaxis()->GetTitle());
    h_err->GetYaxis()->SetTitle(h_sf->GetYaxis()->GetTitle());


    for (unsigned n = 0; n < N; n++)
    {
        histName.Form("SMEAR%i", n);
        h_smr[n] = new TH2D(histName, histName, x_bins, x_edge, y_bins, y_edge);
        h_smr[n]->SetDirectory(0);

        for (unsigned i = 1; i <= x_bins; i++)
        {
            for (unsigned j = 1; j <= y_bins; j++)
            {
                Double_t err = h_sf->GetBinError(i, j);
                Double_t smr = rng.Gaus(0, err);

                h_smr[n]->SetBinContent(i, j, smr);

                if (n == 0)
                    h_err->SetBinContent(i, j, err);
            }
        }
    }
/*
    h_err->SetStats(0);
    h_err->Draw("COLZ");
    h_smr[0]->GetXaxis()->SetTitle(h_sf->GetXaxis()->GetTitle());
    h_smr[0]->GetYaxis()->SetTitle(h_sf->GetYaxis()->GetTitle());
    h_smr[0]->SetStats(0);
    h_smr[0]->Draw("COLZ");
*/


    // Write to file

    inFile->Close();

    TFile *outFile = new TFile("hzz_electron_id_smear.root", "RECREATE");
    for (unsigned n = 0; n < N; n++)
        h_smr[n]->Write();
    outFile->Close();
}
