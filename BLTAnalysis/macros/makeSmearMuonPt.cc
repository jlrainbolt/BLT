#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TRandom3.h"
#include "../interface/RoccoR.h"
#include "TMath.h"

using namespace std;

void makeSmearMuonPt()
{
    gROOT->ProcessLine(".L ../src/RoccoR.cc");

    // Number of histograms to create

    const int N = 100;


    // Choose pt, eta binning
    //      (eta is defined in RoccoR2017v1.txt; pt is chosen to match HZZ + a couple more)

    Int_t x_bins = 14, y_bins = 18;
    Double_t x_edge[] = {    -2.40, -2.10, -1.85, -1.60, -1.20, -0.80, -0.40,
                            0.00,  0.40,  0.80,  1.20,  1.60,  1.85,  2.10,  2.40};
    Double_t y_edge[] = {3, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100, 150, 200};



    // Declare dumb RoccoR object

    RoccoR rc("../data/RoccoR2017v1.txt");



    // Seed RNG (change?)

    TRandom3 rng(13);



    // Create and fill histogram

    TH2D *h_err = new TH2D("ERROR", "ERROR", x_bins, x_edge, y_bins, y_edge);
    h_err->SetDirectory(0);

    TString histName;
    TH2D *h_smr[N];

    for (unsigned n = 0; n < N; n++)
    {
        histName.Form("SMEAR%i", n);
        h_smr[n] = new TH2D(histName, histName, x_bins, x_edge, y_bins, y_edge);
        h_smr[n]->SetDirectory(0);

        for (unsigned i = 1; i <= x_bins; i++)
        {
            // Choose eta as bin center
            Double_t eta = h_smr[n]->GetXaxis()->GetBinCenter(i);

            for (unsigned j = 1; j <= y_bins; j++)
            {
                // Choose pt as bin center
                Double_t pt = h_smr[n]->GetYaxis()->GetBinCenter(j);

                // Randomize phi (maybe not the best idea?)
                Double_t phi = rng.Uniform(-TMath::Pi(), TMath::Pi());
                Int_t q = copysign(1, rng.Uniform(-1, 1));
                Int_t nl = rng.Gaus(13., 3.);
                if (nl < 9)
                    nl = 9;
                if (nl > 18)
                    nl = 18;
//              cout << nl << endl;

                Double_t errDT = rc.kScaleDTerror(q, pt, eta, phi);
                Double_t errMC = rc.kSmearMCerror(q, pt, eta, phi, nl, rng.Uniform(0, 1));
                Double_t err = sqrt(errDT * errDT + errMC * errMC);

//              Double_t smr = rng.Gaus(0, errMC);

                h_smr[n]->SetBinContent(i, j, err);

                if (n == 0)
                    h_err->SetBinContent(i, j, err);
            }
        }
    }

//  h_smr[0]->SetStats(0);
//  h_smr[0]->Draw("COLZ");
    h_err->SetStats(0);
    h_err->Draw("COLZ");



    // Write to file

    TFile *outFile = new TFile("muon_pt_smear.root", "RECREATE");
    for (unsigned n = 0; n < N; n++)
        h_smr[n]->Write();
    outFile->Close();
}
