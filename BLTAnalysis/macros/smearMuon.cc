#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include "TRandom3.h"

using namespace std;

void smearMuon()
{

    //
    //  OPTIONS
    //

    const unsigned N = 300;     // Number of histograms to create



    //
    //  INPUT
    //

    TString inName = "../data/MuonScaleFactors_2011_2012.root";

    TFile   *inFile = TFile::Open(inName);

    TH2     *h_sf;
    inFile->GetObject("TH2D_ALL_2012", h_sf);
    h_sf->SetDirectory(0);
    h_sf->SetNameTitle("FINAL", "FINAL");

    inFile->Close();
    cout << "Got histogram from " << inName << endl;



    //
    //  CREATE HISTOGRAMS
    //

    // Error histogram

    TH2F    *h_err = (TH2F*) h_sf->Clone("ERROR");
    h_err->Reset("ICESM");
    h_err->SetTitle("ERROR");

    for (unsigned i = 1; i <= h_sf->GetNbinsX(); i++)
    {
        for (unsigned j = 1; j <= h_sf->GetNbinsY(); j++)
        {
            float   err = h_sf->GetBinError(i, j);
            h_err->SetBinContent(i, j, err);
        }
    }


    // Smeared histograms

    TString outName = "../data/muon_id_smear_2012.root";

    TFile   *outFile = new TFile(outName, "RECREATE");

    TH2F        *h_smr[N];
    TString     histName;
    TRandom3    rng(0);

    cout << "Creating " << N << " smeared histograms..." << endl;
    for (unsigned n = 0; n < N; n++)
    {
        histName.Form("SMEAR%i", n);

        h_smr[n] = (TH2F*) h_sf->Clone(histName);
        h_smr[n]->Reset("ICESM");
        h_smr[n]->SetTitle(histName);

        float delta = rng.Gaus(0, 0.01);

        for (unsigned i = 1; i <= h_sf->GetNbinsX(); i++)
        {
            for (unsigned j = 1; j <= h_sf->GetNbinsY(); j++)
            {
                float sf  = h_sf->GetBinContent(i, j);
                float err = h_err->GetBinContent(i, j);
                float smr = rng.Gaus(0, err);

                h_smr[n]->SetBinContent(i, j, sf * (1 + delta) + smr);
            }
        }
    }



    //
    //  WRITE
    //

    outFile->cd();

    h_sf->Write();
    h_err->Write();
    for (unsigned n = 0; n < N; n++)
        h_smr[n]->Write();

    outFile->Close();
    cout << "Wrote histograms to " << outName << endl;
}
