#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include "TRandom3.h"

using namespace std;

void smearElectron(const TString YEAR_STR, const TString type, const TString suff = "")
{

    //
    //  OPTIONS
    //

    const unsigned N = 500;     // Number of histograms to create
    // type = "id", "reco"
    // suff = "gap", "lowEt"



    //
    //  INPUT
    //

    TString inName;

    if (suff.IsNull())
        inName = "../data/electron_" + type + "_sf_" + YEAR_STR + ".root";
    else
        inName = "../data/electron_" + type + "_sf_" + YEAR_STR + "_" + suff + ".root";
//  inName = "../data/electron_" + type + "_sf_" + YEAR_STR + "_full.root";

    TFile   *inFile = TFile::Open(inName);

    TH2     *h_sf;
    inFile->GetObject("EGamma_SF2D", h_sf);
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

    TString outName;
    if (suff.IsNull())
        outName = "../data/electron_" + type + "_smear_" + YEAR_STR + ".root";
    else
        outName = "../data/electron_" + type + "_smear_" + YEAR_STR + "_" + suff + ".root";

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
