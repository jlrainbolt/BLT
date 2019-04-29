#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include "TRandom3.h"

using namespace std;

void smearMuon(const TString YEAR_STR)
{

    //
    //  OPTIONS
    //

    const unsigned N = 100;     // Number of histograms to create



    //
    //  INPUT
    //

    TString inName  = "../data/muon_id_sf_" + YEAR_STR + ".root";
    TFile   *inFile = TFile::Open(inName);

    TH2     *h_err, *h_sf;
    inFile->GetObject("ERROR", h_err);
    inFile->GetObject("FINAL", h_sf);

    h_err->SetDirectory(0);
    h_sf->SetDirectory(0);

    inFile->Close();
    cout << "Got histograms from " << inName << endl;


    //
    //  CREATE HISTOGRAMS
    //

    TString outName = "../data/muon_id_smear_" + YEAR_STR + ".root";
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

        for (unsigned i = 1; i <= h_sf->GetNbinsX(); i++)
        {
            for (unsigned j = 1; j <= h_sf->GetNbinsY(); j++)
            {
                float sf  = h_sf->GetBinContent(i, j);
                float err = h_err->GetBinContent(i, j);
                float smr = fabs(rng.Gaus(0, err));

                float var = n < 50 ? sf + smr : sf - smr;

                h_smr[n]->SetBinContent(i, j, var);
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
