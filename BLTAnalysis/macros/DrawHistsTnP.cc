#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"

#include "BLT/BLTAnalysis/bin/TagAndProbeAnalyzer.h"

using namespace std;
using namespace baconhep;



void DrawHistsTnP(const TString selection, const TString tag, const TString suffix)
{
    const unsigned  N_PT = selection.Contains("mumu") ? N_PT_MM : N_PT_EE;
    const unsigned  N_ETA = selection.Contains("mumu") ? N_ETA_MM : N_ETA_EE;

    float   *binsPt = selection.Contains("mumu") ? mumuPt : eePt;
    float   *binsEta = selection.Contains("mumu") ? mumuEta : eeEta;



    //
    //  INPUT FILES
    //

    const TString YEAR_STR = "2017";
    const TString EOS_PATH = "root://cmseos.fnal.gov//store/user/jrainbol/TagAndProbe/";

    TString inFileName  = selection + "_" + tag + "_" + suffix + ".root";
    TFile   *inFile = TFile::Open(EOS_PATH + YEAR_STR + "/" + inFileName);
    cout << endl << endl << "Opened " << inFileName << endl << endl;

    const TString inTreeName = "tree_" + suffix;
    TTreeReader reader(inTreeName, inFile);

    TTreeReaderValue    <Float_t>               genWeight_      (reader,    "genWeight");
    TTreeReaderValue    <Float_t>               ecalWeight_     (reader,    "ECALWeight");
    TTreeReaderValue    <Float_t>               puWeight_       (reader,    "PUWeight");
    TTreeReaderValue    <TLorentzVector>        tagP4_          (reader,    "tagP4");
    TTreeReaderValue    <TLorentzVector>        probeP4_        (reader,    "probeP4");
    TTreeReaderValue    <Bool_t>                probeIsTight_   (reader,    "probeIsTight");
    TTreeReaderValue    <Bool_t>                probeIsIso_     (reader,    "probeIsIsolated");
    TTreeReaderValue    <Float_t>               probeIso_       (reader,    "probeIsolation");



    //
    //  OUTPUT FILE
    //

    TString outFileName = selection + "_" + tag + "_" + suffix + ".root";
    TFile *outFile = new TFile(outFileName, "RECREATE");
    outFile->cd();

//  TDirectory *histDir = (TDirectory*) inFile->Get("hists_" + suffix);
//  histDir->SetName("mll_" + suffix);
//  histDir->Write();



    //
    // HISTOGRAMS
    //

    const int   N_QT = 60,      N_DPHI = 60,        N_ISO = 60;
    const float QT_MIN = 0,     DPHI_MIN = 0,       ISO_MIN = 0;
    const float QT_MAX = 300,   DPHI_MAX = M_PI,    ISO_MAX = tag.Contains("35") ? 0.35 : 1.4;

    // 2d
    outFile->cd();
    TH2 *passed2d, *failed2d;
    inFile->GetObject("hists_" + suffix + "/passed_2d", passed2d);
    inFile->GetObject("hists_" + suffix + "/failed_2d", failed2d);
    passed2d->SetDirectory(outFile);    failed2d->SetDirectory(outFile);

    // Qt
    outFile->mkdir("qt_" + suffix);
    outFile->cd("qt_" + suffix);

    TH1F *passedAllQt = new TH1F("passed_all", "Passed", N_QT, QT_MIN, QT_MAX);
    TH1F *failedAllQt = new TH1F("failed_all", "Failed", N_QT, QT_MIN, QT_MAX);
    passedAllQt->SetMinimum(0);     failedAllQt->SetMinimum(0);
    vector<vector<TH1F*>>   passedBinQt,    failedBinQt;

    // Delta phi (absolute value)
    outFile->mkdir("dphi_" + suffix);
    outFile->cd("dphi_" + suffix);

    TH1F *passedAllDPhi = new TH1F("passed_all", "Passed", N_DPHI, DPHI_MIN, DPHI_MAX);
    TH1F *failedAllDPhi = new TH1F("failed_all", "Failed", N_DPHI, DPHI_MIN, DPHI_MAX);
    passedAllDPhi->SetMinimum(0);   failedAllDPhi->SetMinimum(0);
    vector<vector<TH1F*>>   passedBinDPhi,  failedBinDPhi;

    // Relative isolation (excluding 0)
    outFile->mkdir("iso_" + suffix);
    outFile->cd("iso_" + suffix);

    TH1F *passedAllIso = new TH1F("passed_all", "Passed", N_ISO, ISO_MIN, ISO_MAX);
    TH1F *failedAllIso = new TH1F("failed_all", "Failed", N_ISO, ISO_MIN, ISO_MAX);
    passedAllIso->SetMinimum(0);    failedAllIso->SetMinimum(0);
    vector<vector<TH1F*>>   passedBinIso,   failedBinIso;

    // Loop to construct histograms
    for (unsigned i = 0; i <= N_PT; i++)
    {
        TString ptString;

        if (i < N_PT)
            ptString.Form("(%g,%g)", binsPt[i], binsPt[i+1]);
        else
            ptString.Form("(%g+)", binsPt[i]);

        vector<TH1F*> passedVecQt, failedVecQt, passedVecDPhi, failedVecDPhi, passedVecIso, failedVecIso;

        for (unsigned j = 0; j <= N_ETA; j++)
        {
            TString etaString;

            if (j < N_ETA)
                etaString.Form("(%g,%g)", binsEta[j], binsEta[j+1]);
            else
                etaString.Form("(%g,%g)", binsEta[0], binsEta[j]);

            TString name = "_pt" + ptString + "_eta" + etaString;
            TString title = ": pT " + ptString + ", eta " + etaString;

            outFile->cd("qt_" + suffix);
            passedVecQt.push_back(new TH1F("passed" + name, "Passed" + title, N_QT, QT_MIN, QT_MAX));
            failedVecQt.push_back(new TH1F("failed" + name, "Failed" + title, N_QT, QT_MIN, QT_MAX));

            outFile->cd("dphi_" + suffix);
            passedVecDPhi.push_back(new TH1F("passed" + name, "Passed" + title, N_DPHI, DPHI_MIN, DPHI_MAX));
            failedVecDPhi.push_back(new TH1F("failed" + name, "Failed" + title, N_DPHI, DPHI_MIN, DPHI_MAX));

            outFile->cd("iso_" + suffix);
            passedVecIso.push_back(new TH1F("passed" + name, "Passed" + title, N_ISO, ISO_MIN, ISO_MAX));
            failedVecIso.push_back(new TH1F("failed" + name, "Failed" + title, N_ISO, ISO_MIN, ISO_MAX));
        }
        passedBinQt.push_back(passedVecQt);         failedBinQt.push_back(failedVecQt);
        passedBinDPhi.push_back(passedVecDPhi);     failedBinDPhi.push_back(failedVecDPhi);
        passedBinIso.push_back(passedVecIso);       failedBinIso.push_back(failedVecIso);
    }


    //
    //  EVENT LOOP
    //

    int nEvents = reader.GetEntries(kTRUE);

    cout << "Running over " << nEvents << " total events" << endl;

    while (reader.Next())// && reader.GetCurrentEntry() < 1000000)
    {
        if (reader.GetCurrentEntry() % 10000 == 0)
            cout << "Processed " << reader.GetCurrentEntry() << " of " << nEvents << " events" << endl; 

        TLorentzVector  tagP4 = *tagP4_,    probeP4 = *probeP4_,    totalP4 = tagP4 + probeP4;
        bool    probeIsIso = *probeIsIso_,  probeIsTight = *probeIsTight_;
        float   probeIso = *probeIso_,      weight = (*genWeight_) * (*ecalWeight_) * (*puWeight_);

        unsigned I = passed2d->GetYaxis()->FindBin(probeP4.Pt()) - 1;
        unsigned J = passed2d->GetXaxis()->FindBin(probeP4.Eta()) - 1;

        if (probeIsTight)
        {
            passedAllQt->Fill(totalP4.Pt(), weight);
            passedAllDPhi->Fill(fabs(tagP4.DeltaPhi(probeP4)), weight);
            if (probeIso > 0)
                passedAllIso->Fill(probeIso / probeP4.Pt(), weight);

            if (I <= N_PT)
            {
                passedBinQt[I][N_ETA]->Fill(totalP4.Pt(), weight);
                passedBinDPhi[I][N_ETA]->Fill(fabs(tagP4.DeltaPhi(probeP4)), weight);
                if (probeIso > 0)
                    passedBinIso[I][N_ETA]->Fill(probeIso / probeP4.Pt(), weight);

                if (J<= N_ETA)
                {
                    passedBinQt[I][J]->Fill(totalP4.Pt(), weight);
                    passedBinDPhi[I][J]->Fill(fabs(tagP4.DeltaPhi(probeP4)), weight);
                    if (probeIso > 0)
                        passedBinIso[I][J]->Fill(probeIso / probeP4.Pt(), weight);
                }
            }
        }
        else
        {
            failedAllQt->Fill(totalP4.Pt(), weight);
            failedAllDPhi->Fill(fabs(tagP4.DeltaPhi(probeP4)), weight);
            if (probeIso > 0)
                failedAllIso->Fill(probeIso / probeP4.Pt(), weight);

            if (I <= N_PT)
            {
                failedBinQt[I][N_ETA]->Fill(totalP4.Pt(), weight);
                failedBinDPhi[I][N_ETA]->Fill(fabs(tagP4.DeltaPhi(probeP4)), weight);
                if (probeIso > 0)
                    failedBinIso[I][N_ETA]->Fill(probeIso / probeP4.Pt(), weight);

                if (J<= N_ETA)
                {
                    failedBinQt[I][J]->Fill(totalP4.Pt(), weight);
                    failedBinDPhi[I][J]->Fill(fabs(tagP4.DeltaPhi(probeP4)), weight);
                    if (probeIso > 0)
                        failedBinIso[I][J]->Fill(probeIso / probeP4.Pt(), weight);
                }
            }
        }
    }

    inFile->Close();
    cout << "Closed " << inFileName << endl;

    outFile->Write();
    outFile->Close();
    cout << "Closed " << outFileName << endl;
}
