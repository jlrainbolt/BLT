#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"

using namespace std;
using namespace baconhep;



void printTruthInfo(const TString inFileName)
{

    // Stuff copied from BLTSelector.h
    baconhep::TEventInfo    *fInfo              = 0;
    TClonesArray            *fGenParticleArr    = 0;
//  TClonesArray            *fElectronArr       = 0;
//  TClonesArray            *fMuonArr           = 0;

    TBranch                 *b_Info;
    TBranch                 *b_GenParticleArr;
//  TBranch                 *b_ElectronArr;
//  TBranch                 *b_MuonArr;



    //
    //  BACON NTUPLE
    //

    TFile *inFile = TFile::Open(inFileName);
    cout << "Opened " << inFileName << endl;
    
    const TString inTreeName = "Events";
    TTree *inTree;
    inFile->GetObject(inTreeName, inTree);
    cout << "Got tree " << inTreeName << endl;

    inTree->SetBranchAddress("Info", &fInfo, &b_Info);
    inTree->SetBranchAddress("GenParticle", &fGenParticleArr, &b_GenParticleArr);
//  inTree->SetBranchAddress("Electron", &fElectronArr, &b_ElectronArr);
//  inTree->SetBranchAddress("Muon", &fMuonArr, &b_MuonArr);

    unsigned nEvents = inTree->GetEntries();



    //
    //  SELECTED TREE
    //

    const TString suffix = "zjets_m-50",    YEAR_STR = "2016";
    const TString EOS_PATH = "root://cmseos.fnal.gov//store/user/jrainbol";

    TString inName  = "background_" + suffix + ".root";
    TString inPath  = EOS_PATH + "/Selected/" + YEAR_STR + "/" + inName;
    TFile   *inFile2 = TFile::Open(inPath);

    cout << endl << endl << "Opened " << inPath << endl;

    const TString selection = "4l";
    TTreeReader reader(selection + "_" + suffix, inFile2);

    TTreeReaderValue    <Int_t>                 runNum_         (reader,    "runNum");
    TTreeReaderValue    <Int_t>                 evtNum_         (reader,    "evtNum");
    TTreeReaderValue    <Int_t>                 lumiSec_        (reader,    "lumiSec");
    TTreeReaderValue    <TLorentzVector>        zzp4_           (reader,    "zzp4");
    TTreeReaderValue    <TLorentzVector>        l1p4_           (reader,    "l1p4");
    TTreeReaderValue    <Short_t>               l1pdg_          (reader,    "l1pdg");
    TTreeReaderValue    <TLorentzVector>        l2p4_           (reader,    "l2p4");
    TTreeReaderValue    <Short_t>               l2pdg_          (reader,    "l2pdg");
    TTreeReaderValue    <TLorentzVector>        l3p4_           (reader,    "l3p4");
    TTreeReaderValue    <Short_t>               l3pdg_          (reader,    "l3pdg");
    TTreeReaderValue    <TLorentzVector>        l4p4_           (reader,    "l4p4");
    TTreeReaderValue    <Short_t>               l4pdg_          (reader,    "l4pdg");



    //
    //  EVENT LOOP
    //

    cout << "Reading " << nEvents << " events..." << endl << endl << endl;

    unsigned i = 0;
    while (inTree->GetEntry(i))
    {
        cout << endl;
        cout << "Run number:\t" << fInfo->runNum << endl;
        cout << "Lumi section:\t" << fInfo->lumiSec << endl;
        cout << "Event number:\t" << fInfo->evtNum << endl;
        cout << endl;

        TLorentzVector          zzP4;
        vector<TLorentzVector>  lepP4(4);
        vector<short>           lepPDG(4);


        // Find entry in selected tree
        while (reader.Next())
        {
            if (*runNum_ != fInfo->runNum)
                continue;
            if (*lumiSec_ != fInfo->lumiSec)
                continue;
            if (*evtNum_ != fInfo->evtNum)
                continue;
            // else, event is a match!
            
            zzP4 = *zzp4_;
            lepP4[0] = *l1p4_;      lepPDG[0] = *l1pdg_;
            lepP4[1] = *l2p4_;      lepPDG[1] = *l2pdg_;
            lepP4[2] = *l3p4_;      lepPDG[2] = *l3pdg_;
            lepP4[3] = *l4p4_;      lepPDG[3] = *l4pdg_;
        }

        cout << "4l mass:\t" << zzP4.M() << endl;
        cout << endl;
        cout << "RECO" << endl;
        cout << "Rank" << "\t" << "ID" << "\t" << "Pt" << "\t\t";
        cout << "Match DR" << "\t\t" << "Index" << endl;



        //
        //  MATCHING
        //

        vector<unsigned> matches(fGenParticleArr->GetEntries(), 99);

        for (unsigned k = 0; k < lepP4.size(); k++)  // loop over reco muons
        {
            // Find DeltaR between reco lep and each gen lep
            vector<float> deltaR(fGenParticleArr->GetEntries());

            for (unsigned j = 0; j < fGenParticleArr->GetEntries(); j++)
            {
                TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(j);
                TLorentzVector genP4;
                copy_p4(particle, genP4);

                if (particle->status == 1)
                    deltaR[j] = lepP4[k].DeltaR(genP4);
                else
                    deltaR[j] = 99;
            }

            // Find index of minimum DeltaR
            unsigned m = min_element(deltaR.begin(), deltaR.end()) - deltaR.begin();
            if (matches[m] > 4)
                matches[m] = k;

            // Print
            cout << k + 1 << "\t" << lepPDG[k] << "\t" << lepP4[k].Pt() << "\t\t";
            cout << deltaR[m] << "\t" << m << endl;
        }
        cout << endl << endl;



        //
        //  GEN INFO
        //

        cout << "GEN" << endl;
        cout << "Index" << "\t" << "Parent" << "\t" << "ID" << "\t" << "Status" << "\t";
        cout << "Pt" << "\t\t" << "Match";
        cout << endl;

        for (unsigned j = 0; j < fGenParticleArr->GetEntries(); j++)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(j);
            cout << j << "\t" << particle->parent << "\t" << particle->pdgId << "\t";
            cout << particle->status;

            if      (matches[j] < 99)
                cout << "\t" << particle->pt << "\t\t" << matches[j] + 1;
            else if (particle->status == 1)
                cout << "\t" << particle->pt;

            cout << endl;
        }
        cout << endl << endl;

        i++;
        reader.Restart();
    }

    inFile->Close();
    inFile2->Close();
    cout << "Closed " << inFileName << endl;
    cout << "Closed " << inName << endl;
}
