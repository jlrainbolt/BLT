#include <iostream>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"

#include "BaconAna/DataFormats/interface/TEventInfo.hh"
#include "BaconAna/DataFormats/interface/TGenParticle.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"

using namespace std;
using namespace baconhep;



void printTruthInfo(const TString inFileName, const TString inTreeName = "Events")
{
    // Stuff copied from BLTSelector.h
    baconhep::TEventInfo    *fInfo              = 0;
    TClonesArray            *fGenParticleArr    = 0;
    TClonesArray            *fElectronArr       = 0;
//  TClonesArray            *fMuonArr           = 0;

    TBranch                 *b_Info;
    TBranch                 *b_GenParticleArr;
    TBranch                 *b_ElectronArr;
//  TBranch                 *b_MuonArr;


    // Open ROOT file and get branches
    TFile *inFile = TFile::Open(inFileName);
    cout << "Opened " << inFileName << endl;
    TTree *inTree;
    inFile->GetObject(inTreeName, inTree);
    cout << "Got tree " << inTreeName << endl;

    inTree->SetBranchAddress("Info", &fInfo, &b_Info);
    inTree->SetBranchAddress("GenParticle", &fGenParticleArr, &b_GenParticleArr);
    inTree->SetBranchAddress("Electron", &fElectronArr, &b_ElectronArr);
//  inTree->SetBranchAddress("Muon", &fMuonArr, &b_MuonArr);

    unsigned nEvents = inTree->GetEntries();


    // Loop over events
    cout << "Reading " << nEvents << " events..." << endl << endl << endl;

    unsigned i = 0;
    while (inTree->GetEntry(i))
    {
        cout << endl;
        cout << "Run number:\t" << fInfo->runNum << endl;
        cout << "Lumi section:\t" << fInfo->lumiSec << endl;
        cout << "Event number:\t" << fInfo->evtNum << endl;
        cout << endl;
        cout << "Index" << "\t" << "Parent" << "\t" << "ID" << "\t" << "Status" << "\t";
        cout << "Pt" << "\t\t" << "Match";
        cout << endl;

        for (unsigned j = 0; j < fGenParticleArr->GetEntries(); j++)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(j);
            cout << j << "\t" << particle->parent << "\t" << particle->pdgId << "\t";
            cout << particle->status;

            if (
                    (abs(particle->pdgId) == 11) ||
                    (abs(particle->pdgId) == 13) ||
                    (particle->pdgId == 22)
                )
                cout << "\t" << particle->pt;

            // Try to match photons
            if (
                    (particle->pdgId == 22) &&
                    (particle->pt > 5) &&
                    (particle->status == 1)
                )
            {
                TLorentzVector gammaP4, electronP4;
                copy_p4(particle, 0, gammaP4);

                for (unsigned k = 0; k < fElectronArr->GetEntries(); k++)
                {
                    TElectron* electron = (TElectron*) fElectronArr->At(k);
                    copy_p4(electron, ELE_MASS, electronP4);

                    if (gammaP4.DeltaR(electronP4) < 0.02)
                        cout << "\t\t" << electron->q * (-11);
                }
            }

            cout << endl;
        }
        cout << endl << endl;

        i++;
    }

    inFile->Close();
    cout << "Closed " << inFileName << endl;
}
