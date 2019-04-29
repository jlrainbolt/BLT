// =============================================================================
// A simple analysis on Bacon ntuples
//
// Input arguments:
//   argv[1] => input bacon file name
//   argv[2] => number of entries
//   argv[3] => ...
//
// Users should inherit from BLTSelector and implement the three functions:
//   Begin()
//   Process()
//   Terminate()
// =============================================================================


#ifndef MULTILEPTONANALYZER_HH
#define MULTILEPTONANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"

// BaconAna class definitions (might need to add more)
#include "BaconAna/Utils/interface/TTrigger.hh"

// ROOT headers
#include <TString.h>
#include <TLorentzVector.h>

// C++ headers
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <regex>


class MultileptonAnalyzer: public BLTSelector {
public:
    MultileptonAnalyzer();
    ~MultileptonAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    TTree *outTree;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    std::vector<string>                 muonTriggerNames, electronTriggerNames;



    //
    //  BRANCHES
    //

    // Event
    Int_t       runNumber,          lumiSection,            evtNumber;
    Bool_t      evtMuonTriggered,   evtElectronTriggered;
    Float_t     genWeight,          PUWeight,               nPU;
    UShort_t    nPV;
    Float_t     ECALWeight,         ECALWeightUp,       ECALWeightDown;
    Bool_t      hasTauDecay;

    // Counters
    UShort_t    nLooseMuons,        nLooseElectrons,        nLooseLeptons;
    UShort_t    nTightMuons,        nTightElectrons,        nTightLeptons;
    
    // Muons
    TClonesArray            *muonP4_                = new TClonesArray("TLorentzVector");
    TClonesArray            *muonUncorrectedP4_     = new TClonesArray("TLorentzVector");
    std::vector<Short_t>    muonCharge;
    std::vector<Float_t>    muonEnergySF,           muonEnergySFUp,         muonEnergySFDown;
    std::vector<Float_t>    muonIDSF,               muonIsolation;
    std::vector<Bool_t>     muonIsTight,            muonIsLoose,            muonIsIsolated;
    std::vector<Bool_t>     muonIsPF,               muonIsTrackerHighPt; 
    std::vector<Bool_t>     muonFiredLeg1,          muonFiredLeg2;
    std::vector<Float_t>    muonTrigEffLeg1MC,      muonTrigEffLeg1Data;
    std::vector<Float_t>    muonTrigEffLeg2MC,      muonTrigEffLeg2Data;

    // Electrons
    TClonesArray            *electronP4_            = new TClonesArray("TLorentzVector");
    TClonesArray            *electronUncorrectedP4_ = new TClonesArray("TLorentzVector");
    std::vector<Short_t>    electronCharge;
    std::vector<Float_t>    electronEnergySF,       electronEnergySFUp,     electronEnergySFDown;
    std::vector<Float_t>    electronIDSF,           electronRecoSF;
    std::vector<Float_t>    electronIsolation,      electronScEta;
    std::vector<Bool_t>     electronIsTight,        electronIsLoose,        electronIsV2Iso;
    std::vector<Bool_t>     electronIsGap;
    std::vector<Bool_t>     electronFiredLeg1,      electronFiredLeg2;
    std::vector<Float_t>    electronTrigEffLeg1MC,  electronTrigEffLeg1Data;
    std::vector<Float_t>    electronTrigEffLeg2MC,  electronTrigEffLeg2Data;


    // Gen particles

    // Counters
    UShort_t    nFinalStateMuons,   nFinalStateElectrons,   nFinalStateLeptons;
    UShort_t    nHardProcMuons,     nHardProcElectrons,     nHardProcLeptons;

    // Final state leptons
    TLorentzVector          finalStateLeptonsP4;
    TClonesArray            *finalStateMuonP4_     = new TClonesArray("TLorentzVector");
    TClonesArray            *finalStateElectronP4_ = new TClonesArray("TLorentzVector");
    std::vector<Short_t>    finalStateMuonQ,       finalStateElectronQ;
    std::vector<UShort_t>   finalStateMuonZIndex,  finalStateElectronZIndex;

    // Hard process leptons
    TLorentzVector          hardProcLeptonsP4;
    TClonesArray            *hardProcMuonP4_       = new TClonesArray("TLorentzVector");
    TClonesArray            *hardProcElectronP4_   = new TClonesArray("TLorentzVector");
    std::vector<Short_t>    hardProcMuonQ,         hardProcElectronQ;
    std::vector<UShort_t>   hardProcMuonZIndex,    hardProcElectronZIndex;


    //ClassDef(MultileptonAnalyzer,0);
};


#endif  // MULTILEPTONANALYZER_HH
