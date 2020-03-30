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


#ifndef EFFICIENCYANALYZER_HH
#define EFFICIENCYANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"

// BaconAna class definitions (might need to add more)
#include "BaconAna/DataFormats/interface/TLHEWeight.hh"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>

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
#include <utility>
#include <iterator>
#include <regex>


class EfficiencyAnalyzer: public BLTSelector {
public:
    EfficiencyAnalyzer();
    ~EfficiencyAnalyzer();

    void    Begin(TTree *tree);
//  void    Init(TTree *tree);
    Bool_t  Process(Long64_t entry);
    void    Terminate();
    void    ReportPostBegin();
    void    ReportPostTerminate();

    TFile *outFile;
    TTree *outTree;

    // Params and cuts
    std::unique_ptr<Parameters> params;
    std::unique_ptr<Cuts>               cuts;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<WeightUtils>        weights;

    // Histograms
    TH1D *hSelectedEvents,  *hMatchedEvents;



    //
    //  SELECTION CUTS
    //

    // Phase space requirements
    Float_t     M_MIN = 80,     M_MAX = 100,    MLL_MIN = 4;

    // Fiducial requirements
    Float_t     PT1_MIN = 20,   PT2_MIN = 10,   MU_PT_MIN = 5,  EL_PT_MIN = 7;
    Float_t     MU_ETA_MAX = 2.4,               EL_ETA_MAX = 2.5;

    // FSR recovery, matching
    Float_t     DR_DRESS = 0.1,                 MATCH_DR_MAX = 0.1;



    //
    //  BRANCHES
    //
    
    // Event
    UInt_t      runNumber,  lumiSection,    evtNumber;
    Float_t     genWeight,  PUWeight,       nPU,        nPV;
    UShort_t    decayChannel;
    Bool_t      isMatched;


    // Counters
    UShort_t    nLooseMuons,    nLooseElectrons,    nLooseLeptons;
    UShort_t    nTightMuons,    nTightElectrons,    nTightLeptons;
    UShort_t    nDressedMuons,  nDressedElectrons,  nDressedLeptons;


    // Reco leptons
    TClonesArray    *muonP4_        = new TClonesArray("TLorentzVector");
    std::vector<Float_t>    muonIsolation,          muonDR;
    std::vector<Bool_t>     muonIsTight,            muonIsLoose,            muonIsIsolated;
    std::vector<Bool_t>     muonIsMatched;
    std::vector<Short_t>    muonQ,                  muonMatchIdx;

    TClonesArray    *electronP4_    = new TClonesArray("TLorentzVector");
    std::vector<Float_t>    electronIsolation,      electronDR;
    std::vector<Bool_t>     electronIsTight,        electronIsLoose,        electronIsV2Iso;
    std::vector<Bool_t>     electronIsMatched;
    std::vector<Short_t>    electronQ,              electronMatchIdx;


    // Dressed leptons
    TClonesArray    *dressedMuonP4_     = new TClonesArray("TLorentzVector");
    TClonesArray    *dressedElectronP4_ = new TClonesArray("TLorentzVector");

    std::vector<Short_t>    dressedMuonQ,           dressedMuonMatchIdx;
    std::vector<Bool_t>     dressedMuonIsMatched,   dressedElectronIsMatched;
    std::vector<Short_t>    dressedElectronQ,       dressedElectronMatchIdx;


    //ClassDef(EfficiencyAnalyzer,0);
};


#endif  // EFFICIENCYANALYZER_HH
