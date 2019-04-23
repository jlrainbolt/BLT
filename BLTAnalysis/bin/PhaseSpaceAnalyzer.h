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


#ifndef PHASESPACEANALYZER_HH
#define PHASESPACEANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"

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


class PhaseSpaceAnalyzer: public BLTSelector {
public:
    PhaseSpaceAnalyzer();
    ~PhaseSpaceAnalyzer();

    void    Begin(TTree *tree);
    void    Init(TTree *tree);
    Bool_t  Process(Long64_t entry);
    void    Terminate();
    void    ReportPostBegin();
    void    ReportPostTerminate();

    TFile *outFile;
    TTree *outTree;

    // Params and cuts
    std::unique_ptr<Parameters>         params;

    // Histograms
    TH1D *hPhaseSpaceEvents, *hFiducialEvents;



    //
    //  SELECTION CUTS
    //

    // Phase space requirements
    Float_t     M_MIN = 80,     M_MAX = 100,    MLL_MIN = 4;

    // Fiducial requirements
    Float_t     PT1_MIN = 20,   PT2_MIN = 10,   PT_MIN = 5,     ETA_MAX = 2.5;



    //
    //  BRANCHES
    //
    
    // Event
    Int_t       runNumber,  lumiSection;
    Long64_t    evtNumber;
    Float_t     genWeight,  nomWeight;

    UShort_t    decayChannel;
    Bool_t      isFiducial;

    std::vector<UShort_t>   qcdID,      pdfID;
    std::vector<Float_t>    qcdWeight,  pdfWeight;


    // Counters
    UShort_t    nStatus22Zs,        nFinalStateZs;

    UShort_t    nFinalStateMuons,   nFinalStateElectrons,   nFinalStateLeptons;     // Status 1 & traceable to a Z
    UShort_t    nHardProcMuons,     nHardProcElectrons,     nHardProcLeptons;       // Mother is a Z


    // Final state leptons
    TClonesArray    *finalStateMuonP4           = new TClonesArray("TLorentzVector");
    TClonesArray    &finalStateMuonP4ptr        = *finalStateMuonP4;

    TClonesArray    *finalStateElectronP4       = new TClonesArray("TLorentzVector");
    TClonesArray    &finalStateElectronP4ptr    = *finalStateElectronP4;

    std::vector<Short_t>    finalStateMuonQ,        finalStateElectronQ;
    std::vector<Short_t>    finalStateMuonMother,   finalStateElectronMother;
    std::vector<UShort_t>   finalStateMuonZIndex,   finalStateElectronZIndex;

    TLorentzVector  finalStateLeptonsP4;


    // Hard process leptons
    TClonesArray    *hardProcMuonP4             = new TClonesArray("TLorentzVector");
    TClonesArray    &hardProcMuonP4ptr          = *hardProcMuonP4;

    TClonesArray    *hardProcElectronP4         = new TClonesArray("TLorentzVector");
    TClonesArray    &hardProcElectronP4ptr      = *hardProcElectronP4;

    std::vector<Short_t>    hardProcMuonQ,          hardProcElectronQ;
    std::vector<Short_t>    hardProcMuonStatus,     hardProcElectronStatus;
    std::vector<UShort_t>   hardProcMuonZIndex,     hardProcElectronZIndex;

    TLorentzVector  hardProcLeptonsP4;


    // Z bosons
    TClonesArray    *status22ZP4                = new TClonesArray("TLorentzVector");
    TClonesArray    &status22ZP4ptr             = *status22ZP4;

    TClonesArray    *finalStateZP4              = new TClonesArray("TLorentzVector");
    TClonesArray    &finalStateZP4ptr           = *finalStateZP4;

    std::vector<Short_t>    status22ZMother,        finalStateZStatus;
    std::vector<UShort_t>   status22ZIndex,         finalStateZIndex;

    TLorentzVector  status22ZsP4,       finalStateZsP4;


    //ClassDef(PhaseSpaceAnalyzer,0);
};


#endif  // PHASESPACEANALYZER_HH
