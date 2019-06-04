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
    std::unique_ptr<Parameters> params;

    // Histograms
    TH1D *hPhaseSpaceEvents, *hFiducialEvents;



    //
    //  SELECTION CUTS
    //

    // Phase space requirements
    Float_t     M_MIN = 80,     M_MAX = 100,    MLL_MIN = 4;

    // Fiducial requirements
    Float_t     PT1_MIN = 20,   PT2_MIN = 10,   PT_MIN = 5,     ETA_MAX = 2.5;

    // FSR recovery
    Float_t     DR_DRESS = 0.1;



    //
    //  BRANCHES
    //
    
    // Event
    UInt_t      runNumber,  evtNumber,  lumiSection;
    Float_t     genWeight,  nomWeight;

    UShort_t    decayChannel;
    Bool_t      isFiducial;

    std::vector<UShort_t>   qcdID,      pdfID;
    std::vector<Float_t>    qcdWeight,  pdfWeight;


    // Counters
    UShort_t    nMuons,   nElectrons,   nLeptons,   nZs;


    // Dressed leptons
    TClonesArray    *muonP4_        = new TClonesArray("TLorentzVector");
    TClonesArray    *electronP4_    = new TClonesArray("TLorentzVector");
    TLorentzVector  leptonsP4;

    std::vector<Short_t>    muonQ,        electronQ;
    std::vector<Short_t>    muonMother,   electronMother;
    std::vector<UShort_t>   muonZIndex,   electronZIndex;


    // Z bosons
    TClonesArray    *zP4_             = new TClonesArray("TLorentzVector");
    TLorentzVector  zsP4;

    std::vector<Short_t>    zStatus;
    std::vector<UShort_t>   zIndex;


    //ClassDef(PhaseSpaceAnalyzer,0);
};


#endif  // PHASESPACEANALYZER_HH
