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


#ifndef TAGANDPROBEANALYZER_HH
#define TAGANDPROBEANALYZER_HH

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
#include <TRandom3.h>

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


float frac_diff(float test, float truth) { return fabs(test - truth) / truth; }



//
//  CUTS
//

const unsigned  N_MLL = 60;
const float     MLL_MIN = 60,       MLL_MAX = 120;
const float     MUON_PT_THRESH = 10,ELEC_PT_THRESH = 12;

const float     DR_MIN = 0.02;

const unsigned  N_PT_MM = 18,   N_ETA_MM = 4;
const unsigned  N_PT_EE = 17,   N_ETA_EE = 4;

float   mumuPt[N_PT_MM+1] = {5, 8,  10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 40, 45, 50, 60, 75, 100, 200 };
float   mumuEta[N_ETA_MM+1]={-2.4,  -1.2,       0,  1.2,    2.4 };

float   eePt[N_PT_EE+1] =   {7, 10, 12, 14, 16, 18, 20, 22, 25, 30, 35, 40, 45, 50, 60, 75, 100, 200 };
float   eeEta[N_ETA_EE+1] = {-2.5,  -1.479,     0,  1.479,  2.5 };


class TagAndProbeAnalyzer: public BLTSelector {
public:
    TagAndProbeAnalyzer();
    ~TagAndProbeAnalyzer();

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
    TString                             selection;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    // RNG
    TRandom3 rng;

    // Triggers
    std::vector<string>                 muonTriggerNames, electronTriggerNames;



    //
    //  BRANCHES
    //

    // Event
    Int_t       runNumber,          lumiSection,            evtNumber;
    Bool_t      evtMuonTriggered,   evtElectronTriggered;
    Float_t     genWeight,          PUWeight,               ECALWeight,             weight;
    UShort_t    nPV;
    UShort_t    nMuons,             nElectrons,             nLeptons;
    
    // Leptons
    TLorentzVector          tagP4,                  probeP4,                totalP4;
    Short_t                 tagQ,                   probeQ;
    Short_t                 tagPDG,                 probePDG;
    Float_t                 tagIsolation,           probeIsolation;
    Float_t                 tagEnergySF,            probeEnergySF;
    Bool_t                  tagFiredSingle,         probeFiredSingle;
    Bool_t                  tagIsLoose,             probeIsLoose;
    Bool_t                  tagIsTight,             probeIsTight;

    // Muons
    Bool_t                  tagIsIsolated,          probeIsIsolated;
    Bool_t                  tagIsPF,                probeIsPF;
    Bool_t                  tagIsTrackerHighPt,     probeIsTrackerHighPt;

    // Electrons
    Float_t                 tagScEta,               probeScEta;
    Bool_t                  tagIsGap,               probeIsGap;
    Bool_t                  tagIsV2Iso,             probeIsV2Iso;


    // Gen particles
    UShort_t                nGenMuons,              nGenElectrons,          nGenLeptons;
    TLorentzVector          genTagP4,               genProbeP4,             genZP4;
    Float_t                 genTagDeltaR,           genProbeDeltaR;
    Short_t                 genTagPDG,              genProbePDG;
    Short_t                 genTagMotherPDG,        genProbeMotherPDG;
    UShort_t                nGenTagPhotons,         nGenProbePhotons;
    TClonesArray            *genTagPhotonsP4        = new TClonesArray("TLorentzVector");
    TClonesArray            *genProbePhotonsP4      = new TClonesArray("TLorentzVector");



    //
    //  HISTOGRAMS
    //

    TH2F    *passed2d,      *failed2d;
    TH1F    *passedAll,     *failedAll,     *genPassedAll,  *genFailedAll;
    TH1F    *gapPassedAll,  *gapFailedAll;
    vector<vector<TH1F*>>   passedBin,      failedBin,      genPassedBin,   genFailedBin;
};


#endif  // TAGANDPROBEANALYZER_HH
