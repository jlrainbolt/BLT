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


#ifndef HZGANALYZER_HH
#define HZGANALYZER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/TriggerSelector.hh"
#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"
#include "BLT/BLTAnalysis/interface/WeightUtils.h"

#include "BLT/BLTAnalysis/interface/RoccoR.h"

// BaconAna class definitions (might need to add more)
#include "BaconAna/Utils/interface/TTrigger.hh"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>

// C++ headers
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <memory>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <regex>

class hzgAnalyzer: public BLTSelector {
public:
    hzgAnalyzer();
    ~hzgAnalyzer();

    void   Begin(TTree *tree);
    Bool_t Process(Long64_t entry);
    void   Terminate();
    void   ReportPostBegin();
    void   ReportPostTerminate();

    TFile *outFile;
    TTree *outTree;

    // Lumi mask
    RunLumiRangeMap lumiMask;

    // rochester muon corrections
    RoccoR *muonCorr;
    TRandom3 *rng;

    // Params and cuts
    std::unique_ptr<Parameters>         params;
    std::unique_ptr<Cuts>               cuts;

    // Utilities
    std::unique_ptr<ParticleSelector>   particleSelector;
    std::unique_ptr<baconhep::TTrigger> trigger;
    std::unique_ptr<WeightUtils>        weights;

    std::vector<string> triggerNames;

    // Branches in the output file
    
    // event data
    UInt_t runNumber, lumiSection, nPV, nPartons;
    ULong64_t evtNumber;
    Bool_t triggerStatus;
    Float_t eventWeight, triggerWeight, puWeight, nPU;
    Int_t genWeight;
    Float_t xPV, yPV, zPV;
    UInt_t nJets, nCentralJets, nFwdJets, nBJets, nMuons, nElectrons, nTaus, nPhotons;

    // physics object Lorentz vectors
    Float_t leptonOnePt, leptonOneEta, leptonOnePhi;
    Float_t leptonTwoPt, leptonTwoEta, leptonTwoPhi;

    // Additional lepton data
    Float_t leptonOneIso, leptonTwoIso;
    Int_t leptonOneMother, leptonTwoMother;
    Int_t leptonOneFlavor, leptonTwoFlavor;
    Float_t leptonOneD0, leptonTwoD0;
    Float_t leptonOneDZ, leptonTwoDZ;
    Float_t leptonOneRecoWeight, leptonTwoRecoWeight;
    
    // tau data
    Int_t tauDecayMode;
    Float_t tauMVA; 
    //UInt_t tauPhotonMult, tauChHadMult;

    // photon data
    Float_t photonOnePt, photonOneEta, photonOnePhi;
    Float_t photonOneR9;
    Float_t photonOneMVA;
    Float_t photonOneERes;
    Bool_t passElectronVeto;

    Bool_t isLeptonTag;
    Bool_t isDijetTag;

    // jet data
    Float_t jetOnePt, jetOneEta, jetOnePhi, jetOneM;
    Float_t jetTwoPt, jetTwoEta, jetTwoPhi, jetTwoM;
    Float_t jetOneTag, jetTwoTag, jetThreeTag, jetFourTag;
    Float_t met, metPhi, metNC, metPhiNC, ht, htPhi, htSum;

    // generator level data
    Float_t genLeptonOnePt, genLeptonOneEta, genLeptonOnePhi;
    Float_t genLeptonTwoPt, genLeptonTwoEta, genLeptonTwoPhi;
    Float_t genPhotonPt, genPhotonEta, genPhotonPhi;
    Int_t genLeptonOneId, genLeptonTwoId;
    Bool_t genPhotonFHPFS, genPhotonIPFS;
    
    //Int_t genOneId, genTwoId, genOneMother, genTwoMother, genCategory;
    //TLorentzVector genOneP4, genTwoP4;
    //Bool_t fromHardProcessFinalState, isPromptFinalState, hasPhotonMatch;
    Bool_t vetoDY;

    // dilepton data
    Float_t dileptonPt, dileptonEta, dileptonPhi, dileptonM;
    Float_t dileptonDEta, dileptonDPhi, dileptonDR;
    
    // dilepton vertex data
    //Float_t dileptonVertexOneX, dileptonVertexOneY, dileptonVertexOneZ;
    //Float_t dileptonVertexTwoX, dileptonVertexTwoY, dileptonVertexTwoZ;
    //Float_t dileptonVertexOneXErr, dileptonVertexOneYErr, dileptonVertexOneZErr;
    //Float_t dileptonVertexTwoXErr, dileptonVertexTwoYErr, dileptonVertexTwoZErr;
    //Float_t dileptonVertexChi2One, dileptonVertexDOFOne;
    //Float_t dileptonVertexChi2Two, dileptonVertexDOFTwo;
    
    // dijet data
    Float_t dijetPt, dijetEta, dijetPhi, dijetM;
    Float_t dijetDEta, dijetDPhi, dijetDR;

    // jet, lepton data
    Float_t l1j1DEta, l1j1DPhi, l1j1DR;
    Float_t l1j2DEta, l1j2DPhi, l1j2DR;
    Float_t l2j1DEta, l2j1DPhi, l2j1DR;
    Float_t l2j2DEta, l2j2DPhi, l2j2DR;

    // jet, photon data
    Float_t j1PhotonDEta, j1PhotonDPhi, j1PhotonDR;
    Float_t j2PhotonDEta, j2PhotonDPhi, j2PhotonDR;
    Float_t jPhotonDRMax, jPhotonDRMin;

    // three body
    Float_t llgPt, llgEta, llgPhi, llgM, llgPtOverM;
    Float_t l1PhotonDEta, l1PhotonDPhi, l1PhotonDR;
    Float_t l2PhotonDEta, l2PhotonDPhi, l2PhotonDR;
    Float_t lPhotonDRMax, lPhotonDRMin;
    Float_t dileptonPhotonDEta, dileptonPhotonDPhi, dileptonPhotonDR;
    Float_t ptt;
    Float_t zgBigTheta, zgLittleTheta, zgPhi;
    Float_t genBigTheta, genLittleTheta, genPhi;

    // other
    Float_t llgJJDEta, llgJJDPhi, llgJJDR;
    Float_t zepp;

    float GetMuonIsolation(const baconhep::TMuon*);
    float GetElectronIsolation(const baconhep::TElectron*, float);
    float GetPhotonIsolation(const baconhep::TPhoton*, float);

    //ClassDef(hzgAnalyzer,0);
};

#endif  // HZGANALYZER_HH