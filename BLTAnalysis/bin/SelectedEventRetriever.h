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


#ifndef SELECTEDEVENTRETRIEVER_HH
#define SELECTEDEVENTRETRIEVER_HH

// Analysis tools
#include "BLT/BLTAnalysis/interface/BLTSelector.hh"
#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"

// ROOT headers
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

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


class SelectedEventRetriever: public BLTSelector {
public:
    SelectedEventRetriever();
    ~SelectedEventRetriever();

    void    Begin(TTree *tree);
    Bool_t  Process(Long64_t entry);
    void    Terminate();
    void    ReportPostBegin();
    void    ReportPostTerminate();

    TFile *outFile;
    TTree *outTree;

    // Params
    std::unique_ptr<Parameters>         params;

    // Info about selected events
    UInt_t nEvents;
    std::vector<UInt_t> runNumber, lumiSection, evtNumber;



    //ClassDef(SelectedEventRetriever,0);
};


#endif  // SELECTEDEVENTRETRIEVER_HH
