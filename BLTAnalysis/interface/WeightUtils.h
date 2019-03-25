/*
   Utilities for retrieving weights for PU,etc.
 */

#ifndef _WeightUtils_H
#define _WeightUtils_H

#include "BLT/BLTAnalysis/interface/BLTHelper.hh"

// c++ libraries
#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

// ROOT libraries
#include "TROOT.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

// BaconAna class definitions
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"

using namespace std;



class WeightUtils: public TObject {
    public:
        WeightUtils() {};
        virtual ~WeightUtils() {};
        WeightUtils(string dataPeriod, string selection, bool isRealData);

        void    Initialize();
        void    SetDataBit(bool);
        void    SetDataPeriod(string);
        void    SetSelection(string);

        float   GetPUWeight(float) const;
        float   GetHZZMuonIDSF(const baconhep::TMuon*) const;
        float   GetHZZElectronIDSF(const baconhep::TElectron*) const;
        float   GetElectronRecoSF(const baconhep::TElectron*) const;

        std::pair<float, float> GetDoubleMuonTriggerEff(const baconhep::TMuon*, const int) const;
        std::pair<float, float> GetDoubleElectronTriggerEff(const baconhep::TElectron*, const int) const;

        ClassDef(WeightUtils, 0);

    private:
        // Input parameters
        string _dataPeriod;
        string _sampleName;
        string _selection;
        bool   _isRealData;

        // Pileup
        TH1 *_puReweight;

        // Muon triggers, ID
//      TH2D *h2_MuTriggerSFs[2]; // Good for Mu17_Mu8 or Mu17_TkMu8
        TH2 *_hzz_muIdSF;

        // Electron triggers, ID
//      TH2F *_eff_doubleEle_leg1_DATA, *_eff_doubleEle_leg1_MC;
//      TH2F *_eff_doubleEle_leg2_DATA, *_eff_doubleEle_leg2_MC;
        TH2 *_hzz_eleIdSF;  // Includes ID and reco
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif
