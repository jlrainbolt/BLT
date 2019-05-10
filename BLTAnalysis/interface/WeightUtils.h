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
        void    SetSampleName(string);
        void    SetSelection(string);

        float   GetPUWeight(float) const;
        float   GetSampleWeight() const;
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

        const float XSEC_4x = 0.07691,  XSEC_2x2y = 0.1767,     XSEC_4L = 3 * (XSEC_4x + XSEC_2x2y);

        const float NGEN_4e = 1499093,      NGEN_4mu = 1499064,     NGEN_4tau = 824466;
        const float NGEN_2e2mu = 1497445,   NGEN_2e2tau = 823911,   NGEN_2mu2tau = 823922;
        const float NGEN_4L = NGEN_4e + NGEN_4mu + NGEN_4tau + NGEN_2e2mu + NGEN_2e2tau + NGEN_2mu2tau;
};

#endif

#if !defined(__CINT__)
ClassImp(WeightUtils);
#endif
