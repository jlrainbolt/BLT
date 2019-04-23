#ifndef PARTICLESELECTOR_HH
#define PARTICLESELECTOR_HH


// Bacon header files
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"

#include "BLT/BLTAnalysis/interface/BLTHelper.hh"
#include "BLT/BLTAnalysis/interface/Parameters.hh"
#include "BLT/BLTAnalysis/interface/Cuts.hh"
#include "BLT/BLTAnalysis/interface/RoccoR.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"

#include <string>
#include <vector>
#include <memory>
#include <cassert>

using namespace std;

class ParticleSelector {
public:
    ParticleSelector(const Parameters& parameters, const Cuts& cuts);
    ~ParticleSelector() {}

    // Setters
    void SetRealData(bool isRealData)       { _isRealData = isRealData; }
    void SetPV(const TVector3& pv)          { _pv = pv; }
    void SetNPV(int npv)                    { _npv = npv; }
    void SetRho(float rhoFactor)            { _rhoFactor = rhoFactor; }

    // Muons
    bool PassMuonID(const baconhep::TMuon* mu, const Cuts::muIDCuts& cutLevel) const;
    bool PassMuonIso(const baconhep::TMuon* mu, const Cuts::muIsoCuts& cutLevel) const;
    float GetMuonIso(const baconhep::TMuon* mu) const;
    float GetRochesterCorrection(const baconhep::TMuon* mu) const;

    // Electrons
    bool PassElectronID(const baconhep::TElectron* el, const Cuts::elIDCuts& cutLevel) const;
    bool PassElectronIso(const baconhep::TElectron* el, const Cuts::elIsoCuts& cutLevel) const;
    float GetElectronIso(const baconhep::TElectron* el) const;
    float GetElectronCorrection(const baconhep::TElectron* el) const;

private:
    Parameters  _parameters;
    Cuts        _cuts;
    bool        _isRealData;
    TVector3    _pv;
    int         _npv;
    float       _rhoFactor;
    std::string _dataPeriod;
    TRandom3*   _rng;
    RoccoR*     _rc;
};

#endif  // PARTICLESELECTOR_HH
