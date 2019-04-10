#include "BLT/BLTAnalysis/interface/ParticleSelector.hh"

#include <iostream>

using namespace baconhep;
using namespace std;

bool test_bits(unsigned int bits, unsigned int test) {
        return (bits & test) == test;
}

ParticleSelector::ParticleSelector(const Parameters& parameters, const Cuts& cuts)
{
    this->_parameters = parameters;
    this->_cuts = cuts;

    _dataPeriod = parameters.period;
    _rng = new TRandom3(0);

    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string rcPath = cmssw_base + "/src/BLT/BLTAnalysis/data/RoccoR" + _dataPeriod + ".txt";
    _rc = new RoccoR(rcPath);
}



//
//  MUONS
//

bool ParticleSelector::PassMuonID(const baconhep::TMuon* mu, const Cuts::muIDCuts& cutLevel) const
{
    if      (cutLevel.cutName == "looseHZZMuonID")
    {
        bool isGlobal       = test_bits(mu->typeBits, baconhep::kGlobal)    == cutLevel.IsGLB;
        bool isTracker      = test_bits(mu->typeBits, baconhep::kTracker)   == cutLevel.IsTRK;
        bool isArbitrated   = mu->nMatchStn > cutLevel.NumberOfMatchedStations;

        if (
                    mu->pt          > cutLevel.pt
                &&  fabs(mu->eta)   < cutLevel.eta
                &&  fabs(mu->d0)    < cutLevel.dxy
                &&  fabs(mu->dz)    < cutLevel.dz
                &&  mu->sip3d       < cutLevel.SIP3d
                &&  mu->btt        != cutLevel.BestTrackType
                &&  (isGlobal || (isTracker && isArbitrated))
           )
            return kTRUE;
    }
    else if (cutLevel.cutName == "trackerHighPtMuonID")
    {
//      if      (_dataPeriod == "2016")
//      {
//          if (
//                      mu->pt                  > cutLevel.pt
//                  &&  fabs(mu->d0)            < cutLevel.dxy
//                  &&  fabs(mu->dz)            < cutLevel.dz
//                  &&  (mu->ptErr / mu->pt)    < cutLevel.ptFracError
//                  &&  mu->nMatchStn           > cutLevel.NumberOfMatchedStations
//                  &&  mu->nPixHits            > cutLevel.NumberOfValidPixelHits
//                  &&  mu->nTkLayers           > cutLevel.TrackLayersWithMeasurement
//             )
//              return kTRUE;
//      }
//      else if (_dataPeriod == "2017")
//      {
            bool isHighPt       = mu->pt > cutLevel.pt;
            bool isIdentified   = mu->isTrackerHighPt;

            if (isHighPt && isIdentified)
                return kTRUE;
//      }
    }
    else if (cutLevel.cutName == "tightHZZMuonID")
    {
        bool isLoose            = this->PassMuonID(mu, _cuts.looseHZZMuonID)        == cutLevel.IsLoose;
        bool isPF               = test_bits(mu->typeBits, baconhep::kPFMuon)        == cutLevel.IsPF;
        bool isTrackerHighPt    = this->PassMuonID(mu, _cuts.trackerHighPtMuonID)   == cutLevel.IsTrackerHighPt;
        bool isIsolated         = this->PassMuonIso(mu, _cuts.wpHZZMuonIso)         == cutLevel.IsIsolated;

        if (isLoose && isIsolated && (isPF || isTrackerHighPt))
            return kTRUE;
    }
    return kFALSE;
}

bool ParticleSelector::PassMuonIso(const baconhep::TMuon* mu, const Cuts::muIsoCuts& cutLevel) const
{
    if (cutLevel.cutName == "wpHZZMuonIso")
    {
        float combIso = this->GetMuonIso(mu);

        if ((combIso / mu->pt) < cutLevel.relCombIso03) 
            return kTRUE;
    }
    return kFALSE;
}

float ParticleSelector::GetMuonIso(const baconhep::TMuon* mu) const
{
    return mu->chHadIso03 + std::max(0., (double) mu->neuHadIso03 + mu->gammaIso03 - 0.5 * mu->puIso03);
}

float ParticleSelector::GetRochesterCorrection(const baconhep::TMuon* mu) const
{
    if (_isRealData)
        return _rc->kScaleDT(mu->q, mu->pt, mu->eta, mu->phi, 0, 0);
    else
        return _rc->kSmearMC(mu->q, mu->pt, mu->eta, mu->phi, mu->nTkLayers, _rng->Rndm(), 0, 0);
}



//
//  ELECTRONS
//

bool ParticleSelector::PassElectronID(const baconhep::TElectron* el, const Cuts::elIDCuts& cutLevel) const 
{
    if      (cutLevel.cutName == "looseHZZElectronID")
    {
        if (
                    el->pt          > cutLevel.pt
                &&  fabs(el->eta)   < cutLevel.eta
                &&  fabs(el->d0)    < cutLevel.dxy
                &&  fabs(el->dz)    < cutLevel.dz
                &&  el->sip3d       < cutLevel.SIP3d
           )
            return kTRUE;
    }
    else if (cutLevel.cutName == "tightHZZElectronID")
    {
        bool isLoose    = this->PassElectronID(el, _cuts.looseHZZElectronID)    == cutLevel.IsLoose;
        bool isMVA      = el->pass2017isoV2wpHZZ                                == cutLevel.IsMVA;

        if (isLoose && isMVA)
            return kTRUE;
    }
    return kFALSE;
}


// Deprecated??
bool ParticleSelector::PassElectronMVA(const baconhep::TElectron* el, const Cuts::elMVACuts& cutLevel) const
{
/*
    float bdtVal = -1;
    unsigned ptBin, etaBin;

    if      (cutLevel.cutName == "wpLooseIsoV1")
        bdtVal = el->mvaIso;
    else if (cutLevel.cutName == "wpLooseNoIsoV1")
        bdtVal = el->mva;

    if      (el->pt < cutLevel.pt[0])
        return kFALSE;
    else if (el->pt < cutLevel.pt[1])
        ptBin = 0;
    else
        ptBin = 1;

    if      (fabs(el->scEta) < cutLevel.eta[0])
        etaBin = 0;
    else if (fabs(el->scEta) < cutLevel.eta[1])
        etaBin = 1;
    else if (fabs(el->scEta) < cutLevel.eta[2])
        etaBin = 2;
    else
        return kFALSE;

    if (bdtVal > cutLevel.bdt[ptBin][etaBin])
        return kTRUE;
*/
    return kFALSE;
}


bool ParticleSelector::PassElectronIso(const baconhep::TElectron* el, const Cuts::elIsoCuts& cutLevel) const 
{
    if (cutLevel.cutName == "wpHZZElectronIso")
    {
        float combIso = this->GetElectronIso(el);

        if ((combIso / el->pt) < cutLevel.relCombIso03) 
            return kTRUE;
    }
    return kFALSE;
}

float ParticleSelector::GetElectronIso(const baconhep::TElectron* el) const 
{
    unsigned etaBin = 0;

    for (unsigned i = 0; i < 7; i++)
    {
        if (fabs(el->scEta) >= _cuts.etaBins[i] && fabs(el->scEta) < _cuts.etaBins[i+1])
        {
            etaBin = i;
            break;
        }
    }

    float effArea = 0;

    if      (_dataPeriod == "2016")
        effArea = _cuts.effArea2016[etaBin];
    else if (_dataPeriod == "2017")
        effArea = _cuts.effArea2017[etaBin];

    return el->chHadIso + std::max(0., (double) el->neuHadIso + el->gammaIso - _rhoFactor * effArea);
}

float ParticleSelector::GetElectronCorrection(const baconhep::TElectron* el) const
{
    return el->ecalTrkEnergyPostCorr / el->energy;
}
