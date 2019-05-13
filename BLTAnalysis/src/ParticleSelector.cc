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

    _rc = new rochcor2012();
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
        bool isStandalone   = test_bits(mu->typeBits, baconhep::kStandalone)== cutLevel.IsStandalone;

        if (
                    mu->pt          > cutLevel.pt
                &&  fabs(mu->eta)   < cutLevel.eta
                &&  fabs(mu->d0)    < cutLevel.dxy
                &&  fabs(mu->dz)    < cutLevel.dz
                &&  mu->sip3d       < cutLevel.SIP3d
                &&  (isGlobal || (isTracker && isArbitrated))
                &&  !isStandalone
           )
            return kTRUE;
    }
    else if (cutLevel.cutName == "tightHZZMuonID")
    {
        bool isLoose            = this->PassMuonID(mu, _cuts.looseHZZMuonID)        == cutLevel.IsLoose;
        bool isPF               = test_bits(mu->typeBits, baconhep::kPFMuon)        == cutLevel.IsPF;
        bool isIsolated         = this->PassMuonIso(mu, _cuts.wpHZZMuonIso)         == cutLevel.IsIsolated;

        if (isLoose && isIsolated && isPF)
            return kTRUE;
    }
    return kFALSE;
}

bool ParticleSelector::PassMuonIso(const baconhep::TMuon* mu, const Cuts::muIsoCuts& cutLevel) const
{
    if (cutLevel.cutName == "wpHZZMuonIso")
    {
        float combIso = this->GetMuonIso(mu);

        if ((combIso / mu->pt) < cutLevel.relCombIso04) 
            return kTRUE;
    }
    return kFALSE;
}

float ParticleSelector::GetMuonIso(const baconhep::TMuon* mu) const
{
    return mu->chHadIso04 + std::max(0., (double) mu->neuHadIso04 + mu->gammaIso04 - 0.5 * mu->puIso04);
}

float ParticleSelector::GetRochesterCorrection(const baconhep::TMuon* mu, std::string unc) const
{
    TLorentzVector p4;
    copy_p4(mu, MUON_MASS, p4);

    float corr = 0, qter = 1;

    if (_isRealData)
        _rc->momcor_data(p4, mu->q, 0, qter);
    else
        _rc->momcor_mc(p4, mu->q, 0, qter);

    corr = p4.Pt() / mu->pt;

    _rc->momcor_data(p4, mu->q, 0, qter);  // to get delta/qter/whatever...

    if      (unc == "")
        return corr;
    else if (unc == "up")
        return qter;
    else if (unc == "down")
        return corr - fabs(qter - corr);
    else
        return 0;
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
        bool isMVA      = this->PassElectronMVA(el, _cuts.wpHZZElectronMVA)     == cutLevel.IsMVA;
        bool isIsolated = this->PassElectronIso(el, _cuts.wpHZZElectronIso)     == cutLevel.IsIsolated;

        if (isLoose && isMVA && isIsolated)
            return kTRUE;
    }
    return kFALSE;
}


bool ParticleSelector::PassElectronMVA(const baconhep::TElectron* el, const Cuts::elMVACuts& cutLevel) const
{
    float bdtVal = -1;
    unsigned ptBin, etaBin;

    if (cutLevel.cutName == "wpHZZElectronMVA")
    {
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
    }
    return kFALSE;
}


bool ParticleSelector::PassElectronIso(const baconhep::TElectron* el, const Cuts::elIsoCuts& cutLevel) const 
{
    if (cutLevel.cutName == "wpHZZElectronIso")
    {
        float combIso = this->GetElectronIso(el);

        if ((combIso / el->pt) < cutLevel.relCombIso04) 
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

    effArea = _cuts.effArea2012[etaBin];

    return el->chHadIso04 + std::max(0., (double) el->neuHadIso04 + el->gammaIso04 - _rhoFactor * effArea);
}

float ParticleSelector::GetElectronCorrection(const baconhep::TElectron* el, std::string unc) const
{
    if      (unc == "")
        return el->ptHZZ4l / el->pt;
    else if (unc == "up")
        return (el->ptHZZ4l + el->ptErrHZZ4l) / el->pt;
    else if (unc == "down")
        return (el->ptHZZ4l - el->ptErrHZZ4l) / el->pt;
    else
        return 0;
}
