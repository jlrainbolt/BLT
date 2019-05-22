#include "BLT/BLTAnalysis/interface/Cuts.hh"

#include <iostream>
#include <stdexcept>


Cuts::Cuts() {


    //
    //  MUONS
    //

    // Cut-based ID
    vetoHZZMuonID.cutName                           = "vetoHZZMuonID";
    vetoHZZMuonID.pt                                = 5;
    vetoHZZMuonID.eta                               = 2.4;
    vetoHZZMuonID.dxy                               = 0.5;
    vetoHZZMuonID.dz                                = 1;
    vetoHZZMuonID.SIP3d                             = 4;

    looseHZZMuonID.cutName                          = "looseHZZMuonID";
    looseHZZMuonID.IsVetoed                         = kFALSE;
    looseHZZMuonID.IsGLB                            = kTRUE;
    looseHZZMuonID.IsTRK                            = kTRUE;
    looseHZZMuonID.NumberOfMatchedStations          = 0;
    looseHZZMuonID.BestTrackType                    = 2;

    trackerHighPtMuonID.cutName                     = "trackerHighPtMuonID";
    trackerHighPtMuonID.pt                          = 200;
    trackerHighPtMuonID.ptFracError                 = 0.3;
    trackerHighPtMuonID.dxy                         = 0.2;
    trackerHighPtMuonID.dz                          = 0.5;
    trackerHighPtMuonID.NumberOfValidPixelHits      = 0;
    trackerHighPtMuonID.NumberOfMatchedStations     = 1;
    trackerHighPtMuonID.TrackLayersWithMeasurement  = 5;

    noIsoHZZMuonID.cutName                          = "noIsoHZZMuonID";
    noIsoHZZMuonID.IsLoose                          = kTRUE;
    noIsoHZZMuonID.IsPF                             = kTRUE;
    noIsoHZZMuonID.IsTrackerHighPt                  = kTRUE;

    tightHZZMuonID.cutName                          = "tightHZZMuonID";
    tightHZZMuonID.IsIdentified                     = kTRUE;
    tightHZZMuonID.IsIsolated                       = kTRUE;


    // Isolation
    wpHZZMuonIso.cutName                            = "wpHZZMuonIso";
    wpHZZMuonIso.relCombIso03                       = 0.35;

    wpTightMuonIso.cutName                          = "wpTightMuonIso";
    wpTightMuonIso.relCombIso03                     = 0.15;


    

    //
    //  ELECTRONS
    //

    // Cut-based ID

    looseHZZElectronID.cutName                      = "looseHZZElectronID";
    looseHZZElectronID.pt                           = 7;
    looseHZZElectronID.eta                          = 2.5;
    looseHZZElectronID.dxy                          = 0.5;
    looseHZZElectronID.dz                           = 1;
    looseHZZElectronID.SIP3d                        = 4;

    tightHZZElectronID.cutName                      = "tightHZZElectronID";
    tightHZZElectronID.IsLoose                      = 1;
    tightHZZElectronID.IsMVA                        = 1;


    // Isolation

    wpHZZElectronIso.cutName                        = "wpHZZElectronIso";
    wpHZZElectronIso.relCombIso03                   = 0.35;

    etaBins[0]  = 0;                
    etaBins[1]  = 1;                
    etaBins[2]  = 1.479;            
    etaBins[3]  = 2;                
    etaBins[4]  = 2.2;              
    etaBins[5]  = 2.3;              
    etaBins[6]  = 2.4;              
    etaBins[7]  = 5;
                                                                    // Same as 2017...
    effArea2016[0]  = 0.1703;       effArea2017[0]  = 0.1440;       effArea2018[0]  = 0.1440;
    effArea2016[1]  = 0.1715;       effArea2017[1]  = 0.1562;       effArea2018[1]  = 0.1562;
    effArea2016[2]  = 0.1213;       effArea2017[2]  = 0.1032;       effArea2018[2]  = 0.1032;
    effArea2016[3]  = 0.1230;       effArea2017[3]  = 0.0859;       effArea2018[3]  = 0.0859;
    effArea2016[4]  = 0.1635;       effArea2017[4]  = 0.1116;       effArea2018[4]  = 0.1116;
    effArea2016[5]  = 0.1937;       effArea2017[5]  = 0.1321;       effArea2018[5]  = 0.1321;
    effArea2016[6]  = 0.2393;       effArea2017[6]  = 0.1654;       effArea2018[6]  = 0.1654;
}
