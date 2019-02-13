#include "BLT/BLTAnalysis/interface/Cuts.hh"

#include <iostream>
#include <stdexcept>


Cuts::Cuts() {


    //
    //  MUONS
    //

    // Cut-based ID

    looseHZZMuonID.cutName                          = "looseHZZMuonID";
    looseHZZMuonID.pt                               = 5;
    looseHZZMuonID.eta                              = 2.4;
    looseHZZMuonID.dxy                              = 0.5;
    looseHZZMuonID.dz                               = 1;
    looseHZZMuonID.SIP3d                            = 4;
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

    tightHZZMuonID.cutName                          = "tightHZZMuonID";
    tightHZZMuonID.IsLoose                          = kTRUE;
    tightHZZMuonID.IsPF                             = kTRUE;
    tightHZZMuonID.IsTrackerHighPt                  = kTRUE;
    tightHZZMuonID.IsIsolated                       = kTRUE;


    // Isolation

    wpHZZMuonIso.cutName                            = "wpHZZMuonIso";
    wpHZZMuonIso.relCombIso03                       = 0.35;


    

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
    tightHZZElectronID.IsIsolated                   = 1;

    tightHZZIsoMVAElectronID.cutName                = "tightHZZIsoMVAElectronID";
    tightHZZIsoMVAElectronID.IsLoose                = 1;
    tightHZZIsoMVAElectronID.IsMVA                  = 1;


    // MVA

    wpLooseIsoV1.cutName                            = "wpLooseIsoV1";
    wpLooseIsoV1.pt[0]                              = 5;
    wpLooseIsoV1.pt[1]                              = 10;
    wpLooseIsoV1.eta[0]                             = 0.8;
    wpLooseIsoV1.eta[1]                             = 1.479;
    wpLooseIsoV1.eta[2]                             = 2.5;
    wpLooseIsoV1.bdt[0][0]                          = -0.09564086146419018;
    wpLooseIsoV1.bdt[0][1]                          = -0.28229916981926795;
    wpLooseIsoV1.bdt[0][2]                          = -0.05466682296962322;
    wpLooseIsoV1.bdt[1][0]                          = -0.833466688584422;
    wpLooseIsoV1.bdt[1][1]                          = -0.7677000247570116;
    wpLooseIsoV1.bdt[1][2]                          = -0.6917305995653829;

    wpLooseNoIsoV1.cutName                          = "wpLooseNoIsoV1";
    wpLooseNoIsoV1.pt[0]                            = 5;
    wpLooseNoIsoV1.pt[1]                            = 10;
    wpLooseNoIsoV1.eta[0]                           = 0.8;
    wpLooseNoIsoV1.eta[1]                           = 1.479;
    wpLooseNoIsoV1.eta[2]                           = 2.5;
    wpLooseNoIsoV1.bdt[0][0]                        = -0.13285867293779202;
    wpLooseNoIsoV1.bdt[0][1]                        = -0.31765300958836074;
    wpLooseNoIsoV1.bdt[0][2]                        = -0.0799205914718861;
    wpLooseNoIsoV1.bdt[1][0]                        = -0.856871961305474;
    wpLooseNoIsoV1.bdt[1][1]                        = -0.8107642141584835;
    wpLooseNoIsoV1.bdt[1][2]                        = -0.7179265933023059;


    // Isolation

    wpHZZElectronIso.cutName                        = "wpHZZElectronIso";
    wpHZZElectronIso.relCombIso03                   = 0.35;

    etaBins[0]  = 0;                effArea2017[0]  = 0.1566;           effArea2016[0]  = 0.1703;
    etaBins[1]  = 1;                effArea2017[1]  = 0.1626;           effArea2016[1]  = 0.1715;
    etaBins[2]  = 1.479;            effArea2017[2]  = 0.1073;           effArea2016[2]  = 0.1213;
    etaBins[3]  = 2;                effArea2017[3]  = 0.0854;           effArea2016[3]  = 0.1230;
    etaBins[4]  = 2.2;              effArea2017[4]  = 0.1051;           effArea2016[4]  = 0.1635;
    etaBins[5]  = 2.3;              effArea2017[5]  = 0.1204;           effArea2016[5]  = 0.1937;
    etaBins[6]  = 2.4;              effArea2017[6]  = 0.1524;           effArea2016[6]  = 0.2393;
    etaBins[7]  = 5;
}
