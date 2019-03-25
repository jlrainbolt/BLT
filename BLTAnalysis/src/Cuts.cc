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
    looseHZZMuonID.IsGLB                            = 1;
    looseHZZMuonID.IsTRK                            = 1;
    looseHZZMuonID.IsStandalone                     = 0;
    looseHZZMuonID.NumberOfMatchedStations          = 0;

    tightHZZMuonID.cutName                          = "tightHZZMuonID";
    tightHZZMuonID.IsLoose                          = 1;
    tightHZZMuonID.IsPF                             = 1;
    tightHZZMuonID.IsIsolated                       = 1;


    // Isolation

    wpHZZMuonIso.cutName                            = "wpHZZMuonIso";
    wpHZZMuonIso.relCombIso04                       = 0.4;


    

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


    // MVA

    wpHZZElectronMVA.cutName                        = "wpHZZElectronMVA";
    wpHZZElectronMVA.pt[0]                          = 5;
    wpHZZElectronMVA.pt[1]                          = 10;
    wpHZZElectronMVA.eta[0]                         = 0.8;
    wpHZZElectronMVA.eta[1]                         = 1.479;
    wpHZZElectronMVA.eta[2]                         = 2.5;
    wpHZZElectronMVA.bdt[0][0]                      = 0.47;
    wpHZZElectronMVA.bdt[0][1]                      = 0.004;
    wpHZZElectronMVA.bdt[0][2]                      = 0.295;
    wpHZZElectronMVA.bdt[1][0]                      = -0.34;
    wpHZZElectronMVA.bdt[1][1]                      = -0.65;
    wpHZZElectronMVA.bdt[1][2]                      = 0.6;


    // Isolation

    wpHZZElectronIso.cutName                        = "wpHZZElectronIso";
    wpHZZElectronIso.relCombIso04                   = 0.4;

    etaBins[0]  = 0;                effArea2012[0]  = 0.208;
    etaBins[1]  = 1;                effArea2012[1]  = 0.209;
    etaBins[2]  = 1.479;            effArea2012[2]  = 0.115;
    etaBins[3]  = 2;                effArea2012[3]  = 0.143;
    etaBins[4]  = 2.2;              effArea2012[4]  = 0.183;
    etaBins[5]  = 2.3;              effArea2012[5]  = 0.194;
    etaBins[6]  = 2.4;              effArea2012[6]  = 0.261;
    etaBins[7]  = 5;
}
