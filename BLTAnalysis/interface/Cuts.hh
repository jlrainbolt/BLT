#ifndef CUTS_HH
#define CUTS_HH

#include <string>
#include <vector>
#include <TRandom3.h>
#include "BLT/BLTAnalysis/interface/RoccoR.h"

class Cuts {
public:
    Cuts();
    ~Cuts() {}



    //
    //  MUONS
    //

    struct muIDCuts
    {
        std::string cutName;

        float   pt,     eta,    dxy,    dz,     SIP3d,      ptFracError;
        bool    IsGLB,  IsTRK,  IsPF,   IsTrackerHighPt;
        bool    IsVetoed,       IsLoose,    IsIdentified,   IsIsolated;
        float   NumberOfMatchedStations,    NumberOfValidPixelHits,     TrackLayersWithMeasurement;
        int     BestTrackType;
    }
    vetoHZZMuonID, looseHZZMuonID, tightHZZMuonID, noIsoHZZMuonID, trackerHighPtMuonID;

    struct muIsoCuts
    {
        std::string cutName;
        float relCombIso03;
    }
    wpHZZMuonIso, wpTightMuonIso;



    //
    //  ELECTRONS
    //

    struct elIDCuts
    {
        std::string cutName;

        float   pt,     eta,    dxy,    dz,     SIP3d;
        bool    IsLoose,        IsMVA;
    }
    looseHZZElectronID, tightHZZElectronID;

    struct elIsoCuts
    {
        std::string cutName;
        float relCombIso03;
    }
    wpHZZElectronIso;

    float etaBins[8], effArea2016[7], effArea2017[7], effArea2018[7];
};

#endif  // CUTS_HH
