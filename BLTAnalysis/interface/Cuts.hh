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
        bool    IsGLB,  IsTRK,  IsPF,   IsTrackerHighPt,    IsLoose,    IsIsolated;
        float   NumberOfMatchedStations,    NumberOfValidPixelHits,     TrackLayersWithMeasurement;
        int     BestTrackType;
    }
    looseHZZMuonID, tightHZZMuonID, trackerHighPtMuonID;

    struct muIsoCuts
    {
        std::string cutName;
        float relCombIso03;
    }
    wpHZZMuonIso;



    //
    //  ELECTRONS
    //

    struct elIDCuts
    {
        std::string cutName;

        float   pt,     eta,    dxy,    dz,     SIP3d;
        bool    IsLoose,    IsMVA,  IsIsolated;
    }
    looseHZZElectronID, tightHZZElectronID, tightHZZIsoMVAElectronID;

    struct elMVACuts
    {
        std::string cutName;
        float   pt[2],  eta[3], bdt[2][3];
    }
    wpLooseIsoV1, wpLooseNoIsoV1;

    struct elIsoCuts
    {
        std::string cutName;
        float relCombIso03;
    }
    wpHZZElectronIso;

    float etaBins[8], effArea2016[7], effArea2017[7];
};

#endif  // CUTS_HH
