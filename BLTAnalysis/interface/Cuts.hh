#ifndef CUTS_HH
#define CUTS_HH

#include <string>
#include <vector>

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

        float   pt,     eta,    dxy,    dz,     SIP3d;
        bool    IsGLB,  IsTRK,  IsPF,   IsStandalone,    IsLoose,    IsIsolated;
        unsigned NumberOfMatchedStations;
    }
    looseHZZMuonID, tightHZZMuonID;

    struct muIsoCuts
    {
        std::string cutName;
        float relCombIso04;
    }
    wpHZZMuonIso;



    //
    //  ELECTRONS
    //

    struct elIDCuts
    {
        std::string cutName;

        float   pt,     eta,    dxy,    dz,     SIP3d;
        bool    IsLoose,        IsMVA,          IsIsolated;
    }
    looseHZZElectronID, tightHZZElectronID;

    struct elMVACuts
    {
        std::string cutName;
        float   pt[2],  eta[3], bdt[2][3];
    }
    wpHZZElectronMVA;

    struct elIsoCuts
    {
        std::string cutName;
        float relCombIso04;
    }
    wpHZZElectronIso;

    float etaBins[8], effArea2012[7];
};

#endif  // CUTS_HH
