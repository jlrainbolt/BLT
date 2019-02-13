#include "PhaseSpaceAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;


bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
TLorentzVector GetP4Sum(const vector<TLorentzVector> &p4);


PhaseSpaceAnalyzer::PhaseSpaceAnalyzer() : BLTSelector()
{

}

PhaseSpaceAnalyzer::~PhaseSpaceAnalyzer()
{

}

void PhaseSpaceAnalyzer::Begin(TTree *tree)
{
    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
              std::sregex_token_iterator(), std::back_inserter(options));


    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);


    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");



    //
    //  BRANCHES
    //

    // Event
    outTree->Branch(    "runNumber",                &runNumber);
    outTree->Branch(    "evtNumber",                &evtNumber,                 "eventNumber/l");
    outTree->Branch(    "lumiSection",              &lumiSection);
    outTree->Branch(    "genWeight",                &genWeight);

    outTree->Branch(    "decayChannel",             &decayChannel);

    outTree->Branch(    "foundTauDecay",            &foundTauDecay);


    // Counters
    outTree->Branch(    "nStatus22Zs",              &nStatus22Zs);

    outTree->Branch(    "nFinalStateMuons",         &nFinalStateMuons);
    outTree->Branch(    "nFinalStateElectrons",     &nFinalStateElectrons);
    outTree->Branch(    "nFinalStateLeptons",       &nFinalStateLeptons);
    outTree->Branch(    "nFinalStateZs",            &nFinalStateZs);

    outTree->Branch(    "nHardProcMuons",           &nHardProcMuons);
    outTree->Branch(    "nHardProcElectrons",       &nHardProcElectrons);
    outTree->Branch(    "nHardProcLeptons",         &nHardProcLeptons);
//  outTree->Branch(    "nHardProcZs",              &nHardProcZs);


    // Final state leptons
    outTree->Branch(    "finalStateMuonP4",         &finalStateMuonP4,          32000,      1);
    outTree->Branch(    "finalStateMuonQ",          &finalStateMuonQ);
    outTree->Branch(    "finalStateMuonMother",     &finalStateMuonMother);
    outTree->Branch(    "finalStateMuonZIndex",     &finalStateMuonZIndex);

    outTree->Branch(    "finalStateElectronP4",     &finalStateElectronP4,      32000,      1);
    outTree->Branch(    "finalStateElectronQ",      &finalStateElectronQ);
    outTree->Branch(    "finalStateElectronMother", &finalStateElectronMother);
    outTree->Branch(    "finalStateElectronZIndex", &finalStateElectronZIndex);

    outTree->Branch(    "finalStateLeptonsP4",      &finalStateLeptonsP4);


    // Hard process leptons
    outTree->Branch(    "hardProcMuonP4",           &hardProcMuonP4,            32000,      1);
    outTree->Branch(    "hardProcMuonQ",            &hardProcMuonQ);
    outTree->Branch(    "hardProcMuonStatus",       &hardProcMuonStatus);
    outTree->Branch(    "hardProcMuonZIndex",       &hardProcMuonZIndex);

    outTree->Branch(    "hardProcElectronP4",       &hardProcElectronP4,        32000,      1);
    outTree->Branch(    "hardProcElectronQ",        &hardProcElectronQ);
    outTree->Branch(    "hardProcElectronStatus",   &hardProcElectronStatus);
    outTree->Branch(    "hardProcElectronZIndex",   &hardProcElectronZIndex);

    outTree->Branch(    "hardProcLeptonsP4",        &hardProcLeptonsP4);


    // Z bosons
    outTree->Branch(    "status22ZP4",              &status22ZP4);
    outTree->Branch(    "status22ZMother",          &status22ZMother);
    outTree->Branch(    "status22ZIndex",           &status22ZIndex);

    outTree->Branch(    "status22ZsP4",             &status22ZsP4);

    outTree->Branch(    "finalStateZP4",            &finalStateZP4);
    outTree->Branch(    "finalStateZStatus",        &finalStateZStatus);
    outTree->Branch(    "finalStateZIndex",         &finalStateZIndex);

    outTree->Branch(    "finalStateZsP4",           &finalStateZsP4);



    //
    //  HISTOGRAMS
    //

    // Total event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);

    // Acceptance counters
    outHistName = params->get_output_treename("PhaseSpaceEvents");
    hPhaseSpaceEvents = new TH1D(outHistName.c_str(), "PhaseSpaceEvents", 10, 0.5, 10.5);


    ReportPostBegin();
}


void PhaseSpaceAnalyzer::Init(TTree *tree)
{
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrentFile = tree->GetCurrentFile();

    fInfo           = 0;
    fGenEvtInfo     = 0;
    fGenParticleArr = 0;
    fLHEWeightArr   = 0;

    fChain->SetBranchAddress(   "Info",         &fInfo,             &b_Info);
    fChain->SetBranchAddress(   "GenEvtInfo",   &fGenEvtInfo,       &b_GenEvtInfo);
    fChain->SetBranchAddress(   "GenParticle",  &fGenParticleArr,   &b_GenParticleArr);
    fChain->SetBranchAddress(   "LHEWeight",    &fLHEWeightArr,     &b_LHEWeightArr);
}


Bool_t PhaseSpaceAnalyzer::Process(Long64_t entry)
{

    //
    //  CLEAR CONTAINERS
    //
    
    foundTauDecay = kFALSE;
    decayChannel = 0;

    nFinalStateMuons = 0;               nFinalStateElectrons = 0;           nFinalStateLeptons = 0; 
    nHardProcMuons = 0;                 nHardProcElectrons = 0;             nHardProcLeptons = 0; 
    nStatus22Zs = 0;                    nFinalStateZs = 0;                  nHardProcZs = 0;

    finalStateMuonP4ptr.Delete();       finalStateMuonQ.clear();            finalStateMuonMother.clear();
    finalStateElectronP4ptr.Delete();   finalStateElectronQ.clear();        finalStateElectronMother.clear();
    finalStateMuonZIndex.clear();       finalStateElectronZIndex.clear();
    finalStateLeptonsP4.Clear();

    hardProcMuonP4ptr.Delete();         hardProcMuonQ.clear();              hardProcMuonStatus.clear();
    hardProcElectronP4ptr.Delete();     hardProcElectronQ.clear();          hardProcElectronStatus.clear();
    hardProcMuonZIndex.clear();         hardProcElectronZIndex.clear();
    hardProcLeptonsP4.Clear();

    status22ZP4ptr.Delete();            status22ZMother.clear();            status22ZIndex.clear();
    status22ZsP4.Clear();
    finalStateZP4ptr.Delete();          finalStateZStatus.clear();          finalStateZIndex.clear();
    finalStateZsP4.Clear();



    //
    //  START
    //

    GetEntry(entry, 1);  // load all branches included above
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry << " Run: " << fInfo->runNum 
                  << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum 
                  << std::endl;

    const bool isData = (fInfo->runNum != 1);

    TString dataSetGroup = params->datasetgroup;
    const bool isSignal = dataSetGroup.EqualTo("zz_4l") || dataSetGroup.EqualTo("ZZTo4L");
    const bool isDrellYan = dataSetGroup.EqualTo("zjets_m-50") || dataSetGroup.EqualTo("DYJetsToLL");

    // Reject (accidental?) data events
    if (isData)
        return kTRUE;



    //
    //  EVENT INFO
    //

    genWeight = fGenEvtInfo->weight > 0 ? 1 : -1; 
    if (genWeight < 0)
        hTotalEvents->Fill(10);

    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;






    ////
    ////
    ////    PARTICLE LOOP
    ////
    ////


    hPhaseSpaceEvents->Fill(1, genWeight);

    vector<TLorentzVector>  fsMuonP4,       fsElecP4,       hpMuonP4,       hpElecP4;
    vector<int>             fsMuonQ,        fsElecQ,        hpMuonQ,        hpElecQ;
    vector<int>             fsMuonMother,   fsElecMother,   hpMuonStatus,   hpElecStatus;
    vector<unsigned>        fsMuonZIndex,   fsElecZIndex,   hpMuonZIndex,   hpElecZIndex;

    vector<TLorentzVector>  s22ZP4,         fsZP4;
    vector<int>             s22ZMother,     fsZStatus;
    vector<unsigned>        s22ZIndex;



    for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
    {
        TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);



        //
        //  LEPTONS
        //


        if (((abs(particle->pdgId) == 13) || (abs(particle->pdgId) == 11)) && particle->parent >= 0)
        {
            TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);



            // Hard process

            // ("immediate" mother is a Z)

            if (mother->pdgId == 23)
            {
                TLorentzVector p4;
                int charge = -1 * copysign(1, particle->pdgId);     // Antileptons have negative pdgId...

                if      (abs(particle->pdgId) == 13)
                {
                    p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, MUON_MASS);
                    hpMuonP4.push_back(p4);
                    hpMuonQ.push_back(charge);
                    hpMuonStatus.push_back(particle->status);
                    hpMuonZIndex.push_back(particle->parent);
                    nHardProcMuons++;
                }
                else if (abs(particle->pdgId) == 11)
                {
                    p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, ELE_MASS);
                    hpElecP4.push_back(p4);
                    hpElecQ.push_back(charge);
                    hpElecStatus.push_back(particle->status);
                    hpElecZIndex.push_back(particle->parent);
                    nHardProcElectrons++;
                }
            }



            // Final state

            // (have status 1 & can be traced to a Z)

            if (particle->status == 1)
            {
                int motherID = mother->pdgId;   // Store ID of "immediate" mother
                int motherIndex = particle->parent;

                // Trace the decay chain all the way back, allowing only e or mu as an intermediate state
                while ((abs(mother->pdgId) == 13 || abs(mother->pdgId) == 11) && mother->parent >= 0)
                {
                    motherIndex = mother->parent;
                    mother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                }

                // If the "ultimate mother" is a Z, we have a final-state lepton from a Z decay
                if (mother->pdgId == 23)
                {
                    TLorentzVector p4;
                    int charge = -1 * copysign(1, particle->pdgId);     // Antileptons have negative pdgId...

                    if      (abs(particle->pdgId) == 13)
                    {
                        p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, MUON_MASS);
                        fsMuonP4.push_back(p4);
                        fsMuonQ.push_back(charge);
                        fsMuonMother.push_back(motherID);
                        fsMuonZIndex.push_back(motherIndex);
                        nFinalStateMuons++;
                    }
                    else if (abs(particle->pdgId) == 11)
                    {
                        p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, ELE_MASS);
                        fsElecP4.push_back(p4);
                        fsElecQ.push_back(charge);
                        fsElecMother.push_back(motherID);
                        fsElecZIndex.push_back(motherIndex);
                        nFinalStateElectrons++;
                    }
                }
            }

        } // END lepton case


        else if ((abs(particle->pdgId) == 15) && particle->parent >= 0)
        {
            TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);

            // Trace the decay chain all the way back, allowing any lepton as an intermediate state
            // (maybe it could be done only allowing tau?)
            while (((abs(mother->pdgId) == 15) || (abs(mother->pdgId) == 13) || (abs(mother->pdgId) == 11)) && mother->parent >= 0)
                mother = (TGenParticle*) fGenParticleArr->At(mother->parent);

            // If the "ultimate mother" is a Z, we have a tau from a Z decay
            if (mother->pdgId == 23)
                foundTauDecay = kTRUE;
        }



        //
        //  Z BOSONS
        //

        else if ((particle->pdgId == 23) && (particle->status == 22))   // status 22 is apparently hard process
        {
            TLorentzVector p4;
            p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);

            s22ZP4.push_back(p4);
            s22ZMother.push_back(particle->parent);
            s22ZIndex.push_back(i);
            nStatus22Zs++;
        }

    } // END particle loop


    nFinalStateLeptons  = nFinalStateMuons  + nFinalStateElectrons;
    nHardProcLeptons    = nHardProcMuons    + nHardProcElectrons;






    //
    //  Z COUNTING
    //


    // Status 22
    TLorentzVector s22ZsP4;

    for (unsigned i = 0; i < nStatus22Zs; i++)
        s22ZsP4 += s22ZP4[i];



    // Final state
    vector<unsigned> fsZIndex;

    for (unsigned i = 0; i < nFinalStateMuons; i++)
    {
        // If this muon's mother Z has not already been found
        if (find(fsZIndex.begin(), fsZIndex.end(), fsMuonZIndex[i]) == fsZIndex.end())
        {
            fsZIndex.push_back(fsMuonZIndex[i]);
            nFinalStateZs++;
        }
    }
    for (unsigned i = 0; i < nFinalStateElectrons; i++)
    {
        // If this electron's mother Z has not already been found
        if (find(fsZIndex.begin(), fsZIndex.end(), fsElecZIndex[i]) == fsZIndex.end())
        {
            fsZIndex.push_back(fsElecZIndex[i]);
            nFinalStateZs++;
        }
    }


    // Get the (unique) Zs from GenParticleArray
    for (unsigned i = 0; i < nFinalStateZs; i++)
    {
        TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(fsZIndex[i]);

        TLorentzVector p4;
        p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, particle->mass);

        fsZP4.push_back(p4);
        fsZStatus.push_back(particle->status);
    }


    // Sum them all up
    TLorentzVector fsZsP4;

    for (unsigned i = 0; i < nFinalStateZs; i++)
        fsZsP4 += fsZP4[i];




    //
    //  LEPTON SUMS 
    //


    // Final state

    TLorentzVector fsMuonsP4, fsElecsP4;

    for (unsigned i = 0; i < nFinalStateMuons; i++)
        fsMuonsP4 += fsMuonP4[i];
    for (unsigned i = 0; i < nFinalStateElectrons; i++)
        fsElecsP4 += fsElecP4[i];

    TLorentzVector fsLepsP4 = fsMuonsP4 + fsElecsP4;



    // Hard process

    TLorentzVector hpMuonsP4, hpElecsP4;

    for (unsigned i = 0; i < nHardProcMuons; i++)
        hpMuonsP4 += hpMuonP4[i];
    for (unsigned i = 0; i < nHardProcElectrons; i++)
        hpElecsP4 += hpElecP4[i];

    TLorentzVector hpLepsP4 = hpMuonsP4 + hpElecsP4;






    ////
    ////
    ////    SELECTION
    ////
    ////


    //
    //  FINAL STATE
    //

    if      (params->selection == "final")
    {
        // Require four final-state leptons
        // (because they are what will be detected)

        if (isSignal && (nFinalStateLeptons != 4))
            return kTRUE;
        if (isDrellYan && (nFinalStateLeptons != 2))
            return kTRUE;
        hTotalEvents->Fill(2);



        // 4-lepton mass

        if (fsLepsP4.M() < M_MIN || fsLepsP4.M() > M_MAX)
            return kTRUE;
        hTotalEvents->Fill(3);



        // Dilepton mass

        for (unsigned j = 1; j < nFinalStateMuons; j++)
        {
            for (unsigned i = 0; i < j; i++)
            {
                if (fsMuonQ[i] != fsMuonQ[j])
                {
                    TLorentzVector dimuonP4 = fsMuonP4[i] + fsMuonP4[j];

                    if (dimuonP4.M() < MLL_MIN)
                        return kTRUE;
                }
            }
        }
        for (unsigned j = 1; j < nFinalStateElectrons; j++)
        {
            for (unsigned i = 0; i < j; i++)
            {
                if (fsElecQ[i] != fsElecQ[j])
                {
                    TLorentzVector dielecP4 = fsElecP4[i] + fsElecP4[j];

                    if (dielecP4.M() < MLL_MIN)
                        return kTRUE;
                }
            }
        }
        hTotalEvents->Fill(4);



        // Categorize

        hPhaseSpaceEvents->Fill(1, genWeight);
        unsigned C = 0;                                                 // Index

        if      (nFinalStateMuons == 2 && nFinalStateElectrons == 0)    // mumu = 3
            C = 3;

        else if (nFinalStateMuons == 0 && nFinalStateElectrons == 2)    // ee   = 4
            C = 4;

        else if (nFinalStateMuons == 4 && nFinalStateElectrons == 0)    // 4m   = 6
            C = 6;

        else if (nFinalStateMuons == 2 && nFinalStateElectrons == 2     // 2m2e = 7
                && fsMuonsP4.M() > fsElecsP4.M())
            C = 7;

        else if (nFinalStateMuons == 2 && nFinalStateElectrons == 2     // 2e2m = 8
                && fsMuonsP4.M() < fsElecsP4.M())
            C = 8;

        else if (nFinalStateMuons == 0 && nFinalStateElectrons == 4)    // 4e   = 9
            C = 9;

        unsigned D = (C < 6) ? 2 : 5;
        hPhaseSpaceEvents->Fill(C, genWeight);
        hPhaseSpaceEvents->Fill(D, genWeight);
        decayChannel = C;


    } // END "final" case



    //
    //  HARD PROCESS
    //

    else if (params->selection == "hard")
    {
        // Require four hard-process leptons
        // (because they are what is "real")


        if (isSignal && (nHardProcLeptons != 4))
            return kTRUE;
        if (isDrellYan && (nHardProcLeptons != 2))
            return kTRUE;
        hTotalEvents->Fill(2);



        // 4-lepton mass

        if (hpLepsP4.M() < M_MIN || hpLepsP4.M() > M_MAX)
            return kTRUE;
        hTotalEvents->Fill(3);



        // Dilepton mass

        for (unsigned j = 1; j < nHardProcMuons; j++)
        {
            for (unsigned i = 0; i < j; i++)
            {
                if (hpMuonQ[i] != hpMuonQ[j])
                {
                    TLorentzVector dimuonP4 = hpMuonP4[i] + hpMuonP4[j];

                    if (dimuonP4.M() < MLL_MIN)
                        return kTRUE;
                }
            }
        }
        for (unsigned j = 1; j < nHardProcElectrons; j++)
        {
            for (unsigned i = 0; i < j; i++)
            {
                if (hpElecQ[i] != hpElecQ[j])
                {
                    TLorentzVector dielecP4 = hpElecP4[i] + hpElecP4[j];

                    if (dielecP4.M() < MLL_MIN)
                        return kTRUE;
                }
            }
        }
        hTotalEvents->Fill(4);



        // Categorize

        hPhaseSpaceEvents->Fill(1, genWeight);
        unsigned C = 0;                                             // Index

        if      (nHardProcMuons == 2 && nHardProcElectrons == 0)    // mumu = 3
            C = 3;

        else if (nHardProcMuons == 0 && nHardProcElectrons == 2)    // ee   = 4
            C = 4;

        else if (nHardProcMuons == 4 && nHardProcElectrons == 0)    // 4m   = 6
            C = 6;

        else if (nHardProcMuons == 2 && nHardProcElectrons == 2     // 2m2e = 7
                && hpMuonsP4.M() > hpElecsP4.M())
            C = 7;

        else if (nHardProcMuons == 2 && nHardProcElectrons == 2     // 2e2m = 8
                && hpMuonsP4.M() < hpElecsP4.M())
            C = 8;

        else if (nHardProcMuons == 0 && nHardProcElectrons == 4)    // 4e   = 9
            C = 9;

        unsigned D = (C < 6) ? 2 : 5;
        hPhaseSpaceEvents->Fill(C, genWeight);
        hPhaseSpaceEvents->Fill(D, genWeight);
        decayChannel = C;


    } // END "hard" case

    else if (params->selection != "all")
        return kTRUE;

    hTotalEvents->Fill(9);




    //
    //  DEBUG
    //


    if (nHardProcLeptons != nFinalStateLeptons)
//  if (kTRUE)
    {
        cout << nFinalStateMuons << " fs muons\t" << nFinalStateElectrons << " fs elecs" << endl;
        cout << nHardProcMuons << " hp muons\t" << nHardProcElectrons << " hp elecs" << endl;
        cout << endl;
        cout << "Idx" << "\t" << "ID" << "\t" << "Stat" << "\t" << "Mom" << "\t" << "Pt" << endl;
        for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if (abs(particle->pdgId) == 11 || abs(particle->pdgId) == 13 || abs(particle->pdgId) == 15 || abs(particle->pdgId) == 23)
                cout << i << "\t" << particle->pdgId << "\t" << particle->status << "\t" << particle->parent << "\t" << particle->pt << endl;
        }
        cout << endl;

        cout << "hp muons:\t";
        for (unsigned i = 0; i < nHardProcMuons; i++)
            cout << hpMuonZIndex[i] << ", ";
        cout << endl;
        cout << "hp elecs:\t";
        for (unsigned i = 0; i < nHardProcElectrons; i++)
            cout << hpElecZIndex[i] << ", ";
        cout << endl << endl;
    }



    //
    //  FILL TREE
    //

    // Final state

    // Muons
    for (unsigned i = 0; i < nFinalStateMuons; i++)
    {
        new(finalStateMuonP4ptr[i]) TLorentzVector(fsMuonP4[i]);
        finalStateMuonQ.push_back(fsMuonQ[i]);
        finalStateMuonMother.push_back(fsMuonMother[i]);
        finalStateMuonZIndex.push_back(fsMuonZIndex[i]);
    }

    // Electrons
    for (unsigned i = 0; i < nFinalStateElectrons; i++)
    {
        new(finalStateElectronP4ptr[i]) TLorentzVector(fsElecP4[i]);
        finalStateElectronQ.push_back(fsElecQ[i]);
        finalStateElectronMother.push_back(fsElecMother[i]);
        finalStateElectronZIndex.push_back(fsElecZIndex[i]);
    }

    finalStateLeptonsP4 = fsLepsP4;



    // Hard process

    // Muons
    for (unsigned i = 0; i < nHardProcMuons; i++)
    {
        new(hardProcMuonP4ptr[i]) TLorentzVector(hpMuonP4[i]);
        hardProcMuonQ.push_back(hpMuonQ[i]);
        hardProcMuonStatus.push_back(hpMuonStatus[i]);
        hardProcMuonZIndex.push_back(hpMuonZIndex[i]);
    }

    // Electrons
    for (unsigned i = 0; i < nHardProcElectrons; i++)
    {
        new(hardProcElectronP4ptr[i]) TLorentzVector(hpElecP4[i]);
        hardProcElectronQ.push_back(hpElecQ[i]);
        hardProcElectronStatus.push_back(hpElecStatus[i]);
        hardProcElectronZIndex.push_back(hpElecZIndex[i]);
    }

    hardProcLeptonsP4 = hpLepsP4;



    // Z bosons

    // Status 22
    for (unsigned i = 0; i < nStatus22Zs; i++)
    {
        new(status22ZP4ptr[i]) TLorentzVector(s22ZP4[i]);
        status22ZMother.push_back(s22ZMother[i]);
        status22ZIndex.push_back(s22ZIndex[i]);
    }
    status22ZsP4 = s22ZsP4;

    // "Final state" (traced from FS leptons)
    for (unsigned i = 0; i < nFinalStateZs; i++)
    {
        new(finalStateZP4ptr[i]) TLorentzVector(fsZP4[i]);
        finalStateZStatus.push_back(fsZStatus[i]);
        finalStateZIndex.push_back(fsZIndex[i]);
    }
    finalStateZsP4 = fsZsP4;




    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void PhaseSpaceAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void PhaseSpaceAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void PhaseSpaceAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job =========================================" << std::endl;
    std::cout << "  output   : " << params->get_output_filename("output") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
    std::cout << "  ============================================================" << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<PhaseSpaceAnalyzer> selector(new PhaseSpaceAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}




///////////////////
//    HELPERS    //
///////////////////


TLorentzVector GetP4Sum(const vector<TLorentzVector> &p4)
{
    TLorentzVector p4sum;

    for (unsigned i = 0; i < p4.size(); i++)
        p4sum += p4[i];

    return p4sum;
}
