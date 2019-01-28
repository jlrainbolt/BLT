#include "MultileptonAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 

TLorentzVector GetP4Sum(vector<TLorentzVector> p4)
{
    TLorentzVector sum(0, 0, 0, 0);
    for (auto p4_ = p4.begin(); p4_ != p4.end(); ++p4_)
        sum += *p4_;
    return sum;
}


MultileptonAnalyzer::MultileptonAnalyzer() : BLTSelector()
{

}

MultileptonAnalyzer::~MultileptonAnalyzer()
{

}

void MultileptonAnalyzer::Begin(TTree *tree)
{
    rng = new TRandom3();

    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
              std::sregex_token_iterator(), std::back_inserter(options));


    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);

    TString dataSetGroup    = params->datasetgroup;
    const bool isSignal     = dataSetGroup.EqualTo("zz_4l") || dataSetGroup.EqualTo("ZZTo4L");


    // Set the cuts
//  cuts.reset(new Cuts());
//  particleSelector.reset(new ParticleSelector(*params, *cuts));


    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    // https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2016
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    if (params->selection == "double")
    {
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
        muonTriggerNames.push_back("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        muonTriggerNames.push_back("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

        electronTriggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    }
    else
        cout << "WARNING: No triggers selected!" << endl;


    // Weight utility class
    weights.reset(new WeightUtils(params->period, params->selection, false)); // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    lumiMask.AddJSONFile(jsonFileName);

    // muon momentum corrections
    muonCorr = new RoccoR(cmssw_base + "/src/BLT/BLTAnalysis/data/RoccoR2016.txt");


    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");



    ////  BRANCHES

    // Event
    outTree->Branch(    "runNumber",                &runNumber);
    outTree->Branch(    "evtNumber",                &evtNumber,             "eventNumber/l");
    outTree->Branch(    "lumiSection",              &lumiSection);
    outTree->Branch(    "nPV",                      &nPV);

    outTree->Branch(    "genWeight",                &genWeight);
    outTree->Branch(    "PUWeight",                 &PUWeight);
    outTree->Branch(    "ECALWeight",               &ECALWeight);
    outTree->Branch(    "nPU",                      &nPU);

    outTree->Branch(    "evtMuonTriggered",         &evtMuonTriggered);
    outTree->Branch(    "evtElectronTriggered",     &evtElectronTriggered);


    // Counters
    outTree->Branch(    "nMuons",                   &nMuons);
    outTree->Branch(    "nElectrons",               &nElectrons);
    outTree->Branch(    "nLeptons",                 &nLeptons);

    outTree->Branch(    "nV1Electrons",             &nV1Electrons);
    outTree->Branch(    "nV2Electrons",             &nV2Electrons);
    outTree->Branch(    "nLooseMuons",              &nLooseMuons);

    outTree->Branch(    "nHZZMuons",                &nHZZMuons);
    outTree->Branch(    "nHZZElectrons",            &nHZZElectrons);
    outTree->Branch(    "nHZZLeptons",              &nHZZLeptons);


    // Muons
    outTree->Branch(    "muonP4",                   &muonsP4,               32000, 1);
    outTree->Branch(    "muonQ",                    &muonsQ);
    outTree->Branch(    "muonFiredLeg1",            &muonFiredLeg1);
    outTree->Branch(    "muonFiredLeg2",            &muonFiredLeg2);

    outTree->Branch(    "muonIsGhost",              &muonIsGhost);
    outTree->Branch(    "muonIsLoose",              &muonIsLoose);
    outTree->Branch(    "muonIsHZZ",                &muonIsHZZ);

    outTree->Branch(    "muonEnergySF",             &muonEnergySF);
    outTree->Branch(    "muonHZZIDSF",              &muonHZZIDSF);
    outTree->Branch(    "muonTrigEffLeg1Data",      &muonTrigEffLeg1Data);
    outTree->Branch(    "muonTrigEffLeg1MC",        &muonTrigEffLeg1MC);
    outTree->Branch(    "muonTrigEffLeg2Data",      &muonTrigEffLeg2Data);
    outTree->Branch(    "muonTrigEffLeg2MC",        &muonTrigEffLeg2MC);

    outTree->Branch(    "muonCombIso",              &muonCombIso);
    outTree->Branch(    "muonTrkIso",               &muonsTrkIso);
    outTree->Branch(    "muonD0",                   &muonD0);
    outTree->Branch(    "muonDz",                   &muonDz);
    outTree->Branch(    "muonSIP3d",                &muonSIP3d);
    outTree->Branch(    "muonPtErr",                &muonPtErr);
    outTree->Branch(    "muonNMatchStn",            &muonNMatchStn);
    outTree->Branch(    "muonNPixHits",             &muonNPixHits);
    outTree->Branch(    "muonNTkLayers",            &muonNTkLayers);
    outTree->Branch(    "muonIsPF",                 &muonIsPF);
    outTree->Branch(    "muonIsGlobal",             &muonIsGlobal);
    outTree->Branch(    "muonIsTracker",            &muonIsTracker);
    outTree->Branch(    "muonBestTrackType",        &muonBestTrackType);


    // Electrons
    outTree->Branch(    "electronP4",               &electronsP4,           32000, 1);
    outTree->Branch(    "electronQ",                &electronsQ);
    outTree->Branch(    "electronFiredLeg1",        &electronFiredLeg1);
    outTree->Branch(    "electronFiredLeg2",        &electronFiredLeg2);

    outTree->Branch(    "electronIsGhost",          &electronIsGhost);
    outTree->Branch(    "electronPassV1MVA",        &electronPassV1MVA);
    outTree->Branch(    "electronPassV2MVA",        &electronPassV2MVA);
    outTree->Branch(    "electronIsHZZ",            &electronIsHZZ);

    outTree->Branch(    "electronEnergySF",         &electronEnergySF);
    outTree->Branch(    "electronHZZIDSF",          &electronHZZIDSF);
    outTree->Branch(    "electronTrigEffLeg1Data",  &electronTrigEffLeg1Data);
    outTree->Branch(    "electronTrigEffLeg1MC",    &electronTrigEffLeg1MC);
    outTree->Branch(    "electronTrigEffLeg2Data",  &electronTrigEffLeg2Data);
    outTree->Branch(    "electronTrigEffLeg2MC",    &electronTrigEffLeg2MC);

    outTree->Branch(    "electronCombIso",          &electronCombIso);
    outTree->Branch(    "electronTrkIso",           &electronsTrkIso);
    outTree->Branch(    "electronD0",               &electronD0);
    outTree->Branch(    "electronDz",               &electronDz);
    outTree->Branch(    "electronSIP3d",            &electronSIP3d);
    outTree->Branch(    "electronScEta",            &electronScEta);
    outTree->Branch(    "electronNMissHits",        &electronNMissHits);
    outTree->Branch(    "electronIsGap",            &electronIsGap);


    // Gen particles
    if (isSignal)
    {
        // Counters
        outTree->Branch(    "nFinalStateMuons",         &nFinalStateMuons);
        outTree->Branch(    "nFinalStateElectrons",     &nFinalStateElectrons);
        outTree->Branch(    "nFinalStateLeptons",       &nFinalStateLeptons);

        outTree->Branch(    "nHardProcMuons",           &nHardProcMuons);
        outTree->Branch(    "nHardProcElectrons",       &nHardProcElectrons);
        outTree->Branch(    "nHardProcLeptons",         &nHardProcLeptons);


        // Final state leptons
        outTree->Branch(    "finalStateMuonP4",         &finalStateMuonP4,          32000,      1);
        outTree->Branch(    "finalStateMuonQ",          &finalStateMuonQ);
        outTree->Branch(    "finalStateMuonZIndex",     &finalStateMuonZIndex);

        outTree->Branch(    "finalStateElectronP4",     &finalStateElectronP4,      32000,      1);
        outTree->Branch(    "finalStateElectronQ",      &finalStateElectronQ);
        outTree->Branch(    "finalStateElectronZIndex", &finalStateElectronZIndex);

        outTree->Branch(    "finalStateLeptonsP4",      &finalStateLeptonsP4);


        // Hard process leptons
        outTree->Branch(    "hardProcMuonP4",           &hardProcMuonP4,            32000,      1);
        outTree->Branch(    "hardProcMuonQ",            &hardProcMuonQ);
        outTree->Branch(    "hardProcMuonZIndex",       &hardProcMuonZIndex);

        outTree->Branch(    "hardProcElectronP4",       &hardProcElectronP4,        32000,      1);
        outTree->Branch(    "hardProcElectronQ",        &hardProcElectronQ);
        outTree->Branch(    "hardProcElectronZIndex",   &hardProcElectronZIndex);

        outTree->Branch(    "hardProcLeptonsP4",        &hardProcLeptonsP4);
    }


    // Histograms
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);


    ReportPostBegin();
}


Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{


    ////  CLEAR CONTAINERS

    nMuons = 0;                         nElectrons = 0;                     nLeptons = 0; 
    nLooseMuons = 0;                    nV1Electrons = 0;                   nV2Electrons = 0;
    nHZZMuons = 0;                      nHZZElectrons = 0;                  nHZZLeptons = 0; 
    nGenMuons = 0;                      nGenElectrons = 0;                  nGenLeptons = 0; 

    muonsP4ptr.Delete();                muonsQ.clear();                     muonFiredLeg1.clear();              muonFiredLeg2.clear();
    muonIsGhost.clear();                muonIsLoose.clear();                muonIsHZZ.clear();
    muonEnergySF.clear();               muonHZZIDSF.clear();
    muonTrigEffLeg1Data.clear();        muonTrigEffLeg1MC.clear();          muonTrigEffLeg2Data.clear();        muonTrigEffLeg2MC.clear(); 
    muonCombIso.clear();                muonsTrkIso.clear();                muonD0.clear();                     muonDz.clear();                     
    muonSIP3d.clear();                  muonPtErr.clear();                  muonBestTrackType.clear();
    muonNMatchStn.clear();              muonNPixHits.clear();               muonNTkLayers.clear();
    muonIsPF.clear();                   muonIsGlobal.clear();               muonIsTracker.clear();
                                                                                                        
    electronsP4ptr.Delete();            electronsQ.clear();                 electronFiredLeg1.clear();          electronFiredLeg2.clear();
    electronIsGhost.clear();            electronPassV1MVA.clear();          electronPassV2MVA.clear();          electronIsHZZ.clear();
    electronEnergySF.clear();           electronHZZIDSF.clear();
    electronTrigEffLeg1Data.clear();    electronTrigEffLeg1MC.clear();      electronTrigEffLeg2Data.clear();    electronTrigEffLeg2MC.clear();
    electronCombIso.clear();            electronsTrkIso.clear();            electronD0.clear();                 electronDz.clear();
    electronSIP3d.clear();              electronScEta.clear();
    electronNMissHits.clear();          electronIsGap.clear();

    nFinalStateMuons = 0;               nFinalStateElectrons = 0;           nFinalStateLeptons = 0; 
    nHardProcMuons = 0;                 nHardProcElectrons = 0;             nHardProcLeptons = 0; 

    finalStateMuonP4ptr.Delete();       finalStateMuonQ.clear();            finalStateMuonZIndex.clear();
    finalStateElectronP4ptr.Delete();   finalStateElectronQ.clear();        finalStateElectronZIndex.clear();
    finalStateLeptonsP4.Clear();

    hardProcMuonP4ptr.Delete();         hardProcMuonQ.clear();              hardProcMuonZIndex.clear();
    hardProcElectronP4ptr.Delete();     hardProcElectronQ.clear();          hardProcElectronZIndex.clear();
    hardProcLeptonsP4.Clear();





    /////////////////////
    //      START      //
    /////////////////////
    
    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry << " Run: " << fInfo->runNum 
                  << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum 
                  << std::endl;

    const bool isData = (fInfo->runNum != 1);
//  particleSelector->SetRealData(isData);
    weights->SetDataBit(isData);
   
    TString dataSetGroup    = params->datasetgroup;
    const bool  isSignal    = dataSetGroup.Contains("zz_4l") || dataSetGroup.Contains("ZZTo4L");


    ////  EVENT WEIGHTS

    genWeight   = 1;
    PUWeight    = 1;
    ECALWeight  = 1;
    nPU         = 0;
    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;
    nPV         = fPVArr->GetEntries();


    if (!isData)
    {
        // Save gen weight for amc@nlo Drell-Yan sample
        genWeight = fGenEvtInfo->weight > 0 ? 1 : -1; 
        if (genWeight < 0)
            hTotalEvents->Fill(10);


        // Pileup reweighting
        nPU = fInfo->nPUmean;
        PUWeight = weights->GetPUWeight(nPU);


        // L1 ECAL prefiring weight
        ECALWeight = fInfo->ecalWeight;
    }




    /////////////////////
    //  GEN PARTICLES  //
    /////////////////////

    // Adapted from PhaseSpaceAnalyzer


    vector<TLorentzVector>  fsMuonP4,       fsElecP4,       hpMuonP4,       hpElecP4;
    vector<int>             fsMuonQ,        fsElecQ,        hpMuonQ,        hpElecQ;
    vector<unsigned>        fsMuonZIndex,   fsElecZIndex,   hpMuonZIndex,   hpElecZIndex;

    if (isSignal)
    {
        ////  PARTICLE LOOP

        for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if  ((abs(particle->pdgId) == 13 || abs(particle->pdgId) == 11)  &&  particle->parent >= 0)
            {
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);


                // Hard process

                if (mother->pdgId == 23)
                {
                    TLorentzVector p4;
                    int charge = -1 * copysign(1, particle->pdgId);

                    if      (abs(particle->pdgId) == 13)
                    {
                        p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, MUON_MASS);
                        hpMuonP4.push_back(p4);
                        hpMuonQ.push_back(charge);
                        hpMuonZIndex.push_back(particle->parent);
                        nHardProcMuons++;
                    }
                    else if (abs(particle->pdgId) == 11)
                    {
                        p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, ELE_MASS);
                        hpElecP4.push_back(p4);
                        hpElecQ.push_back(charge);
                        hpElecZIndex.push_back(particle->parent);
                        nHardProcElectrons++;
                    }
                } // END hard process case


                
                // Final state

                if (particle->status == 1)
                {
                    int motherIndex = particle->parent;

                    while ((abs(mother->pdgId) == 13 || abs(mother->pdgId) == 11) && mother->parent >= 0)
                    {
                        motherIndex = mother->parent;
                        mother = (TGenParticle*) fGenParticleArr->At(mother->parent);
                    }

                    if (mother->pdgId == 23)
                    {
                        TLorentzVector p4;
                        int charge = -1 * copysign(1, particle->pdgId);

                        if      (abs(particle->pdgId) == 13)
                        {
                            p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, MUON_MASS);
                            fsMuonP4.push_back(p4);
                            fsMuonQ.push_back(charge);
                            fsMuonZIndex.push_back(motherIndex);
                            nFinalStateMuons++;
                        }
                        else if (abs(particle->pdgId) == 11)
                        {
                            p4.SetPtEtaPhiM(particle->pt, particle->eta, particle->phi, ELE_MASS);
                            fsElecP4.push_back(p4);
                            fsElecQ.push_back(charge);
                            fsElecZIndex.push_back(motherIndex);
                            nFinalStateElectrons++;
                        }
                    }
                } // END final state case
            }

        } // END gen particle loop

        nFinalStateLeptons  = nFinalStateMuons  + nFinalStateElectrons;
        nHardProcLeptons    = nHardProcMuons    + nHardProcElectrons;



        //
        //  TOTAL MOMENTA
        //

        // Hard process
        TLorentzVector hpMuonsP4, hpElecsP4;

        for (unsigned i = 0; i < nHardProcMuons; i++)
            hpMuonsP4 += hpMuonP4[i];
        for (unsigned i = 0; i < nHardProcElectrons; i++)
            hpElecsP4 += hpElecP4[i];

        hardProcLeptonsP4 = hpMuonsP4 + hpElecsP4;


        // Final state
        TLorentzVector fsMuonsP4, fsElecsP4;

        for (unsigned i = 0; i < nFinalStateMuons; i++)
            fsMuonsP4 += fsMuonP4[i];
        for (unsigned i = 0; i < nFinalStateElectrons; i++)
            fsElecsP4 += fsElecP4[i];

        finalStateLeptonsP4 = fsMuonsP4 + fsElecsP4;
    }


    // Lumi mask
    if (isData)
    {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);



    ////  TRIGGER SELECTION

    evtMuonTriggered = kFALSE;
    for (unsigned i = 0; i < muonTriggerNames.size(); i++)
    {
        if (trigger->pass(muonTriggerNames[i], fInfo->triggerBits))
            evtMuonTriggered = kTRUE;
    }

    evtElectronTriggered = kFALSE;
    for (unsigned i = 0; i < electronTriggerNames.size(); i++)
    {
        if (trigger->pass(electronTriggerNames[i], fInfo->triggerBits))
            evtElectronTriggered = kTRUE;
    }

    bool passTrigger = evtMuonTriggered || evtElectronTriggered;

    if (!passTrigger)
        return kTRUE;
    hTotalEvents->Fill(3);




    ////////////////////
    // SELECT OBJECTS //
    ////////////////////


    ////  VERTICES

    if (fInfo->hasGoodPV)
    {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
//      particleSelector->SetPV(pv);
    }
    else
        return kTRUE;

    hTotalEvents->Fill(4);
//  particleSelector->SetNPV(fInfo->nPU + 1);
//  particleSelector->SetRho(fInfo->rhoJet);



    ////  MUONS

    vector<TMuon*> muons;
    for (int i = 0; i < fMuonArr->GetEntries(); i++)
    {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);
        TLorentzVector muonP4;
        copy_p4(muon, MUON_MASS, muonP4);

        if (muon->pt > 4)
        {
            nMuons++;
            muons.push_back(muon);
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);


    // Second pass: store P4 & SF, ghost clean, identify loose & HZZ muons
    for (unsigned i = 0; i < nMuons; i++)
    {
        TMuon* muon = muons[i];
        TLorentzVector muonP4;
        copy_p4(muons[i], MUON_MASS, muonP4);


        // Store P4 before corrections
        new(muonsP4ptr[i]) TLorentzVector(muonP4);


        // Apply Rochester correction
        muonEnergySF.push_back(GetRochesterCorrection(muon, muonCorr, rng, isData));
        muon->pt *= muonEnergySF.back();
        muonP4.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, MUON_MASS);


        // Ghost cleaning
        // Does not include track segment matching...
        float minDeltaR = 0.02;
        bool isGhost = kFALSE;
        for (unsigned j = 0; j < i; j++)
        {
            TLorentzVector tmpMuonP4;
            tmpMuonP4.SetPtEtaPhiM(muons[j]->pt, muons[j]->eta, muons[j]->phi, MUON_MASS);
            float dr = muonP4.DeltaR(tmpMuonP4);
            if (dr < minDeltaR)
            {
                isGhost = kTRUE;
                break;
            }
        }
        muonIsGhost.push_back(isGhost);


        // Remove muon track pt from muon track isolation variable
        for (unsigned j = i + 1; j < nMuons; j++)
        {
            TLorentzVector muon_j;
            copy_p4(muons[j], MUON_MASS, muon_j);

            if (muonP4.DeltaR(muon_j) < 0.3)
            {
                muon->trkIso = max(0., muon->trkIso - muon_j.Pt());
                muons[j]->trkIso = max(0., muons[j]->trkIso - muonP4.Pt());
            }
        }


        // Check muon ID
        if (muon->pogIDBits & baconhep::kPOGLooseMuon)
            nLooseMuons++;
        if (PassMuonHZZTightID(muon))
            nHZZMuons++;
    }



    ////  ELECTRONS

    vector<TElectron*> electrons;
    for (int i = 0; i < fElectronArr->GetEntries(); i++)
    {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (electron->pt > 6)
        {
            nElectrons++;
            electrons.push_back(electron);
        }
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);


    // Second pass: store P4 & SF; identify loose, HZZ, & tight electrons
    for (unsigned i = 0; i < nElectrons; i++)
    {
        TElectron* electron = electrons[i];
        TLorentzVector electronP4;
        copy_p4(electrons[i], ELE_MASS, electronP4);

        // Store P4 before corrections
        new(electronsP4ptr[i]) TLorentzVector(electronP4);

        // Apply electron energy correction
        electronEnergySF.push_back(GetElectronPtSF(electron));
        electronP4 *= electronEnergySF.back();


        // Cross cleaning
        float minDeltaR = 0.05;
        bool isGhost = kFALSE;
        for (unsigned j = 0; j < nMuons; j++)
        {
            TMuon* tmpMuon = muons[j];
            TLorentzVector tmpMuonP4;
            copy_p4(tmpMuon, MUON_MASS, tmpMuonP4);

            float dr = electronP4.DeltaR(tmpMuonP4);
            if (dr < minDeltaR && PassMuonHZZTightID(tmpMuon))
            {
                isGhost = kTRUE;
                break;
            }
        }
        electronIsGhost.push_back(isGhost);


        // Check electron ID
        if (electron->pass2016HZZwpLoose)
            nV1Electrons++;
        if (electron->pass2017isoV2wpHZZ)
            nV2Electrons++;
        if (PassElectronHZZTightID(electron, fInfo->rhoJet))
            nHZZElectrons++;
    }

    nLeptons        = nMuons + nElectrons;
    nHZZLeptons     = nHZZMuons + nHZZElectrons;




    ///////////////
    // SELECTION //
    ///////////////

    // Require at least two same-flavor leptons passing HZZ cuts
    if (nHZZMuons < 2 && nHZZElectrons < 2)
        return kTRUE;
    hTotalEvents->Fill(5);



    /////////////////////
    // FILL CONTAINERS //
    /////////////////////


    ////  MUONS

    for (unsigned i = 0; i < nMuons; i++)
    {
        // Basic quantities
        TMuon* muon = muons[i];
        TLorentzVector muonP4;
        copy_p4(muons[i], MUON_MASS, muonP4);

        muonsQ.push_back(muon->q);


        // Boolean ID
        muonIsLoose.push_back(muon->pogIDBits & baconhep::kPOGLooseMuon);
        muonIsHZZ.push_back(PassMuonHZZTightID(muon));

        EfficiencyContainer effCont;
        pair<float, float> trigEff;

        effCont = weights->GetHZZMuonIDEff(muon);
        muonHZZIDSF.push_back(effCont.GetSF());


        // Trigger bools and SFs
        bool firedLeg1 = kFALSE, firedLeg2 = kFALSE;
        for (unsigned i = 0; i < muonTriggerNames.size(); i++)
        {
            if (trigger->passObj(muonTriggerNames[i], 1, muon->hltMatchBits))
                firedLeg1 = kTRUE;
            if (trigger->passObj(muonTriggerNames[i], 2, muon->hltMatchBits))
                firedLeg2 = kTRUE;
        }
        muonFiredLeg1.push_back(firedLeg1);
        muonFiredLeg2.push_back(firedLeg2);

        // Leptons will automatically fail leg 2 for single trigger selection
//      if (params->selection == "single")
//          effCont = weights->GetSingleMuonTriggerEff(muon);
//      else if (params->selection == "double")
            effCont = weights->GetDoubleMuonTriggerEff(muon, 1);
        trigEff = effCont.GetEff();

        if (!firedLeg1)
        {
            trigEff.first = 0;
            trigEff.second = 0;
        }
        muonTrigEffLeg1Data.push_back(trigEff.first);
        muonTrigEffLeg1MC.push_back(trigEff.second);

        effCont = weights->GetDoubleMuonTriggerEff(muon, 2);
        trigEff = effCont.GetEff();
        if (!firedLeg2)
        {
            trigEff.first = 0;
            trigEff.second = 0;
        }
        muonTrigEffLeg2Data.push_back(trigEff.first);
        muonTrigEffLeg2MC.push_back(trigEff.second);


        // ID criteria for downstream use
        muonCombIso.push_back(GetMuonIsolation(muon));
        muonsTrkIso.push_back(muon->trkIso);
        muonD0.push_back(muon->d0);
        muonDz.push_back(muon->dz);
        muonSIP3d.push_back(muon->sip3d);
        muonPtErr.push_back(muon->ptErr);
        muonNMatchStn.push_back(muon->nMatchStn);
        muonNPixHits.push_back(muon->nPixHits);
        muonNTkLayers.push_back(muon->nTkLayers);
        muonIsPF.push_back(muon->typeBits & baconhep::kPFMuon);
        muonIsGlobal.push_back(muon->typeBits & baconhep::kGlobal);
        muonIsTracker.push_back(muon->typeBits & baconhep::kTracker);
        muonBestTrackType.push_back(muon->btt);
    }



    ////  ELECTRONS

    for (unsigned i = 0; i < nElectrons; i++)
    {
        // Basic quantities
        TElectron* electron = electrons[i];
        TLorentzVector electronP4;
        copy_p4(electrons[i], ELE_MASS, electronP4);

        electronsQ.push_back(electron->q);

        // Boolean ID
        electronPassV1MVA.push_back(electron->pass2016HZZwpLoose);
        electronPassV2MVA.push_back(electron->pass2017isoV2wpHZZ);
        electronIsHZZ.push_back(PassElectronHZZTightID(electron, fInfo->rhoJet));

        EfficiencyContainer effCont;
        pair<float, float> trigEff;

        effCont = weights->GetHZZElectronIDRecoEff(electron);
        electronHZZIDSF.push_back(effCont.GetSF());


        // Trigger
        bool firedLeg1 = kFALSE, firedLeg2 = kFALSE;
        for (unsigned i = 0; i < electronTriggerNames.size(); i++)
        {
            if (trigger->passObj(electronTriggerNames[i], 1, electron->hltMatchBits))
                firedLeg1 = kTRUE;
            if (trigger->passObj(electronTriggerNames[i], 2, electron->hltMatchBits))
                firedLeg2 = kTRUE;
        }
        electronFiredLeg1.push_back(firedLeg1);
        electronFiredLeg2.push_back(firedLeg2);

        effCont = weights->GetDoubleElectronTriggerEff(electron, 1);
        trigEff = effCont.GetEff();
        if (!firedLeg1)
        {
            trigEff.first = 0;
            trigEff.second = 0;
        }
        electronTrigEffLeg1Data.push_back(trigEff.first);
        electronTrigEffLeg1MC.push_back(trigEff.second);

        effCont = weights->GetDoubleElectronTriggerEff(electron, 2);
        trigEff = effCont.GetEff();
        if (!firedLeg2)
        {
            trigEff.first = 0;
            trigEff.second = 0;
        }
        electronTrigEffLeg2Data.push_back(trigEff.first);
        electronTrigEffLeg2MC.push_back(trigEff.second);


        // ID criteria
        electronCombIso.push_back(GetElectronIsolation(electron, fInfo->rhoJet));
        electronsTrkIso.push_back(electron->trkIso);
        electronD0.push_back(electron->d0);
        electronDz.push_back(electron->dz);
        electronSIP3d.push_back(electron->sip3d);
        electronScEta.push_back(electron->scEta);
        electronNMissHits.push_back(electron->nMissingHits);
        electronIsGap.push_back(electron->fiducialBits & kIsGap);
    }



    ////  GEN PARTICLES

    if (isSignal)
    {
        // Final state
        for (unsigned i = 0; i < nFinalStateMuons; i++)
        {
            new(finalStateMuonP4ptr[i]) TLorentzVector(fsMuonP4[i]);
            finalStateMuonQ.push_back(fsMuonQ[i]);
            finalStateMuonZIndex.push_back(fsMuonZIndex[i]);
        }
        for (unsigned i = 0; i < nFinalStateElectrons; i++)
        {
            new(finalStateElectronP4ptr[i]) TLorentzVector(fsElecP4[i]);
            finalStateElectronQ.push_back(fsElecQ[i]);
            finalStateElectronZIndex.push_back(fsElecZIndex[i]);
        }


        // Hard process
        for (unsigned i = 0; i < nHardProcMuons; i++)
        {
            new(hardProcMuonP4ptr[i]) TLorentzVector(hpMuonP4[i]);
            hardProcMuonQ.push_back(hpMuonQ[i]);
            hardProcMuonZIndex.push_back(hpMuonZIndex[i]);
        }
        for (unsigned i = 0; i < nHardProcElectrons; i++)
        {
            new(hardProcElectronP4ptr[i]) TLorentzVector(hpElecP4[i]);
            hardProcElectronQ.push_back(hpElecQ[i]);
            hardProcElectronZIndex.push_back(hpElecZIndex[i]);
        }
    }


    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void MultileptonAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void MultileptonAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void MultileptonAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<MultileptonAnalyzer> selector(new MultileptonAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}



// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Muons
float MultileptonAnalyzer::GetMuonIsolation(const baconhep::TMuon* mu)
{
    float combIso = (mu->chHadIso03 + std::max(0.,(double)mu->neuHadIso03 + mu->gammaIso03 - 0.5*mu->puIso03));
    return combIso;
}



// roccor.Run2.v1.tgz: https://twiki.cern.ch/twiki/bin/view/CMS/RochcorMuon
// From README:
// double mcSF = rc.kSpreadMC(Q, pt, eta, phi, genPt, s=0, m=0); //(recommended), MC scale and resolution correction when matched gen muon is available
// double mcSF = rc.kSmearMC(Q, pt, eta, phi, nl, u, s=0, m=0); //MC scale and extra smearing when matched gen muon is not available
float MultileptonAnalyzer::GetRochesterCorrection(const baconhep::TMuon* muon, RoccoR* muonCorr, TRandom3* rng, bool isData)
{   
    // https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Muon_scale_and_resolution_correc
    if (muon->pt < 200 && muon->btt == 1 && muon->nTkLayers > 5)
    {
        if (isData)
            return muonCorr->kScaleDT(muon->q, muon->pt, muon->eta, muon->phi, 0, 0);
        else
            return muonCorr->kSmearMC(muon->q, muon->pt, muon->eta, muon->phi, muon->nTkLayers, rng->Rndm(), 0, 0);
    }
    else
        return 1;
}



// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Muons
// No change since 2016
bool MultileptonAnalyzer::PassMuonHZZTightID(const baconhep::TMuon* muon)
{
    // AFTER Rochester correction!

    if (    muon->pt > 5
        &&  fabs(muon->eta) < 2.4
        &&  fabs(muon->d0) < 0.5        // These are, indeed, corrected: https://github.com/NWUHEP/BaconProd/blob/jbueghly_2017/Ntupler/src/FillerMuon.cc
        &&  fabs(muon->dz) < 1.0
        &&  ((muon->typeBits & baconhep::kGlobal) || ((muon->typeBits & baconhep::kTracker) && muon->nMatchStn > 0)) // Global muon or (arbitrated) tracker muon
        &&  muon->btt != 2
        &&  fabs(muon->sip3d) < 4.0                 
        &&  GetMuonIsolation(muon)/muon->pt < 0.35)
    {
        if (muon->pogIDBits & baconhep::kPOGLooseMuon)
            return kTRUE;
        else if (   muon->pt > 200
                 && muon->nMatchStn > 1
                 && (muon->ptErr/muon->pt) < 0.3
                 && fabs(muon->d0) < 0.2
                 && fabs(muon->dz) < 0.5
                 && muon->nPixHits > 0
                 && muon->nTkLayers > 5)
            return kTRUE;
        else
            return kFALSE;
    }
    else
        return kFALSE;
}


// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Electrons
// seems to have typos?  So this is unchanged from 2016
//
// Effective area from
// https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt
// and
// https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_92X.txt
// which is NOT the page linked on the HZZ twiki!!!
float MultileptonAnalyzer::GetElectronIsolation(const baconhep::TElectron* el, float rho)
{
    int iEta = 0;
    float etaBins[8] = {0., 1., 1.479, 2.0, 2.2, 2.3, 2.4, 5.0};
    float effArea[7] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};
//  float effArea2017[7] = {0.1566, 0.1626, 0.1073, 0.0854, 0.1051, 0.1204, 0.1524};
    for (unsigned i = 0; i < 7; i++)
    {
        if (fabs(el->scEta) > etaBins[i] && fabs(el->scEta) < etaBins[i+1])
            iEta = i;
    }

    // These seem to be defined in CMSSW for 0.3 iso cone (see HZZ twiki) and they are pulled directly in the ntupler
    float combIso = el->chHadIso + std::max(0., (double)el->neuHadIso + el->gammaIso - rho*effArea[iEta]);
    return combIso;
}


// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaMiniAODV2#Energy_Scale_and_Smearing
float MultileptonAnalyzer::GetElectronPtSF(baconhep::TElectron* electron)
{
    return electron->ecalTrkEnergyPostCorr / electron->energy;
}


// https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Electrons
bool MultileptonAnalyzer::PassElectronHZZTightID(const baconhep::TElectron* electron, float rho)
{
    // FIXME to use corrected pt
    if (    electron->pt > 7
        &&  fabs(electron->scEta) < 2.5
//      &&  (particleSelector->PassElectronMVA(electron, cuts->hzzIsoV1) || particleSelector->PassElectronMVA(electron, cuts->hzzMVANoIsoV1))   // Needs to be checked separately!!
        &&  fabs(electron->d0) < 0.5
        &&  fabs(electron->dz) < 1
        &&  fabs(electron->sip3d) < 4.0
        &&  GetElectronIsolation(electron, rho)/electron->pt < 0.35)   // Not sure if this needs to be included for iso MVA?
        return kTRUE;
    else
        return kFALSE;
}
