#include "MultileptonAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 


MultileptonAnalyzer::MultileptonAnalyzer() : BLTSelector()
{

}

MultileptonAnalyzer::~MultileptonAnalyzer()
{

}

void MultileptonAnalyzer::Begin(TTree *tree)
{
    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
              std::sregex_token_iterator(), std::back_inserter(options));

    const std::string cmssw_base = getenv("CMSSW_BASE");
    const std::string data_dir = "/src/BLT/BLTAnalysis/data/";



    //
    //  SETUP
    //

    // Parameters
    params.reset(new Parameters());
    params->setup(options);

    const bool isSignal     = params->dataset == "ZZTo4L" || params->dataset == "ZZTo4L_aMC" || params->dataset == "ZZTo4L_M-1";
    const bool isDrellYan   = params->dataset == "DYJetsToLL_M-50";

    // Particle selector, cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Weight utilities
    weights.reset(new WeightUtils(params->period, params->selection, false));

    // Lumi mask
    lumiMask = RunLumiRangeMap();
    string jsonFileName;
    if      (params->period == "2016")
        jsonFileName = "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt";
    else if (params->period == "2017")
        jsonFileName = "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt";
    else if (params->period == "2018")
        jsonFileName = "Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt";
    lumiMask.AddJSONFile(cmssw_base + data_dir + jsonFileName);


    // Triggers
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if      (params->period == "2016")
    {
        doubleMuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        doubleMuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        doubleMuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        doubleMuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
        doubleMuonTriggers.push_back("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        doubleMuonTriggers.push_back("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
        singleMuonTriggers.push_back("HLT_IsoMu24_v*");
        singleMuonTriggers.push_back("HLT_IsoTkMu24_v*");

        doubleElecTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
        singleElecTriggers.push_back("HLT_Ele27_WPTight_Gsf_v*");

        MUON_LEG1_PT = 17;      MUON_LEG2_PT = 8;       MUON_SINGLE_PT = 24;
        ELEC_LEG1_PT = 23;      ELEC_LEG2_PT = 12;      ELEC_SINGLE_PT = 27;
    }

    else if (params->period == "2017")
    {
        doubleMuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
        doubleMuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*");
        singleMuonTriggers.push_back("HLT_IsoMu27_v*");

        doubleElecTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
        singleElecTriggers.push_back("HLT_Ele35_WPTight_Gsf_v*");

        MUON_LEG1_PT = 17;      MUON_LEG2_PT = 8;       MUON_SINGLE_PT = 27;
        ELEC_LEG1_PT = 23;      ELEC_LEG2_PT = 12;      ELEC_SINGLE_PT = 35;
    }

    else if (params->period == "2018")
    {
        doubleMuonTriggers.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
        singleMuonTriggers.push_back("HLT_IsoMu24_v*");

        doubleElecTriggers.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
        singleElecTriggers.push_back("HLT_Ele35_WPTight_Gsf_v*");

        MUON_LEG1_PT = 17;      MUON_LEG2_PT = 8;       MUON_SINGLE_PT = 24;
        ELEC_LEG1_PT = 23;      ELEC_LEG2_PT = 12;      ELEC_SINGLE_PT = 35;
    }



    //
    //  OUTPUT TREE
    //

    string outFileName = params->get_output_filename("output");
    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    string outTreeName = params->get_output_treename("tree");
    outTree = new TTree(outTreeName.c_str(), "bltTree");


    // Branches
    outTree->Branch(    "runNumber",                &runNumber);
    outTree->Branch(    "evtNumber",                &evtNumber);
    outTree->Branch(    "lumiSection",              &lumiSection);
    outTree->Branch(    "evtMuonTriggered",         &evtMuonTriggered);
    outTree->Branch(    "evtDoubleMuTriggered",     &evtDoubleMuTriggered);
    outTree->Branch(    "evtSingleMuTriggered",     &evtSingleMuTriggered);
    outTree->Branch(    "evtElectronTriggered",     &evtElectronTriggered);
    outTree->Branch(    "evtDoubleElTriggered",     &evtDoubleElTriggered);
    outTree->Branch(    "evtSingleElTriggered",     &evtSingleElTriggered);

    outTree->Branch(    "genWeight",                &genWeight);
    outTree->Branch(    "ECALWeight",               &ECALWeight);
    outTree->Branch(    "ECALWeightUp",             &ECALWeightUp);
    outTree->Branch(    "ECALWeightDown",           &ECALWeightDown);
    outTree->Branch(    "PUWeight",                 &PUWeight);
    outTree->Branch(    "nPU",                      &nPU);
    outTree->Branch(    "nPV",                      &nPV);
    outTree->Branch(    "hasTauDecay",              &hasTauDecay);

    outTree->Branch(    "nLooseMuons",              &nLooseMuons);
    outTree->Branch(    "nLooseElectrons",          &nLooseElectrons);
    outTree->Branch(    "nLooseLeptons",            &nLooseLeptons);
    outTree->Branch(    "nTightMuons",              &nTightMuons);
    outTree->Branch(    "nTightElectrons",          &nTightElectrons);
    outTree->Branch(    "nTightLeptons",            &nTightLeptons);

    outTree->Branch(    "muonP4",                   &muonP4_,               32000,      1);
    outTree->Branch(    "muonUncorrectedP4",        &muonUncorrP4_,         32000,      1);
    outTree->Branch(    "muonQ",                    &muonCharge);
    outTree->Branch(    "muonEnergySF",             &muonEnergySF);
    outTree->Branch(    "muonEnergySFUp",           &muonEnergySFUp);
    outTree->Branch(    "muonEnergySFDown",         &muonEnergySFDown);
    outTree->Branch(    "muonIDSF",                 &muonIDSF);
    outTree->Branch(    "muonIsolation",            &muonIsolation);
    outTree->Branch(    "muonIsTight",              &muonIsTight);
    outTree->Branch(    "muonIsLoose",              &muonIsLoose);
    outTree->Branch(    "muonIsIsolated",           &muonIsIsolated);
    outTree->Branch(    "muonIsPF",                 &muonIsPF);
    outTree->Branch(    "muonIsTrackerHighPt",      &muonIsTrackerHighPt);
    outTree->Branch(    "muonFiredLeg1",            &muonFiredLeg1);
    outTree->Branch(    "muonFiredLeg2",            &muonFiredLeg2);
    outTree->Branch(    "muonFiredSingle",          &muonFiredSingle);

    outTree->Branch(    "electronP4",               &electronP4_,           32000,      1);
    outTree->Branch(    "electronUncorrectedP4",    &electronUncorrP4_,     32000,      1);
    outTree->Branch(    "electronQ",                &electronCharge);
    outTree->Branch(    "electronEnergySF",         &electronEnergySF);
    outTree->Branch(    "electronEnergySFUp",       &electronEnergySFUp);
    outTree->Branch(    "electronEnergySFDown",     &electronEnergySFDown);
    outTree->Branch(    "electronIDSF",             &electronIDSF);
    outTree->Branch(    "electronRecoSF",           &electronRecoSF);
    outTree->Branch(    "electronIsolation",        &electronIsolation);
    outTree->Branch(    "electronScEta",            &electronScEta);
    outTree->Branch(    "electronScEt",             &electronScEt);
    outTree->Branch(    "electronIsTight",          &electronIsTight);
    outTree->Branch(    "electronIsLoose",          &electronIsLoose);
    outTree->Branch(    "electronIsV2Iso",          &electronIsV2Iso);
    outTree->Branch(    "electronIsGap",            &electronIsGap);
    outTree->Branch(    "electronFiredLeg1",        &electronFiredLeg1);
    outTree->Branch(    "electronFiredLeg2",        &electronFiredLeg2);
    outTree->Branch(    "electronFiredSingle",      &electronFiredSingle);

    if (isSignal || isDrellYan)
    {
        outTree->Branch(    "nDressedMuons",            &nDressedMuons);
        outTree->Branch(    "nDressedElectrons",        &nDressedElectrons);
        outTree->Branch(    "nDressedLeptons",          &nDressedLeptons);
        outTree->Branch(    "dressedLeptonsP4",         &dressedLeptonsP4);
        outTree->Branch(    "dressedMuonP4",            &dressedMuonP4_,        32000,      1);
        outTree->Branch(    "dressedMuonQ",             &dressedMuonQ);
        outTree->Branch(    "dressedMuonZIndex",        &dressedMuonZIndex);
        outTree->Branch(    "dressedElectronP4",        &dressedElectronP4_,    32000,      1);
        outTree->Branch(    "dressedElectronQ",         &dressedElectronQ);
        outTree->Branch(    "dressedElectronZIndex",    &dressedElectronZIndex);
    }

    // Histograms
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);


    ReportPostBegin();
}


Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{
    //
    //  CLEAR CONTAINERS
    //
  
    evtMuonTriggered = kFALSE;      evtDoubleMuTriggered = kFALSE;  evtSingleMuTriggered = kFALSE;
    evtElectronTriggered = kFALSE;  evtDoubleElTriggered = kFALSE;  evtSingleElTriggered = kFALSE;

    nLooseMuons     = 0;            nLooseElectrons = 0;            nLooseLeptons   = 0;
    nTightMuons     = 0;            nTightElectrons = 0;            nTightLeptons   = 0; 
    genWeight       = 1;            PUWeight        = 1;            nPU             = 0;
    ECALWeight      = 1;            ECALWeightUp    = 1;            ECALWeightDown  = 1;
    hasTauDecay     = kFALSE;

    muonP4_->Delete();              muonUncorrP4_->Delete();        muonCharge.clear();
    muonEnergySF.clear();           muonEnergySFUp.clear();         muonEnergySFDown.clear();
    muonIDSF.clear();               muonIsolation.clear();
    muonIsTight.clear();            muonIsLoose.clear();            muonIsIsolated.clear();
    muonIsPF.clear();               muonIsTrackerHighPt.clear();
    muonFiredLeg1.clear();          muonFiredLeg2.clear();          muonFiredSingle.clear();
                                                                                                        
    electronP4_->Delete();          electronUncorrP4_->Delete();    electronCharge.clear();
    electronEnergySF.clear();       electronEnergySFUp.clear();     electronEnergySFDown.clear();
    electronIDSF.clear();           electronRecoSF.clear();         electronIsolation.clear();
    electronScEta.clear();          electronScEt.clear();           electronIsGap.clear();
    electronIsTight.clear();        electronIsLoose.clear();        electronIsV2Iso.clear();
    electronFiredLeg1.clear();      electronFiredLeg2.clear();      electronFiredSingle.clear();

    nDressedMuons = 0;              nDressedElectrons = 0;          nDressedLeptons = 0; 
    dressedMuonP4_->Delete();       dressedMuonQ.clear();           dressedMuonZIndex.clear();
    dressedElectronP4_->Delete();   dressedElectronQ.clear();       dressedElectronZIndex.clear();
    dressedLeptonsP4.Clear();






    ////
    ////
    ////    START
    ////
    ////
    

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry << " Run: " << fInfo->runNum 
                  << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum 
                  << std::endl;

    // Parameters
    const bool isData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isData);
    weights->SetDataBit(isData);

    TString sampleName = params->dataset;
 
    const bool isSignal     = sampleName.EqualTo("ZZTo4L") || sampleName.EqualTo("ZZTo4L_aMC") || sampleName.EqualTo("ZZTo4L_M-1");
    const bool isDrellYan   = sampleName.Contains("DYJetsToLL_M-50");



    //
    //  EVENT INFO
    //

    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;
    nPV         = fPVArr->GetEntries();

    if (!isData)
    {
        // Gen weight
        genWeight = fGenEvtInfo->weight > 0 ? 1 : -1; 
        if (genWeight < 0)
            hTotalEvents->Fill(10);

        // Pileup weight
        nPU = fInfo->nPUmean;
        PUWeight = weights->GetPUWeight(nPU);

        // ECAL weight
        ECALWeight      = fInfo->ecalWeight;
        ECALWeightUp    = fInfo->ecalWeightUp;
        ECALWeightDown  = fInfo->ecalWeightDown;
    }
    else
    {
        // Lumi mask
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);



    //
    //  GEN PARTICLES
    //

    if (isSignal || isDrellYan)
    {
        dressedLeptonsP4 = TLorentzVector();

        // Particle loop
        for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
            int motherIndex = particle->parent;

            if (motherIndex < 0)            // seg faults are bad!
                continue;

            TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(motherIndex);

            // Look for taus from a Z
            if ((abs(particle->pdgId) == 15) && (mother->pdgId == 23))
                hasTauDecay = kTRUE;

            if (!isSignal && !isDrellYan)   // don't waste time if we aren't storing the leptons
                continue;

            // Now look for electrons and muons
            if ((abs(particle->pdgId) != 13) && (abs(particle->pdgId) != 11))
                continue;
            if (particle->status != 1)      // no point in saving Born leptons anymore...
                continue;

            // Try to trace back to a Z
            while (((abs(mother->pdgId) == 11) || (abs(mother->pdgId) == 13)) && (mother->parent >= 0))
            {
                motherIndex = mother->parent;
                mother = (TGenParticle*) fGenParticleArr->At(motherIndex);
            }

            if (mother->pdgId != 23)
                continue;

            // FSR recovery
            TLorentzVector lepP4;
            copy_p4(particle, lepP4);

            for (int j = 0; j < fGenParticleArr->GetEntries(); j++)
            {
                TGenParticle* gamma = (TGenParticle*) fGenParticleArr->At(i);

                if ((gamma->pdgId != 22) || (gamma->status != 1))
                    continue;

                TLorentzVector gammaP4;
                copy_p4(gamma, 0, gammaP4);

                if (lepP4.DeltaR(gammaP4) < 0.1)
                    lepP4 = lepP4 + gammaP4;
            }

            particle->pt = lepP4.Pt();
            particle->eta = lepP4.Eta();
            particle->phi = lepP4.Phi();
            particle->mass = lepP4.M();

            int charge = -1 * copysign(1, particle->pdgId);

            if      (abs(particle->pdgId) == 13)
            {
                TLorentzVector *p4 = (TLorentzVector*)dressedMuonP4_->ConstructedAt(nDressedMuons);
                copy_p4(particle, MUON_MASS, p4);
                dressedLeptonsP4 = dressedLeptonsP4 + *p4;
                dressedMuonQ.push_back(charge);
                dressedMuonZIndex.push_back(motherIndex);
                nDressedMuons++;
            }
            else if (abs(particle->pdgId) == 11)
            {
                TLorentzVector *p4 = (TLorentzVector*)dressedElectronP4_->ConstructedAt(nDressedElectrons);
                copy_p4(particle, ELE_MASS, p4);
                dressedLeptonsP4 = dressedLeptonsP4 + *p4;
                dressedElectronQ.push_back(charge);
                dressedElectronZIndex.push_back(motherIndex);
                nDressedElectrons++;
            }
        }
        nDressedLeptons = nDressedMuons + nDressedElectrons;
    }



    //
    //  TRIGGER SELECTION
    //

    std::vector<string> firedDoubleMuTriggers, firedDoubleElTriggers;
    std::vector<string> firedSingleMuTriggers, firedSingleElTriggers;

    for (unsigned i = 0; i < doubleMuonTriggers.size(); i++)            // double muon
    {
        if (trigger->pass(doubleMuonTriggers[i], fInfo->triggerBits))
            firedDoubleMuTriggers.push_back(doubleMuonTriggers[i]);
    }
    evtDoubleMuTriggered = firedDoubleMuTriggers.size() > 0;

    for (unsigned i = 0; i < singleMuonTriggers.size(); i++)            // single muon
    {
        if (trigger->pass(singleMuonTriggers[i], fInfo->triggerBits))
            firedSingleMuTriggers.push_back(singleMuonTriggers[i]);
    }
    evtSingleMuTriggered = firedSingleMuTriggers.size() > 0;

    for (unsigned i = 0; i < doubleElecTriggers.size(); i++)            // double electron
    {
        if (trigger->pass(doubleElecTriggers[i], fInfo->triggerBits))
            firedDoubleElTriggers.push_back(doubleElecTriggers[i]);
    }
    evtDoubleElTriggered = firedDoubleElTriggers.size() > 0;

    for (unsigned i = 0; i < singleElecTriggers.size(); i++)            // single electron
    {
        if (trigger->pass(singleElecTriggers[i], fInfo->triggerBits))
            firedSingleElTriggers.push_back(singleElecTriggers[i]);
    }
    evtSingleElTriggered = firedSingleElTriggers.size() > 0;
    
    evtMuonTriggered = evtSingleMuTriggered || evtDoubleMuTriggered;
    evtElectronTriggered = evtSingleElTriggered || evtDoubleElTriggered;

    bool passTrigger = evtMuonTriggered || evtElectronTriggered;

    if (!passTrigger && !isSignal && !isDrellYan)
        return kTRUE;
    hTotalEvents->Fill(3);



    //
    //  OBJECTS
    //

    // Vertices
    if (fInfo->hasGoodPV)
    {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    }
    else
        return kTRUE;
    hTotalEvents->Fill(4);

    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);


    // Muons
    vector<TMuon*> muons;
    for (int i = 0; i < fMuonArr->GetEntries(); i++)
    {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        // Store only (preliminary) loose muons
        if (particleSelector->PassMuonID(muon, cuts->looseHZZMuonID))
            muons.push_back(muon);
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    nLooseMuons = muons.size();


    // Electrons
    vector<TElectron*> electrons;
    for (int i = 0; i < fElectronArr->GetEntries(); i++)
    {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        // Store only (preliminary) loose electrons
        if (particleSelector->PassElectronID(electron, cuts->looseHZZElectronID))
            electrons.push_back(electron);
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);
    nLooseElectrons = electrons.size();






    ////
    ////
    ////    PRESELECTION
    ////
    ////


    if (nLooseMuons < 2 && nLooseElectrons < 2)
        return kTRUE;






    ////
    ////
    ////    CONTAINERS
    ////
    ////


    nLooseMuons = 0;
    nLooseElectrons = 0;


    //
    //  MUONS
    //

    for (unsigned i = 0; i < muons.size(); i++)
    {
        TMuon* muon = muons[i];

        // Initialize entries in TClonesArrays
        TLorentzVector *uncorrP4 = (TLorentzVector*)muonUncorrP4_->ConstructedAt(i);
        TLorentzVector *corrP4 = (TLorentzVector*)muonP4_->ConstructedAt(i);

        // Uncorrected momentum
        copy_p4(muon, MUON_MASS, uncorrP4);

        // Charge
        muonCharge.push_back(muon->q);

        // Rochester correction
        muonEnergySFUp.push_back(particleSelector->GetRochesterCorrection(muon, "up"));
        muonEnergySFDown.push_back(particleSelector->GetRochesterCorrection(muon, "down"));

        float corr = particleSelector->GetRochesterCorrection(muon);
        muon->pt *= corr;
        copy_p4(muon, MUON_MASS, corrP4);
        muonEnergySF.push_back(corr);

        // ID scale factor
        muonIDSF.push_back(weights->GetHZZMuonIDSF(muon));

        // Isolation
        muonIsolation.push_back(particleSelector->GetMuonIso(muon));

        // ID/iso bools
        muonIsTight.push_back(particleSelector->PassMuonID(muon, cuts->tightHZZMuonID));
        muonIsLoose.push_back(particleSelector->PassMuonID(muon, cuts->looseHZZMuonID));
        muonIsIsolated.push_back(particleSelector->PassMuonIso(muon, cuts->wpHZZMuonIso));
        muonIsPF.push_back(muon->typeBits & baconhep::kPFMuon);
        muonIsTrackerHighPt.push_back(particleSelector->PassMuonID(muon, cuts->trackerHighPtMuonID));

        if (muonIsLoose.back())
            nLooseMuons++;
        if (muonIsTight.back())
            nTightMuons++;
    }



    //
    //  ELECTRONS
    //

    for (unsigned i = 0; i < electrons.size(); i++)
    {
        TElectron* electron = electrons[i];

        // Initialize entries in TClonesArrays
        TLorentzVector *uncorrP4 = (TLorentzVector*)electronUncorrP4_->ConstructedAt(i);
        TLorentzVector *corrP4 = (TLorentzVector*)electronP4_->ConstructedAt(i);

        // Uncorrected momentum
        copy_p4(electrons[i], ELE_MASS, uncorrP4);
        electronScEta.push_back(electron->scEta);
        electronScEt.push_back(electron->scEt);

        // Charge
        electronCharge.push_back(electron->q);

        // Energy correction
        electronEnergySFUp.push_back(particleSelector->GetElectronCorrection(electron, "up"));
        electronEnergySFDown.push_back(particleSelector->GetElectronCorrection(electron, "down"));

        float corr = particleSelector->GetElectronCorrection(electron);
        TLorentzVector electronP4(*uncorrP4);
        electronP4 *= corr;

        electron->pt = electronP4.Pt();
        electron->eta = electronP4.Eta();
        electron->phi = electronP4.Phi();
        copy_p4(electrons[i], ELE_MASS, corrP4);

        electronEnergySF.push_back(corr);

        // ID scale factor
        electronIDSF.push_back(weights->GetHZZElectronIDSF(electron));
        electronRecoSF.push_back(weights->GetElectronRecoSF(electron));

        // Isolation
        electronIsolation.push_back(particleSelector->GetElectronIso(electron));

        // ID/iso bools
        electronIsTight.push_back(particleSelector->PassElectronID(electron, cuts->tightHZZElectronID));
        electronIsLoose.push_back(particleSelector->PassElectronID(electron, cuts->looseHZZElectronID));
        electronIsV2Iso.push_back(electron->pass2017isoV2wpHZZ);
        electronIsGap.push_back(electron->fiducialBits & kIsGap);

        if (electronIsLoose.back())
            nLooseElectrons++;
        if (electronIsTight.back())
            nTightElectrons++;
    }

    nLooseLeptons = nLooseMuons + nLooseElectrons;
    nTightLeptons = nTightMuons + nTightElectrons;





    ////
    ////
    ////    TRIGGER MATCHING
    ////
    ////


    std::vector<string> matchedDoubleMuTriggers, matchedDoubleElTriggers;
    std::vector<string> matchedSingleMuTriggers, matchedSingleElTriggers;



    //
    //  DOUBLE MUON
    //

    for (unsigned i = 0; i < firedDoubleMuTriggers.size(); i++)
    {
        bool matched = kFALSE;
        for (unsigned j = 0; j < muons.size(); j++)
        {
            if (!muonIsTight[j])
                continue;

            // Did not match leg 1
            if (muons[j]->pt < MUON_LEG1_PT)
                continue;
            if (!trigger->passObj(firedDoubleMuTriggers[i], 1, muons[j]->hltMatchBits))
                continue;

            // Look for a second lepton passing leg 2
            for (unsigned k = 0; k < muons.size(); k++)
            {
                if (j == k)
                    continue;
                if (!muonIsTight[k])
                    continue;

                if (muons[k]->pt < MUON_LEG2_PT)
                    continue;
                if (trigger->passObj(firedDoubleMuTriggers[i], 2, muons[k]->hltMatchBits))
                    matched = kTRUE;
            }

        }
        if (matched)
            matchedDoubleMuTriggers.push_back(firedDoubleMuTriggers[i]);
    }
    evtDoubleMuTriggered = matchedDoubleMuTriggers.size() > 0;



    //
    //  SINGLE MUON
    //

    for (unsigned i = 0; i < firedSingleMuTriggers.size(); i++)
    {
        bool matched = kFALSE;
        for (unsigned j = 0; j < muons.size(); j++)
        {
            if (!muonIsTight[j])
                continue;

            if (muons[j]->pt < MUON_SINGLE_PT)
                continue;
            if (trigger->passObj(firedSingleMuTriggers[i], 1, muons[j]->hltMatchBits))
                matched = kTRUE;
        }
        if (matched)
            matchedSingleMuTriggers.push_back(firedSingleMuTriggers[i]);
    }
    evtSingleMuTriggered = matchedSingleMuTriggers.size() > 0;

    evtMuonTriggered = evtSingleMuTriggered || evtDoubleMuTriggered;



    //
    //  DOUBLE ELECTRON
    //

    for (unsigned i = 0; i < firedDoubleElTriggers.size(); i++)
    {
        bool matched = kFALSE;
        for (unsigned j = 0; j < electrons.size(); j++)
        {
            if (!electronIsTight[j])
                continue;

            // Did not match leg 1
            if (electrons[j]->pt < ELEC_LEG1_PT)
                continue;
            if (!trigger->passObj(firedDoubleElTriggers[i], 1, electrons[j]->hltMatchBits))
                continue;

            // Look for a second lepton passing leg 2
            for (unsigned k = 0; k < electrons.size(); k++)
            {
                if (j == k)
                    continue;
                if (!electronIsTight[k])
                    continue;

                if (electrons[k]->pt < ELEC_LEG2_PT)
                    continue;
                if (trigger->passObj(firedDoubleElTriggers[i], 2, electrons[k]->hltMatchBits))
                    matched = kTRUE;
            }

        }
        if (matched)
            matchedDoubleElTriggers.push_back(firedDoubleElTriggers[i]);
    }
    evtDoubleElTriggered = matchedDoubleElTriggers.size() > 0;



    //
    //  SINGLE ELECTRON
    //

    for (unsigned i = 0; i < firedSingleElTriggers.size(); i++)
    {
        bool matched = kFALSE;
        for (unsigned j = 0; j < electrons.size(); j++)
        {
            if (!electronIsTight[j])
                continue;

            if (electrons[j]->pt < ELEC_SINGLE_PT)
                continue;
            if (trigger->passObj(firedSingleElTriggers[i], 1, electrons[j]->hltMatchBits))
                matched = kTRUE;
        }
        if (matched)
            matchedSingleElTriggers.push_back(firedSingleElTriggers[i]);
    }
    evtSingleElTriggered = matchedSingleElTriggers.size() > 0;

    evtElectronTriggered = evtSingleElTriggered || evtDoubleElTriggered;



    //
    //  DATA HANDLING
    //

    if (isData)
    {
        if      (evtDoubleMuTriggered)
            passTrigger = sampleName.Contains("DoubleMuon");
        else if (evtSingleMuTriggered)
            passTrigger = sampleName.Contains("SingleMuon");
        else if (evtDoubleElTriggered)
            passTrigger = sampleName.Contains("DoubleEG") || sampleName.Contains("EGamma");
        else if (evtSingleElTriggered)
            passTrigger = sampleName.Contains("SingleElectron") || sampleName.Contains("EGamma");
        else
            passTrigger = kFALSE;
    }
    else
        passTrigger = evtMuonTriggered || evtElectronTriggered;

    if (!passTrigger && !isSignal && !isDrellYan)
        return kTRUE;
    hTotalEvents->Fill(5);



    //
    //  FILL OBJECT INFO
    //


    // Muons
    for (unsigned j = 0; j < muons.size(); j++)
    {
        TMuon* muon = muons[j];

        bool firedLeg1 = kFALSE, firedLeg2 = kFALSE, firedSingle = kFALSE;

        for (unsigned i = 0; i < matchedDoubleMuTriggers.size(); i++)
        {
            if (trigger->passObj(matchedDoubleMuTriggers[i], 1, muon->hltMatchBits)
                    && (muon->pt > MUON_LEG1_PT))
                firedLeg1 = kTRUE;
            if (trigger->passObj(matchedDoubleMuTriggers[i], 2, muon->hltMatchBits)
                    && (muon->pt > MUON_LEG2_PT))
                firedLeg2 = kTRUE;
        }
        for (unsigned i = 0; i < matchedSingleMuTriggers.size(); i++)
        {
            if (trigger->passObj(matchedSingleMuTriggers[i], 1, muon->hltMatchBits)
                    && (muon->pt > MUON_SINGLE_PT))
                firedSingle = kTRUE;
        }

        muonFiredLeg1.push_back(firedLeg1);
        muonFiredLeg2.push_back(firedLeg2);
        muonFiredSingle.push_back(firedSingle);
    }


    // Electrons
    for (unsigned j = 0; j < electrons.size(); j++)
    {
        TElectron* electron = electrons[j];

        bool firedLeg1 = kFALSE, firedLeg2 = kFALSE, firedSingle = kFALSE;

        for (unsigned i = 0; i < matchedDoubleElTriggers.size(); i++)
        {
            if (trigger->passObj(matchedDoubleElTriggers[i], 1, electron->hltMatchBits)
                    && (electron->pt > ELEC_LEG1_PT))
                firedLeg1 = kTRUE;
            if (trigger->passObj(matchedDoubleElTriggers[i], 2, electron->hltMatchBits)
                    && (electron->pt > ELEC_LEG2_PT))
                firedLeg2 = kTRUE;
        }
        for (unsigned i = 0; i < matchedSingleElTriggers.size(); i++)
        {
            if (trigger->passObj(matchedSingleElTriggers[i], 1, electron->hltMatchBits)
                    && (electron->pt > ELEC_SINGLE_PT))
                firedSingle = kTRUE;
        }

        electronFiredLeg1.push_back(firedLeg1);
        electronFiredLeg2.push_back(firedLeg2);
        electronFiredSingle.push_back(firedSingle);
    }





    ////
    ////
    ////    SELECTION
    ////
    ////


    if      (params->selection == "loose")
    {
        if (nLooseMuons < 2 && nLooseElectrons < 2)
            return kTRUE;
    }
    else if (params->selection == "tight")
    {
        if (nTightMuons < 2 && nTightElectrons < 2)
            return kTRUE;
    }
    else
    {
        cout << "Invalid selection criterion" << endl;
        return kTRUE;
    }
    hTotalEvents->Fill(6);


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
