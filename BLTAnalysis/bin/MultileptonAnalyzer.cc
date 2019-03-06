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

    const bool isSignal = params->datasetgroup == "zz_4l";

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
    lumiMask.AddJSONFile(cmssw_base + data_dir + jsonFileName);


    // Trigger
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if      (params->period == "2016")
    {
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*");
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*");
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");
        muonTriggerNames.push_back("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*");
        muonTriggerNames.push_back("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*");

        electronTriggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*");
    }
    else if (params->period == "2017")
    {
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*");
        muonTriggerNames.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*");

        electronTriggerNames.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*");
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
    outTree->Branch(    "evtElectronTriggered",     &evtElectronTriggered);
    outTree->Branch(    "genWeight",                &genWeight);
    outTree->Branch(    "ECALWeight",               &ECALWeight);
    outTree->Branch(    "PUWeight",                 &PUWeight);
    outTree->Branch(    "nPU",                      &nPU);
    outTree->Branch(    "nPV",                      &nPV);

    outTree->Branch(    "nLooseMuons",              &nLooseMuons);
    outTree->Branch(    "nLooseElectrons",          &nLooseElectrons);
    outTree->Branch(    "nLooseLeptons",            &nLooseLeptons);
    outTree->Branch(    "nTightMuons",              &nTightMuons);
    outTree->Branch(    "nTightElectrons",          &nTightElectrons);
    outTree->Branch(    "nTightLeptons",            &nTightLeptons);

    outTree->Branch(    "muonP4",                   &muonP4_,                       32000,      1);
    outTree->Branch(    "muonUncorrectedP4",        &muonUncorrectedP4_,            32000,      1);
    outTree->Branch(    "muonQ",                    &muonCharge);
    outTree->Branch(    "muonEnergySF",             &muonEnergySF);
    outTree->Branch(    "muonIDSF",                 &muonIDSF);
    outTree->Branch(    "muonIsolation",            &muonIsolation);
    outTree->Branch(    "muonIsTight",              &muonIsTight);
    outTree->Branch(    "muonIsLoose",              &muonIsLoose);
    outTree->Branch(    "muonIsIsolated",           &muonIsIsolated);
    outTree->Branch(    "muonIsPF",                 &muonIsPF);
    outTree->Branch(    "muonIsTrackerHighPt",      &muonIsTrackerHighPt);
    outTree->Branch(    "muonFiredLeg1",            &muonFiredLeg1);
    outTree->Branch(    "muonFiredLeg2",            &muonFiredLeg2);
    outTree->Branch(    "muonTrigEffLeg1Data",      &muonTrigEffLeg1Data);
    outTree->Branch(    "muonTrigEffLeg1MC",        &muonTrigEffLeg1MC);
    outTree->Branch(    "muonTrigEffLeg2Data",      &muonTrigEffLeg2Data);
    outTree->Branch(    "muonTrigEffLeg2MC",        &muonTrigEffLeg2MC);

    outTree->Branch(    "electronP4",               &electronP4_,                   32000,      1);
    outTree->Branch(    "electronUncorrectedP4",    &electronUncorrectedP4_,        32000,      1);
    outTree->Branch(    "electronQ",                &electronCharge);
    outTree->Branch(    "electronEnergySF",         &electronEnergySF);
    outTree->Branch(    "electronIDSF",             &electronIDSF);
    outTree->Branch(    "electronRecoSF",           &electronRecoSF);
    outTree->Branch(    "electronIsolation",        &electronIsolation);
    outTree->Branch(    "electronScEta",            &electronScEta);
    outTree->Branch(    "electronIsTight",          &electronIsTight);
    outTree->Branch(    "electronIsLoose",          &electronIsLoose);
    outTree->Branch(    "electronIsV2Iso",          &electronIsV2Iso);
    outTree->Branch(    "electronIsGap",            &electronIsGap);
    outTree->Branch(    "electronFiredLeg1",        &electronFiredLeg1);
    outTree->Branch(    "electronFiredLeg2",        &electronFiredLeg2);
    outTree->Branch(    "electronTrigEffLeg1Data",  &electronTrigEffLeg1Data);
    outTree->Branch(    "electronTrigEffLeg1MC",    &electronTrigEffLeg1MC);
    outTree->Branch(    "electronTrigEffLeg2Data",  &electronTrigEffLeg2Data);
    outTree->Branch(    "electronTrigEffLeg2MC",    &electronTrigEffLeg2MC);

    if (isSignal)
    {
        outTree->Branch(    "nFinalStateMuons",         &nFinalStateMuons);
        outTree->Branch(    "nFinalStateElectrons",     &nFinalStateElectrons);
        outTree->Branch(    "nFinalStateLeptons",       &nFinalStateLeptons);
        outTree->Branch(    "nHardProcMuons",           &nHardProcMuons);
        outTree->Branch(    "nHardProcElectrons",       &nHardProcElectrons);
        outTree->Branch(    "nHardProcLeptons",         &nHardProcLeptons);

        outTree->Branch(    "finalStateLeptonsP4",      &finalStateLeptonsP4);
        outTree->Branch(    "finalStateMuonP4",         &finalStateMuonP4_,         32000,      1);
        outTree->Branch(    "finalStateMuonQ",          &finalStateMuonQ);
        outTree->Branch(    "finalStateMuonZIndex",     &finalStateMuonZIndex);
        outTree->Branch(    "finalStateElectronP4",     &finalStateElectronP4_,     32000,      1);
        outTree->Branch(    "finalStateElectronQ",      &finalStateElectronQ);
        outTree->Branch(    "finalStateElectronZIndex", &finalStateElectronZIndex);

        outTree->Branch(    "hardProcLeptonsP4",        &hardProcLeptonsP4);
        outTree->Branch(    "hardProcMuonP4",           &hardProcMuonP4_,           32000,      1);
        outTree->Branch(    "hardProcMuonQ",            &hardProcMuonQ);
        outTree->Branch(    "hardProcMuonZIndex",       &hardProcMuonZIndex);
        outTree->Branch(    "hardProcElectronP4",       &hardProcElectronP4_,       32000,      1);
        outTree->Branch(    "hardProcElectronQ",        &hardProcElectronQ);
        outTree->Branch(    "hardProcElectronZIndex",   &hardProcElectronZIndex);
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

    nLooseMuons = 0;                nLooseElectrons = 0;            nLooseLeptons = 0;
    nTightMuons = 0;                nTightElectrons = 0;            nTightLeptons = 0; 

    muonP4_->Delete();              muonUncorrectedP4_->Delete();   muonCharge.clear();
    muonEnergySF.clear();           muonIDSF.clear();               muonIsolation.clear();
    muonIsTight.clear();            muonIsLoose.clear();            muonIsIsolated.clear();
    muonIsPF.clear();               muonIsTrackerHighPt.clear();
    muonFiredLeg1.clear();          muonFiredLeg2.clear();          muonTrigEffLeg1Data.clear();
    muonTrigEffLeg1MC.clear();      muonTrigEffLeg2Data.clear();    muonTrigEffLeg2MC.clear(); 
                                                                                                        
    electronP4_->Delete();          electronUncorrectedP4_->Delete();   electronCharge.clear();
    electronEnergySF.clear();       electronIDSF.clear();               electronRecoSF.clear();
    electronIsolation.clear();      electronScEta.clear();          electronIsGap.clear();
    electronIsTight.clear();        electronIsLoose.clear();        electronIsV2Iso.clear();
    electronFiredLeg1.clear();      electronFiredLeg2.clear();      electronTrigEffLeg1Data.clear();
    electronTrigEffLeg1MC.clear();  electronTrigEffLeg2Data.clear();electronTrigEffLeg2MC.clear();

    nFinalStateMuons = 0;           nFinalStateElectrons = 0;       nFinalStateLeptons = 0; 
    nHardProcMuons = 0;             nHardProcElectrons = 0;         nHardProcLeptons = 0; 
    finalStateLeptonsP4.Clear();    hardProcLeptonsP4.Clear();
    finalStateMuonP4_->Delete();    finalStateMuonQ.clear();        finalStateMuonZIndex.clear();
    finalStateElectronP4_->Delete();finalStateElectronQ.clear();    finalStateElectronZIndex.clear();
    hardProcMuonP4_->Delete();      hardProcMuonQ.clear();          hardProcMuonZIndex.clear();
    hardProcElectronP4_->Delete();  hardProcElectronQ.clear();      hardProcElectronZIndex.clear();






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
 
    const bool isSignal = params->datasetgroup == "zz_4l";



    //
    //  EVENT INFO
    //

    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;
    genWeight   = 1;
    ECALWeight  = 1;
    PUWeight    = 1;
    nPU         = 0;
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
        ECALWeight = fInfo->ecalWeight;
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

    if (isSignal)
    {
        finalStateLeptonsP4 = TLorentzVector();
        hardProcLeptonsP4 = TLorentzVector();

        // Particle loop
        for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if  ((abs(particle->pdgId) == 13 || abs(particle->pdgId) == 11)  &&  particle->parent >= 0)
            {
                TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(particle->parent);

                // Hard process
                if (mother->pdgId == 23)
                {
                    int charge = -1 * copysign(1, particle->pdgId);

                    if      (abs(particle->pdgId) == 13)
                    {
                        TLorentzVector *p4 = (TLorentzVector*)hardProcMuonP4_->ConstructedAt(nHardProcMuons);
                        copy_p4(particle, MUON_MASS, p4);
                        hardProcLeptonsP4 = hardProcLeptonsP4 + *p4;
                        hardProcMuonQ.push_back(charge);
                        hardProcMuonZIndex.push_back(particle->parent);
                        nHardProcMuons++;
                    }
                    else if (abs(particle->pdgId) == 11)
                    {
                        TLorentzVector *p4 = (TLorentzVector*)hardProcElectronP4_->ConstructedAt(nHardProcElectrons);
                        copy_p4(particle, ELE_MASS, p4);
                        hardProcLeptonsP4 = hardProcLeptonsP4 + *p4;
                        hardProcElectronQ.push_back(charge);
                        hardProcElectronZIndex.push_back(particle->parent);
                        nHardProcElectrons++;
                    }
                }

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
                            TLorentzVector *p4 = (TLorentzVector*)finalStateMuonP4_->ConstructedAt(nFinalStateMuons);
                            copy_p4(particle, MUON_MASS, p4);
                            finalStateLeptonsP4 = finalStateLeptonsP4 + *p4;
                            finalStateMuonQ.push_back(charge);
                            finalStateMuonZIndex.push_back(motherIndex);
                            nFinalStateMuons++;
                        }
                        else if (abs(particle->pdgId) == 11)
                        {
                            TLorentzVector *p4 = (TLorentzVector*)finalStateElectronP4_->ConstructedAt(nFinalStateElectrons);
                            copy_p4(particle, ELE_MASS, p4);
                            finalStateLeptonsP4 = finalStateLeptonsP4 + *p4;
                            finalStateElectronQ.push_back(charge);
                            finalStateElectronZIndex.push_back(motherIndex);
                            nFinalStateElectrons++;
                        }
                    }
                }
            }
        }

        nFinalStateLeptons      = nFinalStateMuons + nFinalStateElectrons;
        nHardProcLeptons        = nHardProcMuons + nHardProcElectrons;
    }



    //
    //  TRIGGER SELECTION
    //

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
        TLorentzVector *uncorrP4 = (TLorentzVector*)muonUncorrectedP4_->ConstructedAt(i);
        TLorentzVector *corrP4 = (TLorentzVector*)muonP4_->ConstructedAt(i);

        // Uncorrected momentum
        copy_p4(muon, MUON_MASS, uncorrP4);

        // Charge
        muonCharge.push_back(muon->q);

        // Rochester correction
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

        pair<float, float> trigEff;
        trigEff = weights->GetDoubleMuonTriggerEff(muon, 1);
        muonTrigEffLeg1Data.push_back(trigEff.first);
        muonTrigEffLeg1MC.push_back(trigEff.second);

        trigEff = weights->GetDoubleMuonTriggerEff(muon, 2);
        muonTrigEffLeg2Data.push_back(trigEff.first);
        muonTrigEffLeg2MC.push_back(trigEff.second);
    }



    //
    //  ELECTRONS
    //

    for (unsigned i = 0; i < electrons.size(); i++)
    {
        TElectron* electron = electrons[i];

        // Initialize entries in TClonesArrays
        TLorentzVector *uncorrP4 = (TLorentzVector*)electronUncorrectedP4_->ConstructedAt(i);
        TLorentzVector *corrP4 = (TLorentzVector*)electronP4_->ConstructedAt(i);

        // Uncorrected momentum
        copy_p4(electrons[i], ELE_MASS, uncorrP4);
        electronScEta.push_back(electron->scEta);

        // Charge
        electronCharge.push_back(electron->q);

        // Energy correction
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

        // Trigger bools and SFs
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

        pair<float, float> trigEff;
        trigEff = weights->GetDoubleElectronTriggerEff(electron, 1);
        electronTrigEffLeg1Data.push_back(trigEff.first);
        electronTrigEffLeg1MC.push_back(trigEff.second);

        trigEff = weights->GetDoubleElectronTriggerEff(electron, 2);
        electronTrigEffLeg2Data.push_back(trigEff.first);
        electronTrigEffLeg2MC.push_back(trigEff.second);
    }

    nLooseLeptons = nLooseMuons + nLooseElectrons;
    nTightLeptons = nTightMuons + nTightElectrons;






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
    hTotalEvents->Fill(5);


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
