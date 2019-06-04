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

    const bool isSignal     = params->datasetgroup == "zz_4l";

    // Particle selector, cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Weight utilities
    weights.reset(new WeightUtils(params->period, params->selection, false));

    // Lumi mask
    lumiMask = RunLumiRangeMap();
    string jsonFileName;
    if (params->period == "2012")
        jsonFileName = "Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
    lumiMask.AddJSONFile(cmssw_base + data_dir + jsonFileName);


    // Trigger
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_v2";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if (params->period == "2012")
    {
        doubleMuonTriggers.push_back("HLT_Mu17_Mu8_v*");
        doubleMuonTriggers.push_back("HLT_Mu17_TkMu8_v*");
        singleMuonTriggers.push_back("HLT_IsoMu24_eta2p1_v*");

        doubleElecTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
        singleElecTriggers.push_back("HLT_Ele27_WP80_v*");
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
    if (isSignal)
        outTree->Branch(    "sampleName",               &sampleName);
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
    outTree->Branch(    "PUWeight",                 &PUWeight);
    outTree->Branch(    "PUWeightUp",               &PUWeightUp);
    outTree->Branch(    "PUWeightDown",             &PUWeightDown);
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
    outTree->Branch(    "muonUncorrP4",             &muonUncorrP4_,         32000,      1);
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
//  outTree->Branch(    "muonTrigEffLeg1Data",      &muonTrigEffLeg1Data);
//  outTree->Branch(    "muonTrigEffLeg1MC",        &muonTrigEffLeg1MC);
//  outTree->Branch(    "muonTrigEffLeg2Data",      &muonTrigEffLeg2Data);
//  outTree->Branch(    "muonTrigEffLeg2MC",        &muonTrigEffLeg2MC);

    outTree->Branch(    "electronP4",               &electronP4_,           32000,      1);
    outTree->Branch(    "electronUncorrP4",         &electronUncorrP4_,     32000,      1);
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
    outTree->Branch(    "electronIsMVA",            &electronIsMVA);
    outTree->Branch(    "electronIsIsolated",       &electronIsIsolated);
    outTree->Branch(    "electronIsGap",            &electronIsGap);
    outTree->Branch(    "electronFiredLeg1",        &electronFiredLeg1);
    outTree->Branch(    "electronFiredLeg2",        &electronFiredLeg2);
    outTree->Branch(    "electronFiredSingle",      &electronFiredSingle);
//  outTree->Branch(    "electronTrigEffLeg1Data",  &electronTrigEffLeg1Data);
//  outTree->Branch(    "electronTrigEffLeg1MC",    &electronTrigEffLeg1MC);
//  outTree->Branch(    "electronTrigEffLeg2Data",  &electronTrigEffLeg2Data);
//  outTree->Branch(    "electronTrigEffLeg2MC",    &electronTrigEffLeg2MC);

    if (isSignal)
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

    sampleName.Clear();
    evtMuonTriggered = kFALSE;      evtDoubleMuTriggered = kFALSE;  evtSingleMuTriggered = kFALSE;
    evtElectronTriggered = kFALSE;  evtDoubleElTriggered = kFALSE;  evtSingleElTriggered = kFALSE;

    nLooseMuons     = 0;            nLooseElectrons = 0;            nLooseLeptons   = 0;
    nTightMuons     = 0;            nTightElectrons = 0;            nTightLeptons   = 0; 
    genWeight       = 1;            ECALWeight      = 1;            nPU             = 0;
    PUWeight        = 1;            PUWeightUp      = 1;            PUWeightDown    = 1;
    hasTauDecay     = kFALSE;

    muonP4_->Delete();              muonUncorrP4_->Delete();        muonCharge.clear();
    muonEnergySF.clear();           muonEnergySFUp.clear();         muonEnergySFDown.clear();
    muonIDSF.clear();               muonIsolation.clear();
    muonIsTight.clear();            muonIsLoose.clear();            muonIsIsolated.clear();
    muonIsPF.clear();               muonIsTrackerHighPt.clear();
    muonFiredLeg1.clear();          muonFiredLeg2.clear();          muonFiredSingle.clear();
//  muonTrigEffLeg1Data.clear();
//  muonTrigEffLeg1MC.clear();      muonTrigEffLeg2Data.clear();    muonTrigEffLeg2MC.clear(); 
                                                                                                        
    electronP4_->Delete();          electronUncorrP4_->Delete();    electronCharge.clear();
    electronEnergySF.clear();       electronEnergySFUp.clear();     electronEnergySFDown.clear();
    electronIDSF.clear();           electronRecoSF.clear();         electronIsolation.clear();      
    electronScEta.clear();          electronScEt.clear();           electronIsGap.clear();
    electronIsIsolated.clear();
    electronIsTight.clear();        electronIsLoose.clear();        electronIsMVA.clear();
    electronFiredLeg1.clear();      electronFiredLeg2.clear();      electronFiredSingle.clear();
//  electronTrigEffLeg1Data.clear();
//  electronTrigEffLeg1MC.clear();  electronTrigEffLeg2Data.clear();electronTrigEffLeg2MC.clear();

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

    weights->SetSampleName(params->dataset);
    sampleName = params->dataset;
 
    const bool isSignal     = params->datasetgroup == "zz_4l";
    const bool isDrellYan   = params->datasetgroup == "zjets_m-50";



    //
    //  EVENT INFO
    //

    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;
    nPV         = fPVArr->GetEntries();

    if (!isData)
    {
        // MC sample weight for ZZtoXXXX samples
        if (isSignal)
            genWeight = weights->GetSampleWeight();


        // Pileup weight
        nPU = fInfo->nPUmean;
        float nPUUp = fInfo->nPUmean + fabs((float) fInfo->nPUp - (float) fInfo->nPU);
        float nPUDown = fInfo->nPUmean - fabs((float) fInfo->nPUm - (float) fInfo->nPU);

        PUWeight = weights->GetPUWeight(nPU);
        PUWeightUp = weights->GetPUWeight(nPUUp);
        PUWeightDown = weights->GetPUWeight(nPUDown);
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

            if (!isSignal)                  // don't waste time if we aren't storing the leptons
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

    for (unsigned i = 0; i < doubleMuonTriggers.size(); i++)            // double muon
    {
        if (trigger->pass(doubleMuonTriggers[i], fInfo->triggerBits))
            evtDoubleMuTriggered = kTRUE;
    }
    for (unsigned i = 0; i < singleMuonTriggers.size(); i++)            // single muon
    {
        if (trigger->pass(singleMuonTriggers[i], fInfo->triggerBits))
            evtSingleMuTriggered = kTRUE;
    }
    for (unsigned i = 0; i < doubleElecTriggers.size(); i++)            // double electron
    {
        if (trigger->pass(doubleElecTriggers[i], fInfo->triggerBits))
            evtDoubleElTriggered = kTRUE;
    }
    for (unsigned i = 0; i < singleElecTriggers.size(); i++)            // single electron
    {
        if (trigger->pass(singleElecTriggers[i], fInfo->triggerBits))
            evtSingleElTriggered = kTRUE;
    }

    evtMuonTriggered = evtSingleMuTriggered || evtDoubleMuTriggered;
    evtElectronTriggered = evtSingleElTriggered || evtDoubleElTriggered;

    bool passTrigger;
    if (isData)
    {
        if      (evtDoubleMuTriggered)
            passTrigger = sampleName.Contains("DoubleMu");
        else if (evtSingleMuTriggered)
            passTrigger = sampleName.Contains("SingleMu");
        else if (evtDoubleElTriggered)
            passTrigger = sampleName.Contains("DoubleElectron");
        else if (evtSingleElTriggered)
            passTrigger = sampleName.Contains("SingleElectron");
        else
            passTrigger = kFALSE;
    }
    else
        passTrigger = evtMuonTriggered || evtElectronTriggered;

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
        TLorentzVector *uncorrP4 = (TLorentzVector*) muonUncorrP4_->ConstructedAt(i);
        TLorentzVector *corrP4 = (TLorentzVector*) muonP4_->ConstructedAt(i);

        // Uncorr momentum
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
        muonIsTrackerHighPt.push_back(kFALSE);

        if (muonIsLoose.back())
            nLooseMuons++;
        if (muonIsTight.back())
            nTightMuons++;

        // Trigger bools and SFs
        bool firedLeg1 = kFALSE, firedLeg2 = kFALSE, firedSingle = kFALSE;

        for (unsigned i = 0; i < doubleMuonTriggers.size(); i++)
        {
            if (trigger->passObj(doubleMuonTriggers[i], 1, muon->hltMatchBits))
                firedLeg1 = kTRUE;
            if (trigger->passObj(doubleMuonTriggers[i], 2, muon->hltMatchBits))
                firedLeg2 = kTRUE;
        }
        for (unsigned i = 0; i < singleMuonTriggers.size(); i++)
        {
            if (trigger->passObj(singleMuonTriggers[i], 1, muon->hltMatchBits))
                firedSingle = kTRUE;
        }

        muonFiredLeg1.push_back(firedLeg1);
        muonFiredLeg2.push_back(firedLeg2);
        muonFiredSingle.push_back(firedSingle);
/*
        pair<float, float> trigEff;
        trigEff = weights->GetDoubleMuonTriggerEff(muon, 1);
        muonTrigEffLeg1Data.push_back(trigEff.first);
        muonTrigEffLeg1MC.push_back(trigEff.second);

        trigEff = weights->GetDoubleMuonTriggerEff(muon, 2);
        muonTrigEffLeg2Data.push_back(trigEff.first);
        muonTrigEffLeg2MC.push_back(trigEff.second);
*/
    }



    //
    //  ELECTRONS
    //

    for (unsigned i = 0; i < electrons.size(); i++)
    {
        TElectron* electron = electrons[i];

        // Initialize entries in TClonesArrays
        TLorentzVector *uncorrP4 = (TLorentzVector*) electronUncorrP4_->ConstructedAt(i);
        TLorentzVector *corrP4 = (TLorentzVector*) electronP4_->ConstructedAt(i);

        // Uncorr momentum
        copy_p4(electrons[i], ELE_MASS, uncorrP4);
        electronScEta.push_back(electron->scEta);
        electronScEt.push_back(electron->scEt);

        // Charge
        electronCharge.push_back(electron->q);

        // Energy correction
        electronEnergySFUp.push_back(particleSelector->GetElectronCorrection(electron, "up"));
        electronEnergySFDown.push_back(particleSelector->GetElectronCorrection(electron, "down"));
        electronEnergySF.push_back(particleSelector->GetElectronCorrection(electron));
        electron->pt = electron->ptHZZ4l;
        copy_p4(electrons[i], ELE_MASS, corrP4);

        // ID scale factor
        electronIDSF.push_back(weights->GetHZZElectronIDSF(electron));
        electronRecoSF.push_back(1);

        // Isolation
        electronIsolation.push_back(particleSelector->GetElectronIso(electron));

        // ID/iso bools
        electronIsTight.push_back(particleSelector->PassElectronID(electron, cuts->tightHZZElectronID));
        electronIsLoose.push_back(particleSelector->PassElectronID(electron, cuts->looseHZZElectronID));
        electronIsMVA.push_back(particleSelector->PassElectronMVA(electron, cuts->wpHZZElectronMVA));
        electronIsIsolated.push_back(particleSelector->PassElectronIso(electron, cuts->wpHZZElectronIso));
        electronIsGap.push_back(electron->fiducialBits & kIsGap);

        if (electronIsLoose.back())
            nLooseElectrons++;
        if (electronIsTight.back())
            nTightElectrons++;

        // Trigger bools and SFs
        bool firedLeg1 = kFALSE, firedLeg2 = kFALSE, firedSingle = kFALSE;

        for (unsigned i = 0; i < doubleElecTriggers.size(); i++)
        {
            if (trigger->passObj(doubleElecTriggers[i], 1, electron->hltMatchBits))
                firedLeg1 = kTRUE;
            if (trigger->passObj(doubleElecTriggers[i], 2, electron->hltMatchBits))
                firedLeg2 = kTRUE;
        }
        for (unsigned i = 0; i < singleElecTriggers.size(); i++)
        {
            if (trigger->passObj(singleElecTriggers[i], 1, electron->hltMatchBits))
                firedSingle = kTRUE;
        }

        electronFiredLeg1.push_back(firedLeg1);
        electronFiredLeg2.push_back(firedLeg2);
        electronFiredSingle.push_back(firedSingle);
/*
        pair<float, float> trigEff;
        trigEff = weights->GetDoubleElectronTriggerEff(electron, 1);
        electronTrigEffLeg1Data.push_back(trigEff.first);
        electronTrigEffLeg1MC.push_back(trigEff.second);

        trigEff = weights->GetDoubleElectronTriggerEff(electron, 2);
        electronTrigEffLeg2Data.push_back(trigEff.first);
        electronTrigEffLeg2MC.push_back(trigEff.second);
*/
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
