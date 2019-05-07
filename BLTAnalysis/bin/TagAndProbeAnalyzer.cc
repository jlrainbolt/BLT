#include "TagAndProbeAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 


TagAndProbeAnalyzer::TagAndProbeAnalyzer() : BLTSelector()
{

}

TagAndProbeAnalyzer::~TagAndProbeAnalyzer()
{

}

void TagAndProbeAnalyzer::Begin(TTree *tree)
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

    const bool storeGenInfo = params->datasetgroup == "zjets_m-50";

    // Particle selector, cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Weight utilities
    weights.reset(new WeightUtils(params->period, params->selection, false));

    // RNG
    rng.SetSeed(0);

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


    // Trigger
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_25ns";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if      (params->period == "2016")
    {
        muonTriggerNames.push_back("HLT_IsoMu24_v*");
        muonTriggerNames.push_back("HLT_IsoTkMu24_v*");
        electronTriggerNames.push_back("HLT_Ele27_WPTight_Gsf_v*");
    }
    else if (params->period == "2017")
    {
        muonTriggerNames.push_back("HLT_IsoMu27_v*");
        electronTriggerNames.push_back("HLT_Ele35_WPTight_Gsf_v*");
    }
    else if (params->period == "2018")
    {
        muonTriggerNames.push_back("HLT_IsoMu24_v*");
        electronTriggerNames.push_back("HLT_Ele32_WPTight_Gsf_v*");
    }



    //
    //  OUTPUT TREE
    //

    string outFileName = params->get_output_filename(params->selection);
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
    outTree->Branch(    "nPV",                      &nPV);
    outTree->Branch(    "nMuons",                   &nMuons);
    outTree->Branch(    "nElectrons",               &nElectrons);
    outTree->Branch(    "nLeptons",                 &nLeptons);


    outTree->Branch(    "tagP4",                    &tagP4);
    outTree->Branch(    "tagQ",                     &tagQ);
    outTree->Branch(    "tagPDG",                   &tagPDG);
    outTree->Branch(    "tagIsolation",             &tagIsolation);
    outTree->Branch(    "tagEnergySF",              &tagEnergySF);
    outTree->Branch(    "tagIsLoose",               &tagIsLoose);
    outTree->Branch(    "tagIsTight",               &tagIsTight);
    outTree->Branch(    "tagFiredSingle",           &tagFiredSingle);

    if      (params->selection == "mumu")
    {
        outTree->Branch(    "tagIsIsolated",            &tagIsIsolated);
        outTree->Branch(    "tagIsPF",                  &tagIsPF);
        outTree->Branch(    "tagIsTrackerHighPt",       &tagIsTrackerHighPt);
    }
    else if (params->selection == "ee")
    {
        outTree->Branch(    "tagScEta",                 &tagScEta);
        outTree->Branch(    "tagIsGap",                 &tagIsGap);
        outTree->Branch(    "tagIsV2Iso",               &tagIsV2Iso);
    }


    outTree->Branch(    "probeP4",                  &probeP4);
    outTree->Branch(    "probeQ",                   &probeQ);
    outTree->Branch(    "probePDG",                 &probeQ);
    outTree->Branch(    "probeIsolation",           &probeIsolation);
    outTree->Branch(    "probeEnergySF",            &probeEnergySF);
    outTree->Branch(    "probeIsLoose",             &probeIsLoose);
    outTree->Branch(    "probeIsTight",             &probeIsTight);
    outTree->Branch(    "probeFiredSingle",         &probeFiredSingle);

    if      (params->selection == "mumu")
    {
        outTree->Branch(    "probeIsIsolated",          &probeIsIsolated);
        outTree->Branch(    "probeIsPF",                &probeIsPF);
        outTree->Branch(    "probeIsTrackerHighPt",     &probeIsTrackerHighPt);
    }
    else if (params->selection == "ee")
    {
        outTree->Branch(    "probeScEta",               &probeScEta);
        outTree->Branch(    "probeIsGap",               &probeIsGap);
        outTree->Branch(    "probeIsV2Iso",             &probeIsV2Iso);
    }


    if (storeGenInfo)
    {
        outTree->Branch(    "nGenMuons",                &nGenMuons);
        outTree->Branch(    "nGenElectrons",            &nGenElectrons);
        outTree->Branch(    "nGenLeptons",              &nGenLeptons);

        outTree->Branch(    "genTagP4",                 &genTagP4);
        outTree->Branch(    "genTagQ",                  &genTagQ);
        outTree->Branch(    "genTagPDG",                &genTagPDG);
        outTree->Branch(    "genTagDeltaR",             &genTagDeltaR);
        outTree->Branch(    "genTagStatus",             &genTagStatus);
        outTree->Branch(    "genTagMotherPDG",          &genTagMotherPDG);
        outTree->Branch(    "nGenTagPhotons",           &nGenTagPhotons);
        outTree->Branch(    "genTagPhotonsP4",          &genTagPhotonsP4,       32000,      1);

        outTree->Branch(    "genProbeP4",               &genProbeP4);
        outTree->Branch(    "genProbeQ",                &genProbeQ);
        outTree->Branch(    "genProbePDG",              &genProbePDG);
        outTree->Branch(    "genProbeDeltaR",           &genProbeDeltaR);
        outTree->Branch(    "genProbeStatus",           &genProbeStatus);
        outTree->Branch(    "genProbeMotherPDG",        &genProbeMotherPDG);
        outTree->Branch(    "nGenProbePhotons",         &nGenTagPhotons);
        outTree->Branch(    "genProbePhotonsP4",        &genProbePhotonsP4,     32000,      1);

        outTree->Branch(    "genZP4",                   &genZP4);
    }



    //
    //  HISTOGRAMS
    //

    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);

    string outDirName = params->get_output_treename("hists");
    outFile->mkdir(outDirName.c_str());
    outFile->cd(outDirName.c_str());

    const unsigned  N_PT = (params->selection == "mumu") ? N_PT_MM : N_PT_EE;
    const unsigned  N_ETA = (params->selection == "mumu") ? N_ETA_MM : N_ETA_EE;

    float   *binsPt = (params->selection == "mumu") ? mumuPt : eePt;
    float   *binsEta = (params->selection == "mumu") ? mumuEta : eeEta;

    passed2d = new TH2F("passed_2d", "Passed", N_ETA, binsEta, N_PT, binsPt);
    failed2d = new TH2F("failed_2d", "Failed", N_ETA, binsEta, N_PT, binsPt);

    passedAll = new TH1F("passed_all", "Passed", N_MLL, MLL_MIN, MLL_MAX);
    failedAll = new TH1F("failed_all", "Failed", N_MLL, MLL_MIN, MLL_MAX);

    if (storeGenInfo)
    {
        genPassedAll = new TH1F("gen_passed_all", "Gen Passed", N_MLL, MLL_MIN, MLL_MAX);
        genFailedAll = new TH1F("gen_failed_all", "Gen Failed", N_MLL, MLL_MIN, MLL_MAX);
    }

    for (unsigned i = 0; i <= N_PT; i++)
    {
        TString ptString;
        
        if (i < N_PT)
            ptString.Form("(%g,%g)", binsPt[i], binsPt[i+1]);
        else
            ptString.Form("(%g+)", binsPt[i]);

        vector<TH1F*> passedVec, failedVec, genPassedVec, genFailedVec;

        for (unsigned j = 0; j <= N_ETA; j++)
        {
            TString etaString;

            if (j < N_ETA)
                etaString.Form("(%g,%g)", binsEta[j], binsEta[j+1]);
            else
                etaString.Form("(%g,%g)", binsEta[0], binsEta[j]);

            TString name = "_pt" + ptString + "_eta" + etaString;
            TString title = ": pT " + ptString + ", eta " + etaString;

            passedVec.push_back(new TH1F("passed" + name, "Passed" + title, N_MLL, MLL_MIN, MLL_MAX));
            failedVec.push_back(new TH1F("failed" + name, "Failed" + title, N_MLL, MLL_MIN, MLL_MAX));

            if (storeGenInfo)
            {
                genPassedVec.push_back(new TH1F("gen_passed" + name, "Gen Passed" + title, N_MLL, MLL_MIN, MLL_MAX));
                genFailedVec.push_back(new TH1F("gen_failed" + name, "Gen Failed" + title, N_MLL, MLL_MIN, MLL_MAX));
            }
        }
        passedBin.push_back(passedVec);
        failedBin.push_back(failedVec);

        if (storeGenInfo)
        {
            genPassedBin.push_back(genPassedVec);
            genFailedBin.push_back(genFailedVec);
        }
    }


    ReportPostBegin();
}


Bool_t TagAndProbeAnalyzer::Process(Long64_t entry)
{


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

    const bool storeGenInfo = params->datasetgroup == "zjets_m-50";



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
        float nPU   = fInfo->nPUmean;
        PUWeight    = weights->GetPUWeight(nPU);

        // ECAL weight
        ECALWeight  = fInfo->ecalWeight;
    }
    else
    {
        // Lumi mask
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;

        genWeight   = 1;
        PUWeight    = 1;
        ECALWeight  = 1;
    }
    hTotalEvents->Fill(2);
    weight = genWeight * PUWeight * ECALWeight;



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

        // Kinematic requirements
        if (muon->pt < MUON_PT_MIN)
            continue;

        if (fabs(muon->eta) > MUON_ETA_MAX)
            continue;

        muons.push_back(muon);
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);
    nMuons = muons.size();


    // Electrons
    vector<TElectron*> electrons;
    for (int i = 0; i < fElectronArr->GetEntries(); i++)
    {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        // Kinematic requirements
        if (electron->pt < ELEC_PT_MIN)
            continue;

        if (fabs(electron->eta) > ELEC_ETA_MAX)
            continue;

        electrons.push_back(electron);
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);
    nElectrons = electrons.size();

    nLeptons = nMuons + nElectrons;





    ////
    ////
    ////    PRESELECTION
    ////
    ////

    // Require exactly two leptons of given flavor

    if      (params->selection == "mumu")
    {
        if (nMuons != 2)
            return kTRUE;
        if (nElectrons != 0)
            return kTRUE;
    }
    else if (params->selection == "electron")
    {
        if (nElectrons != 2)
            return kTRUE;
        if (nMuons != 0)
            return kTRUE;
    }
    else
        return kTRUE;

    hTotalEvents->Fill(5);





    ////
    ////
    ////    SELECTION
    ////
    ////


    TLorentzVector totalP4;



    //
    //  MUONS
    //

    if (params->selection == "mumu")
    {
        // Charge requirement
        if (muons[0]->q == muons[1]->q)
            return kTRUE;

        TLorentzVector  muon1P4,  muon2P4;


        // Energy corrections
        for (unsigned i = 0; i < 2; i++)
            muons[i]->pt *= particleSelector->GetRochesterCorrection(muons[i]);

        copy_p4(muons[0], MUON_MASS, muon1P4);
        copy_p4(muons[1], MUON_MASS, muon2P4);

        totalP4 = muon1P4 + muon2P4;


        // Mass requirement
        if ((totalP4.M() < MLL_MIN) || (totalP4.M() > MLL_MAX))
            return kTRUE;


        // Randomly pick muon 1 or 2 to try as tag
        unsigned T = round(rng.Rndm());

        if (muons[T]->pt < MUON_PT_THRESH)  // if pt is too low, switch to the other one
            T = 1 - T;

        TMuon *tag = muons[T];
        
        // Check trigger
        tagFiredSingle = kFALSE;
        for (unsigned i = 0; i < muonTriggerNames.size(); i++)
        {
            if (trigger->passObj(muonTriggerNames[i], 1, tag->hltMatchBits))
                tagFiredSingle = kTRUE;
        }
        if (!tagFiredSingle)
            return kTRUE;

        // Check tight ID
        tagIsTight = particleSelector->PassMuonID(tag, cuts->tightHZZMuonID);
        if (!tagIsTight)
            return kTRUE;


        // Fill containers
        TMuon *probe = muons[1-T];

        copy_p4(tag, MUON_MASS, tagP4);
        tagQ                    = tag->q;
        tagPDG                  = -13 * tagQ;
        tagIsolation            = particleSelector->GetMuonIso(tag);
        tagEnergySF             = particleSelector->GetRochesterCorrection(tag);
        tagIsLoose              = particleSelector->PassMuonID(tag, cuts->looseHZZMuonID);
        tagIsIsolated           = particleSelector->PassMuonIso(tag, cuts->wpHZZMuonIso);
        tagIsPF                 = tag->typeBits & baconhep::kPFMuon;
        tagIsTrackerHighPt      = particleSelector->PassMuonID(tag, cuts->trackerHighPtMuonID);

        copy_p4(probe, MUON_MASS, probeP4);
        probeQ                  = probe->q;
        probePDG                = -13 * probeQ;
        probeIsolation          = particleSelector->GetMuonIso(probe);
        probeEnergySF           = particleSelector->GetRochesterCorrection(probe);
        probeIsLoose            = particleSelector->PassMuonID(probe, cuts->looseHZZMuonID);
        probeIsTight            = particleSelector->PassMuonID(probe, cuts->tightHZZMuonID);
        probeIsIsolated         = particleSelector->PassMuonIso(probe, cuts->wpHZZMuonIso);
        probeIsPF               = probe->typeBits & baconhep::kPFMuon;
        probeIsTrackerHighPt    = particleSelector->PassMuonID(probe, cuts->trackerHighPtMuonID);
        
        probeFiredSingle = kFALSE;
        for (unsigned i = 0; i < muonTriggerNames.size(); i++)
        {
            if (trigger->passObj(muonTriggerNames[i], 1, probe->hltMatchBits))
                probeFiredSingle = kTRUE;
        }
    }



    if (storeGenInfo)
    {
        hTotalEvents->Fill(6);

        genTagP4.Clear();       genProbeP4.Clear();     genZP4.Clear();
        genTagQ = 0;            genProbeQ = 0;          genTagPDG = 0;          genProbePDG = 0;
        genTagDeltaR = 0;       genProbeDeltaR = 0;     genTagStatus = 0;       genProbeStatus = 0;
        genTagMotherPDG = 0;    genProbeMotherPDG = 0;  nGenTagPhotons = 0;     nGenProbePhotons = 0;
        genTagPhotonsP4->Delete();                      genProbePhotonsP4->Delete();



        //
        //  GEN PARTICLES
        //

        vector<TGenParticle*> genMuons, genElectrons, genPhotons;

        for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if ((particle->pdgId == 23) && (particle->status == 62))    // this is the Z
                copy_p4(particle, genZP4);

            if (particle->status != 1)          // use final-state leptons
                continue;

            if      (abs(particle->pdgId) == 13)
                genMuons.push_back(particle);

            else if (abs(particle->pdgId) == 11)
                genElectrons.push_back(particle);

            else if (abs(particle->pdgId) == 22)
                genPhotons.push_back(particle);
        }
        sort(genMuons.begin(), genMuons.end(), sort_by_higher_pt<TGenParticle>);
        sort(genElectrons.begin(), genElectrons.end(), sort_by_higher_pt<TGenParticle>);

        nGenMuons       = genMuons.size();
        nGenElectrons   = genElectrons.size();
        nGenLeptons     = nGenMuons + nGenElectrons;



        //
        //  MATCHING
        //

        TGenParticle *genTag, *genProbe;
        vector<TGenParticle*> genFlavor;

        if      (params->selection == "mumu")
            genFlavor = genMuons;
        else if (params->selection == "ee")
            genFlavor = genElectrons;

        unsigned nGenFlavor = genFlavor.size();

        if (nGenFlavor == 0)   // don't try to match if there are no gen leptons
            goto fill;


        // Find DeltaR for all gen leptons (of selected flavor)
        vector<float> tagDeltaR(nGenFlavor), probeDeltaR(nGenFlavor);

        for (unsigned i = 0; i < nGenFlavor; i++)
        {
            TLorentzVector genP4;
            copy_p4(genFlavor[i], genP4);
            tagDeltaR[i] = tagP4.DeltaR(genP4);
            probeDeltaR[i] = probeP4.DeltaR(genP4);
        }

        // Find indices of mimimum DeltaR for each lepton
        unsigned T = min_element(tagDeltaR.begin(), tagDeltaR.end()) - tagDeltaR.begin();
        unsigned P = min_element(probeDeltaR.begin(), probeDeltaR.end()) - probeDeltaR.begin();


        // Fill tag info
        genTag = genFlavor[T];
        copy_p4(genTag, genTagP4);
        genTagDeltaR    = tagP4.DeltaR(genTagP4);
        genTagPDG       = genTag->pdgId;
        genTagQ         = genTagPDG < 0 ? 1 : -1;
        genTagStatus    = genTag->status;

        TGenParticle *tagMother = genTag;
        do
            tagMother = (TGenParticle*) fGenParticleArr->At(tagMother->parent);
        while ((tagMother->pdgId == genTagPDG) && (tagMother->parent >= 0));
        genTagMotherPDG = tagMother->pdgId;


        // Fill probe info
        genProbe = genFlavor[P];
        copy_p4(genProbe, genProbeP4);
        genProbeDeltaR  = probeP4.DeltaR(genProbeP4);
        genProbePDG     = genProbe->pdgId;
        genProbeQ       = genProbePDG < 0 ? 1 : -1;
        genProbeStatus  = genProbe->status;

        TGenParticle *probeMother = genProbe;
        do
            probeMother = (TGenParticle*) fGenParticleArr->At(probeMother->parent);
        while ((probeMother->pdgId == genProbePDG) && (probeMother->parent >= 0));
        genProbeMotherPDG = probeMother->pdgId;


        // Find FSR photons
        for (unsigned i = 0; i < genPhotons.size(); i++)
        {
            if (genPhotons[i]->parent == genTag->parent)
            {
                TLorentzVector *photonP4 = (TLorentzVector*) genTagPhotonsP4->ConstructedAt(nGenTagPhotons);
                copy_p4(genPhotons[i], photonP4);
                nGenTagPhotons++;
            }

            if (genPhotons[i]->parent == genProbe->parent)
            {
                TLorentzVector *photonP4 = (TLorentzVector*) genProbePhotonsP4->ConstructedAt(nGenProbePhotons);
                copy_p4(genPhotons[i], photonP4);
                nGenProbePhotons++;
            }
        }
    }



fill:

    //
    //  FILL
    //

    const unsigned  N_ETA = (params->selection == "mumu") ? N_ETA_MM : N_ETA_EE;

    unsigned I = passed2d->GetYaxis()->FindBin(probeP4.Pt()) - 1;
    unsigned J = passed2d->GetXaxis()->FindBin(probeP4.Eta()) - 1;

    if (probeIsTight)
    {
        passed2d->Fill(probeP4.Eta(), probeP4.Pt(), weight);
        passedAll->Fill(totalP4.M(), weight);
        passedBin[I][J]->Fill(totalP4.M(), weight);
        passedBin[I][N_ETA]->Fill(totalP4.M(), weight);

        if (storeGenInfo)
        {
            genPassedAll->Fill(genZP4.M(), weight);
            genPassedBin[I][J]->Fill(genZP4.M(), weight);
            genPassedBin[I][N_ETA]->Fill(genZP4.M(), weight);
        }
    }
    else
    {
        failed2d->Fill(probeP4.Eta(), probeP4.Pt(), weight);
        failedAll->Fill(totalP4.M(), weight);
        failedBin[I][J]->Fill(totalP4.M(), weight);
        failedBin[I][N_ETA]->Fill(totalP4.M(), weight);

        if (storeGenInfo)
        {
            genFailedAll->Fill(genZP4.M(), weight);
            genFailedBin[I][J]->Fill(genZP4.M(), weight);
            genFailedBin[I][N_ETA]->Fill(genZP4.M(), weight);
        }
    }


    if (isData)
        hTotalEvents->Fill(6);
    else
        hTotalEvents->Fill(7);


    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void TagAndProbeAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void TagAndProbeAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void TagAndProbeAnalyzer::ReportPostTerminate()
{
    std::cout << "  ==== Terminate Job =========================================" << std::endl;
    std::cout << "  output   : " << params->get_output_filename(params->selection) << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;
    std::cout << "  ============================================================" << std::endl;
}


// _____________________________________________________________________________
// Main function

int main(int argc, char **argv)
{
    std::unique_ptr<TagAndProbeAnalyzer> selector(new TagAndProbeAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
