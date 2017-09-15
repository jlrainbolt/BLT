#include "MultileptonAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 
bool sort_gen_pt (TGenParticle *i, TGenParticle *j) { return (i->pt > j->pt); }

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

    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);

    // Set the cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Trigger bits mapping file
    const std::string cmssw_base = getenv("CMSSW_BASE");
    std::string trigfilename = cmssw_base + "/src/BaconAna/DataFormats/data/HLTFile_v2";
    trigger.reset(new baconhep::TTrigger(trigfilename));

    if (params->selection == "mumu" || params->selection == "emu" || params->selection == "4mu" || params->selection == "2e2mu" || params->selection == "2e2muSign") {
        triggerNames.push_back("HLT_IsoMu24_eta2p1_v*");
    } else if (params->selection == "ee") {
        triggerNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
    }
//  triggerNames.push_back("HLT_PFJet320_v*");

    // Weight utility class
    if (params->selection == "4mu")
        weights.reset(new WeightUtils(params->period, "mumu", false));
    else if (params->selection == "2e2mu" || params->selection == "2e2muSign")
        weights.reset(new WeightUtils(params->period, "emu", false));
    else
        weights.reset(new WeightUtils(params->period, params->selection, false));

    // Lumi mask
    // Set up object to handle good run-lumi filtering if necessary
    lumiMask = RunLumiRangeMap();
    if (true) { // this will need to be turned off for MC
        string jsonFileName = cmssw_base + "/src/BLT/BLTAnalysis/test/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
        lumiMask.AddJSONFile(jsonFileName);
    }

    // muon momentum corrections
    muonCorr = new rochcor2012();

    // Prepare the output tree
    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");

    // event data
    outTree->Branch("runNumber", &runNumber);
    outTree->Branch("evtNumber", &evtNumber, "eventNumber/l");
    outTree->Branch("lumiSection", &lumiSection);
    outTree->Branch("triggerStatus", &triggerStatus);
    outTree->Branch("eventWeight", &eventWeight);
    outTree->Branch("nPU", &nPU);
    outTree->Branch("nPartons", &nPartons);

    outTree->Branch("met", &met);
    outTree->Branch("metPhi", &metPhi);

    // leptons
    outTree->Branch("leptonOneP4", &leptonOneP4);
    outTree->Branch("leptonOneIso", &leptonOneIso);
    outTree->Branch("leptonOneQ", &leptonOneQ);
    outTree->Branch("leptonOneFlavor", &leptonOneFlavor);
    outTree->Branch("leptonOneTrigger", &leptonOneTrigger);
    outTree->Branch("leptonTwoP4", &leptonTwoP4);
    outTree->Branch("leptonTwoIso", &leptonTwoIso);
    outTree->Branch("leptonTwoQ", &leptonTwoQ);
    outTree->Branch("leptonTwoFlavor", &leptonTwoFlavor);
    outTree->Branch("leptonTwoTrigger", &leptonTwoTrigger);
    outTree->Branch("leptonThreeP4", &leptonThreeP4);
    outTree->Branch("leptonThreeIso", &leptonThreeIso);
    outTree->Branch("leptonThreeQ", &leptonThreeQ);
    outTree->Branch("leptonThreeFlavor", &leptonThreeFlavor);
    outTree->Branch("leptonThreeTrigger", &leptonThreeTrigger);
    outTree->Branch("leptonFourP4", &leptonFourP4);
    outTree->Branch("leptonFourIso", &leptonFourIso);
    outTree->Branch("leptonFourQ", &leptonFourQ);
    outTree->Branch("leptonFourFlavor", &leptonFourFlavor);
    outTree->Branch("leptonFourTrigger", &leptonFourTrigger);

    outTree->Branch("GENleptonOneP4", &muon1);
    outTree->Branch("GENleptonTwoP4", &muon2);
    outTree->Branch("GENleptonThreeP4", &muon3);
    outTree->Branch("GENleptonFourP4", &muon4);

    // jets
    outTree->Branch("jetP4", &jetP4);
    outTree->Branch("jetD0", &jetD0);
    outTree->Branch("jetTag", &jetTag);
    //outTree->Branch("jetPUID", &jetPUID);
    outTree->Branch("jetFlavor", &jetFlavor);
    outTree->Branch("bjetP4", &bjetP4);
    outTree->Branch("bjetTag", &bjetTag);
    //outTree->Branch("bjetPUID", &bjetPUID);
    outTree->Branch("bjetFlavor", &bjetFlavor);
    outTree->Branch("bjetD0", &bjetD0);
    outTree->Branch("genBJetP4", &genBJetP4);
    outTree->Branch("genBJetTag", &genBJetTag);

    // object counters
    outTree->Branch("nMuons", &nMuons);
    outTree->Branch("nElectrons", &nElectrons);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nFwdJets", &nFwdJets);
    outTree->Branch("nBJets", &nBJets);

    // Event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(),"TotalEvents",10,0.5,10.5);

    // Histograms
    hZEta = new TH1D("z_eta", "z_eta", 100, -10, 10);
    hZPt = new TH1D("z_pt", "z_pt", 100, 0, 100);
    hZMass = new TH1D("z_mass", "z_mass", 100, 80, 100);
    hMu1Eta = new TH1D("mu1_eta", "mu1_eta", 100, -10, 10);
    hMu1Pt = new TH1D("mu1_pt", "mu1_pt", 100, 0, 100);
    hMu2Eta = new TH1D("mu2_eta", "mu2_eta", 100, -10, 10);
    hMu2Pt = new TH1D("mu2_pt", "mu2_pt", 100, 0, 100);


    ReportPostBegin();
}

Bool_t MultileptonAnalyzer::Process(Long64_t entry)
{
    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;
    hTotalEvents->Fill(1);

    //if (entry%1==0)  std::cout << "... Processing event: " << entry << "." << std::endl;
    if (entry%10000==0) { 
        std::cout << "... Processing event: " << entry << " Run: " << fInfo->runNum 
            << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum << "." << std::endl;
    }

    const bool isRealData = (fInfo->runNum != 1);
    particleSelector->SetRealData(isRealData);

    /////////////////////////
    // Gen-level selection //
    /////////////////////////
    
//    bool isProcess = false;
    bool isAccepted = false;

    if (!isRealData)
    {
        vector<TGenParticle*> leptons, zbosons;
//        float massZ = -1;

        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i)
        {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            if (abs(particle->pdgId) == 23)
//                massZ = particle->mass;
                zbosons.push_back(particle);
            else if (particle->status == 1 && abs(particle->pdgId) == 13)
                leptons.push_back(particle);
        }
        std::sort(leptons.begin(), leptons.end(), sort_gen_pt);
        if (zbosons.size() > 1)
        {
            float zmass_diff = 1000;
            for (unsigned int i = 0; i < zbosons.size(); i++)
            {
                if (abs(zbosons[i]->mass - 91.2) < zmass_diff)
                {
                    zmass_diff = abs(zbosons[i]->mass - 91.2);
                    the_z.SetPtEtaPhiM(zbosons[i]->pt, zbosons[i]->eta, zbosons[i]->phi, zbosons[i]->mass);
                }
            }
        }

        if (params->selection == "mumu")
        {
            if (leptons.size() >= 2)
            {
                muon1.SetPtEtaPhiM(leptons[0]->pt, leptons[0]->eta, leptons[0]->phi, leptons[0]->mass);
                for (unsigned i = 1; i < leptons.size(); ++i)
                {
                    if (leptons[0]->pdgId != leptons[i]->pdgId)
                    {
                        muon2.SetPtEtaPhiM(leptons[i]->pt, leptons[i]->eta, leptons[i]->phi, leptons[i]->mass);
                        break;
                    }
                }
                hZEta->Fill(the_z.Eta());
                hZPt->Fill(the_z.Pt());
                hZMass->Fill(the_z.M());
                hMu1Eta->Fill(muon1.Eta());
                hMu1Pt->Fill(muon1.Pt());
                hMu2Eta->Fill(muon2.Eta());
                hMu2Pt->Fill(muon2.Pt());
            }
        }
/*
        if (params->selection == "mumu")
        {
            if (leptons.size() > 1)
            {
                for (unsigned int j = 0; j < leptons.size(); j++)
                {
                    for (unsigned int i = 0; i < j; i++)
                    {
                        float Mll = sqrt(2 * leptons[i]->pt * leptons[j]->pt * (cosh(leptons[i]->eta - leptons[j]->eta) - cos(leptons[i]->phi - leptons[j]->phi)));
                        if (massZ > 0 && abs(Mll - massZ) <= 3 && leptons[i]->pdgId != leptons[j]->pdgId)
                        {
                            this->genEvents++;
                            isProcess = true;

                            if (leptons[i]->pt > 25 && leptons[j]->pt > 25
                                    && abs(leptons[i]->eta) < 2.1 && abs(leptons[j]->eta) < 2.1
                                    && Mll > 60 && Mll < 120)
                            {
                                this->genAcceptedEvents++;
                                isAccepted = true;
                            }
                        }

                        if (isProcess)
                            break;
                    }
                    if (isProcess)
                        break;
                }
            }
        }
*/        else if (params->selection == "4mu")
        {
            if (leptons.size() > 3)
            {
                muon1.SetPtEtaPhiM(leptons[0]->pt, leptons[0]->eta, leptons[0]->phi, leptons[0]->mass);
                muon2.SetPtEtaPhiM(leptons[1]->pt, leptons[1]->eta, leptons[1]->phi, leptons[1]->mass);
                muon3.SetPtEtaPhiM(leptons[2]->pt, leptons[2]->eta, leptons[2]->phi, leptons[2]->mass);
                muon4.SetPtEtaPhiM(leptons[3]->pt, leptons[3]->eta, leptons[3]->phi, leptons[3]->mass);
            }
        }
    }

    // Apply lumi mask
    if (isRealData) {
        RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(!lumiMask.HasRunLumi(rl)) 
            return kTRUE;
    }
    hTotalEvents->Fill(2);

    /* Trigger selection */
    bool passTrigger = false;
//    bool passTrigger2 = false;
    for (unsigned i = 0; i < triggerNames.size(); ++i) {
        passTrigger |= trigger->pass(triggerNames[i], fInfo->triggerBits);
    }
//        passTrigger2 |= trigger->pass("HLT_PFJet320_v*", fInfo->triggerBits);
    if (!passTrigger)// || !passTrigger2)
        return kTRUE;

    hTotalEvents->Fill(3);

    /////////////////////
    // Fill event info //
    /////////////////////

    eventWeight   = 1;
    runNumber     = fInfo->runNum;
    evtNumber     = fInfo->evtNum;
    lumiSection   = fInfo->lumiSec;
    triggerStatus = passTrigger;
    nPU           = fPVArr->GetEntries();
    if (!isRealData) {
        eventWeight *= weights->GetPUWeight(fInfo->nPUmean); // pileup reweighting
    }

    ///////////////////////
    // Generator objects //
    ///////////////////////

    if (!isRealData) {
        unsigned count = 0;
        for (int i = 0; i < fGenParticleArr->GetEntries(); ++i) {
            TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);

            //cout << particle->status << ", "
            //     << particle->pdgId  << ", "
            //     << particle->parent
            //     << endl;

            if (particle->status == 3 && (abs(particle->pdgId) < 6 || particle->pdgId == 21)) {
                ++count;
            }
        }
        nPartons = count-4; // This is saved for reweighting inclusive DY and combining it with parton binned DY
        //cout << nPartons << "\n" << endl;
    } else {
        nPartons = 0;
    }

    ////////////////////
    // Select objects //
    ////////////////////

    /* Vertices */
    if (fInfo->hasGoodPV) {
        assert(fPVArr->GetEntries() != 0);
        TVector3 pv;
        copy_xyz((TVertex*) fPVArr->At(0), pv);
        particleSelector->SetPV(pv);
    } else {
        return kTRUE;
    }
    hTotalEvents->Fill(4);
    particleSelector->SetNPV(fInfo->nPU + 1);
    particleSelector->SetRho(fInfo->rhoJet);

    /* MUONS */
    /* Apply a preselection so we can make a collection of muons to clean against */
    vector<TMuon*> tmp_muons;
    for (int i=0; i < fMuonArr->GetEntries(); i++) {
        TMuon* muon = (TMuon*) fMuonArr->At(i);
        assert(muon);

        if (
                muon->pt > 5 
                && fabs(muon->eta) < 2.4
                // tight muon ID
                //&& (muon->typeBits & baconhep::kPFMuon) 
                && (muon->typeBits & baconhep::kGlobal) 
                && muon->muNchi2    < 10.
                && muon->nMatchStn  > 1
                && muon->nPixHits   > 0
                && fabs(muon->d0)   < 0.2
                && fabs(muon->dz)   < 0.5
                && muon->nTkLayers  > 5 
                && muon->nValidHits > 0
           ) {
            tmp_muons.push_back(muon);
        }
    }
    sort(tmp_muons.begin(), tmp_muons.end(), sort_by_higher_pt<TMuon>);

    // Second pass
    //int trigger_muon_index = -1;
    vector<TLorentzVector> veto_muons;
    vector<TMuon*> muons;
    for (unsigned i = 0; i < tmp_muons.size(); i++) {
        TMuon* muon = tmp_muons[i];

        TLorentzVector muonP4;
        copy_p4(tmp_muons[i], MUON_MASS, muonP4);

        // Remove muon track pt from muon track isolation variable
        for (unsigned j = i+1; j < tmp_muons.size(); j++) {
            TLorentzVector muon_j;
            copy_p4(tmp_muons[j], MUON_MASS, muon_j);

            if (muonP4.DeltaR(muon_j) < 0.3) {
                muon->trkIso03 = max(0., muon->trkIso03 - muon_j.Pt());
                tmp_muons[j]->trkIso03 = max(0., tmp_muons[j]->trkIso03 - muonP4.Pt());
            }
        }

        // Apply rochester muon momentum corrections
        float qter = 1.;
        if (isRealData) {
            muonCorr->momcor_data(muonP4, muon->q, 0, qter);
        } else {
            muonCorr->momcor_mc(muonP4, muon->q, 0, qter);
        }

        // Fill containers
        if (muon->trkIso03/muonP4.Pt() < 0.1) {
            if (muonP4.Pt() > 20)
                veto_muons.push_back(muonP4);

            if (muonP4.Pt() > 5) {
                muon->pt  = muonP4.Pt();
                muon->eta = muonP4.Eta();
                muon->phi = muonP4.Phi();
                muons.push_back(muon);
                ++nMuons;
            }
        }
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);

    /* ELECTRONS */
    std::vector<TElectron*> electrons;
    for (int i=0; i<fElectronArr->GetEntries(); i++) {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        if (
                electron->pt > 5 //7 
                && fabs(electron->eta) < 2.5
                && particleSelector->PassElectronID(electron, cuts->tightElID)
                && particleSelector->PassElectronIso(electron, cuts->tightElIso, cuts->EAEl)
           ) {
            electrons.push_back(electron);
            ++nElectrons;

            // trigger matching
            //bool triggered = false;
            //for (unsigned i = 0; i < triggerNames.size(); ++i) {
            //    triggered |= trigger->passObj(triggerNames[i], 1, electron->hltMatchBits);
            //}
        }
    }

    std::sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

    /* JETS */
    TClonesArray* jetCollection;
    jetCollection = fAK5Arr;

    std::vector<TJet*> jets;
    std::vector<TJet*> fwdjets;
    std::vector<TJet*> bjets;
    std::vector<TJet*> genjets;
    nJets    = 0;
    nFwdJets = 0;
    nBJets   = 0;
    for (int i=0; i < jetCollection->GetEntries(); i++) {
        TJet* jet = (TJet*) jetCollection->At(i);
        assert(jet);

        // Prevent overlap of muons and jets
        TLorentzVector vJet; 
        vJet.SetPtEtaPhiM(jet->pt, jet->eta, jet->phi, jet->mass);
        bool muOverlap = false;
        for (const auto& mu: veto_muons) {
            if (vJet.DeltaR(mu) < 0.5) {
                muOverlap = true;
                break;
            }
        }
        bool elOverlap = false;
        //for (const auto& el: electrons) {
        //    if (vJet.DeltaR(el) < 0.5) {
        //        elOverlap = true;
        //        break;
        //    }
        //}

        if (!isRealData && abs(jet->mcFlavor) == 5) genjets.push_back(jet);

        if (
                jet->pt > 30 
                && fabs(jet->eta < 4.7)
                && particleSelector->PassJetID(jet, cuts->looseJetID)
           ) {

            if (fabs(jet->eta) <= 2.4) { 
                if (
                        particleSelector->PassJetPUID(jet, cuts->looseJetID)
                        && !muOverlap 
                        && !elOverlap
                   ) { 
                    if (isRealData) {
                        if (jet->csv > 0.898) {
                        //if (jet->bmva > 0.783) {
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            ++nJets;
                            jets.push_back(jet);
                        }
                    } else {
                        if (particleSelector->BTagModifier(jet, "CSVT")) { 
                            bjets.push_back(jet);
                            ++nBJets;
                        } else {
                            ++nJets;
                            jets.push_back(jet);
                        }
                    }
                }
            } else {
                if (jet->pt > 30) {
                    fwdjets.push_back(jet);
                    ++nFwdJets;
                }
            }
        }
    }
    std::sort(fwdjets.begin(), fwdjets.end(), sort_by_higher_pt<TJet>);
    std::sort(bjets.begin(), bjets.end(), sort_by_higher_pt<TJet>);

    // Add additional b jets to the central jet collection
    if (bjets.size() > 1) {
        jets.insert(jets.end(), bjets.begin()+1, bjets.end());
    }
    std::sort(jets.begin(), jets.end(), sort_by_higher_pt<TJet>);

    // MET //
    met    = fInfo->pfMETC;
    metPhi = fInfo->pfMETCphi;

    ///////////////////////////////
    // Apply analysis selections //
    ///////////////////////////////

    nMuons     = muons.size();
    nElectrons = electrons.size();

    if (params->selection == "mumu") {

        if (muons.size() < 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Find leading positive and negatively charged muons and convert to TLorentzVectors
        TLorentzVector muonOneP4, muonTwoP4;
        unsigned muonTwoIndex = 1;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
/*
        if (muons.size() == 2) {
            muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
        } else if (muons.size() > 2) {
            if (muons[0]->q != muons[1]->q) {
                muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
            } else if (muons[0]->q != muons[2]->q) {
                muonTwoIndex = 2;
                muonTwoP4.SetPtEtaPhiM(muons[2]->pt, muons[2]->eta, muons[2]->phi, MUON_MASS);
            }
        }
*/
        muonTwoP4.SetPtEtaPhiM(muons[muonTwoIndex]->pt, muons[muonTwoIndex]->eta, muons[muonTwoIndex]->phi, MUON_MASS);
        if (muons[0]->q * muons[1]->q > 0)
            return kTRUE;

        if (
            muonOneP4.Pt() < 25 || muonTwoP4.Pt() < 25
            || fabs(muonOneP4.Eta()) > 2.1 || fabs(muonTwoP4.Eta()) > 2.1
           )
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector dimuon;
        dimuon = muonOneP4 + muonTwoP4;
        if (dimuon.M() < 12)
            return kTRUE;
        hTotalEvents->Fill(7);

        leptonOneP4      = muonOneP4;
        leptonOneIso     = muons[0]->trkIso03;
        leptonOneQ       = muons[0]->q;
        leptonOneFlavor  = 13;

        leptonTwoP4      = muonTwoP4;
        leptonTwoIso     = muons[muonTwoIndex]->trkIso03;
        leptonTwoQ       = muons[muonTwoIndex]->q;
        leptonTwoFlavor  = 13;

        if (!passTrigger)
            return kTRUE;
        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            leptonOneTrigger |= trigger->passObj(triggerNames[i], 1, muons[0]->hltMatchBits);
            leptonTwoTrigger |= trigger->passObj(triggerNames[i], 1, muons[muonTwoIndex]->hltMatchBits);
        }
/*
        if (!isRealData) {
            eventWeight *= weights->GetMuonRecoEff(muonOneP4);
            eventWeight *= weights->GetMuonRecoEff(muonTwoP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonTwoP4);
            eventWeight *= (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
        }
*/
    } else if (params->selection == "4mu") {

        if (muons.size() < 4)
            return kTRUE;
//        hTotalEvents->Fill(5);

        // Check muon charges and convert to TLorentzVectors
        if (muons[0]->q * muons[1]->q > 0 || muons[2]->q * muons[3]->q > 0)
            return kTRUE;

        TLorentzVector muonOneP4, muonTwoP4, muonThreeP4, muonFourP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
        muonThreeP4.SetPtEtaPhiM(muons[2]->pt, muons[2]->eta, muons[2]->phi, MUON_MASS);
        muonFourP4.SetPtEtaPhiM(muons[3]->pt, muons[3]->eta, muons[3]->phi, MUON_MASS);

        if (
            muonOneP4.Pt() < 25 || muonTwoP4.Pt() < 25
            || fabs(muonOneP4.Eta()) > 2.1 || fabs(muonTwoP4.Eta()) > 2.1
           )
            return kTRUE;
        hTotalEvents->Fill(6);

        TLorentzVector dimuon;
        dimuon = muonOneP4 + muonTwoP4;
        if (dimuon.M() < 12)
            return kTRUE;
        hTotalEvents->Fill(7);

        leptonOneP4       = muonOneP4;
        leptonOneIso      = muons[0]->trkIso03;
        leptonOneQ        = muons[0]->q;
        leptonOneFlavor   = 13;

        leptonTwoP4       = muonTwoP4;
        leptonTwoIso      = muons[1]->trkIso03;
        leptonTwoQ        = muons[1]->q;
        leptonTwoFlavor   = 13;

        leptonThreeP4     = muonThreeP4;
        leptonThreeIso    = muons[2]->trkIso03;
        leptonThreeQ      = muons[2]->q;
        leptonThreeFlavor = 13;

        leptonFourP4      = muonFourP4;
        leptonFourIso     = muons[3]->trkIso03;
        leptonFourQ       = muons[3]->q;
        leptonFourFlavor  = 13;

        if (!passTrigger)
            return kTRUE;
        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            leptonOneTrigger |= trigger->passObj(triggerNames[i], 1, muons[0]->hltMatchBits);
            leptonTwoTrigger |= trigger->passObj(triggerNames[i], 1, muons[1]->hltMatchBits);
            leptonThreeTrigger |= trigger->passObj(triggerNames[i], 1, muons[2]->hltMatchBits);
            leptonFourTrigger |= trigger->passObj(triggerNames[i], 1, muons[3]->hltMatchBits);
        }
/*
        if (!isRealData) {
            eventWeight *= weights->GetMuonRecoEff(muonOneP4);
            eventWeight *= weights->GetMuonRecoEff(muonTwoP4);
            eventWeight *= weights->GetMuonRecoEff(muonThreeP4);
            eventWeight *= weights->GetMuonRecoEff(muonThreeP4);

            // trigger weight
            // Not adjusted for four leptons
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonTwoP4);
            eventWeight *= (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
        }
*/
    } else if (params->selection == "ee") {

        if (electrons.size() != 2)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Find leading positive and negatively charged electrons and convert to TLorentzVectors
        TLorentzVector electronOneP4, electronTwoP4;
        unsigned electronTwoIndex = 1;
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
        if (electrons.size() == 2) {
            electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);
        } else if (electrons.size() > 2) {
            if (electrons[0]->q != electrons[1]->q) {
                electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);
            } else if (electrons[0]->q != electrons[2]->q) {
                electronTwoIndex = 2;
                electronTwoP4.SetPtEtaPhiM(electrons[2]->pt, electrons[2]->eta, electrons[2]->phi, ELE_MASS);
            }
        }

        TLorentzVector dielectron;
        dielectron = electronOneP4 + electronTwoP4;
        if (dielectron.M() < 12. || dielectron.M() > 70.)
            return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4      = electronOneP4;
        leptonOneIso     = electrons[0]->trkIso03;
        leptonOneQ       = electrons[0]->q;
        leptonOneFlavor  = 11;

        leptonTwoP4      = electronTwoP4;
        leptonTwoIso     = electrons[electronTwoIndex]->trkIso03;
        leptonTwoQ       = electrons[electronTwoIndex]->q;
        leptonTwoFlavor  = 11;
    } else if (params->selection == "emu") {

        if (muons.size() != 1 || electrons.size() != 1)
            return kTRUE;
        hTotalEvents->Fill(5);

        // Convert leading leptons to TLorentzVectors
        TLorentzVector electronP4, muonP4;
        electronP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
        muonP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);

        TLorentzVector dilepton;
        dilepton = muonP4 + electronP4;
        if (dilepton.M() > 12 && dilepton.M() < 70)
            return kTRUE;
        hTotalEvents->Fill(6);

        leptonOneP4      = muonP4;
        leptonOneIso     = muons[0]->trkIso03;
        leptonOneQ       = muons[0]->q;
        leptonOneFlavor  = 13;

        leptonTwoP4      = electronP4;
        leptonTwoIso     = electrons[0]->trkIso03;
        leptonTwoQ       = electrons[0]->q;
        leptonTwoFlavor  = 11;

        if (!isRealData) {
            eventWeight *= weights->GetMuonRecoEff(muonP4);

            // trigger efficiency
            pair<float, float> trigEff = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonP4);
            eventWeight *= trigEff.first/trigEff.second;
        }
    } else if (params->selection == "2e2mu") {

        if (muons.size() < 2 || electrons.size() < 2)
            return kTRUE;

        // Check lepton charges and convert to TLorentzVectors
        if (muons[0]->q * muons[1]->q > 0 || electrons[0]->q * electrons[1]->q > 0)
            return kTRUE;

        TLorentzVector muonOneP4, muonTwoP4, electronOneP4, electronTwoP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);

        if (
            muonOneP4.Pt() < 25 || muonTwoP4.Pt() < 25
//            || electronOneP4.Pt() > muonTwoP4.Pt() || electronTwoP4.Pt() > muonTwoP4.Pt()
//            || electronOneP4.DeltaR(electronTwoP4) < 0.05
            )
            return kTRUE;

        TLorentzVector dimuon;//, dielectron;
        dimuon = muonOneP4 + muonTwoP4;
//        dielectron = electronOneP4 + electronTwoP4;
        if (dimuon.M() < 12)
            return kTRUE;

        leptonOneP4       = muonOneP4;
        leptonOneIso      = muons[0]->trkIso03;
        leptonOneQ        = muons[0]->q;
        leptonOneFlavor   = 13;

        leptonTwoP4       = muonTwoP4;
        leptonTwoIso      = muons[1]->trkIso03;
        leptonTwoQ        = muons[1]->q;
        leptonTwoFlavor   = 13;

        leptonThreeP4     = electronOneP4;
        leptonThreeIso    = electrons[0]->trkIso03;
        leptonThreeQ      = electrons[0]->q;
        leptonThreeFlavor = 11;

        leptonFourP4      = electronTwoP4;
        leptonFourIso     = electrons[1]->trkIso03;
        leptonFourQ       = electrons[1]->q;
        leptonFourFlavor  = 11;

        if (!passTrigger)
            return kTRUE;
/*        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            leptonOneTrigger |= trigger->passObj(triggerNames[i], 1, muons[0]->hltMatchBits);
            leptonTwoTrigger |= trigger->passObj(triggerNames[i], 1, muons[1]->hltMatchBits);
            leptonThreeTrigger |= trigger->passObj(triggerNames[i], 1, muons[2]->hltMatchBits);
            leptonFourTrigger |= trigger->passObj(triggerNames[i], 1, muons[3]->hltMatchBits);
        }

        if (!isRealData) {
            eventWeight *= weights->GetMuonRecoEff(muonOneP4);
            eventWeight *= weights->GetMuonRecoEff(muonTwoP4);
            eventWeight *= weights->GetMuonRecoEff(muonThreeP4);
            eventWeight *= weights->GetMuonRecoEff(muonThreeP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonTwoP4);
            eventWeight *= (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
        }*/
    } else if (params->selection == "2e2muSign") {

        if (muons.size() < 2 || electrons.size() < 2)
            return kTRUE;

        // Check lepton charges and convert to TLorentzVectors
        if (muons[0]->q * muons[1]->q < 0 || electrons[0]->q * electrons[1]->q < 0)
            return kTRUE;

        TLorentzVector muonOneP4, muonTwoP4, electronOneP4, electronTwoP4;
        muonOneP4.SetPtEtaPhiM(muons[0]->pt, muons[0]->eta, muons[0]->phi, MUON_MASS);
        muonTwoP4.SetPtEtaPhiM(muons[1]->pt, muons[1]->eta, muons[1]->phi, MUON_MASS);
        electronOneP4.SetPtEtaPhiM(electrons[0]->pt, electrons[0]->eta, electrons[0]->phi, ELE_MASS);
        electronTwoP4.SetPtEtaPhiM(electrons[1]->pt, electrons[1]->eta, electrons[1]->phi, ELE_MASS);

        if (
            muonOneP4.Pt() < 25 || muonTwoP4.Pt() < 25
//            || electronOneP4.Pt() > muonTwoP4.Pt() || electronTwoP4.Pt() > muonTwoP4.Pt()
//            || electronOneP4.DeltaR(electronTwoP4) < 0.05
            )
            return kTRUE;

        TLorentzVector dimuon;//, dielectron;
        dimuon = muonOneP4 + muonTwoP4;
//        dielectron = electronOneP4 + electronTwoP4;
        if (dimuon.M() < 12)
            return kTRUE;

        leptonOneP4       = muonOneP4;
        leptonOneIso      = muons[0]->trkIso03;
        leptonOneQ        = muons[0]->q;
        leptonOneFlavor   = 13;

        leptonTwoP4       = muonTwoP4;
        leptonTwoIso      = muons[1]->trkIso03;
        leptonTwoQ        = muons[1]->q;
        leptonTwoFlavor   = 13;

        leptonThreeP4     = electronOneP4;
        leptonThreeIso    = electrons[0]->trkIso03;
        leptonThreeQ      = electrons[0]->q;
        leptonThreeFlavor = 11;

        leptonFourP4      = electronTwoP4;
        leptonFourIso     = electrons[1]->trkIso03;
        leptonFourQ       = electrons[1]->q;
        leptonFourFlavor  = 11;

        if (!passTrigger)
            return kTRUE;
/*        for (unsigned i = 0; i < triggerNames.size(); ++i) {
            leptonOneTrigger |= trigger->passObj(triggerNames[i], 1, muons[0]->hltMatchBits);
            leptonTwoTrigger |= trigger->passObj(triggerNames[i], 1, muons[1]->hltMatchBits);
            leptonThreeTrigger |= trigger->passObj(triggerNames[i], 1, muons[2]->hltMatchBits);
            leptonFourTrigger |= trigger->passObj(triggerNames[i], 1, muons[3]->hltMatchBits);
        }

        if (!isRealData) {
            eventWeight *= weights->GetMuonRecoEff(muonOneP4);
            eventWeight *= weights->GetMuonRecoEff(muonTwoP4);
            eventWeight *= weights->GetMuonRecoEff(muonThreeP4);
            eventWeight *= weights->GetMuonRecoEff(muonThreeP4);

            // trigger weight
            pair<float, float> trigEff1 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonOneP4);
            pair<float, float> trigEff2 = weights->GetTriggerEffWeight("HLT_IsoMu24_eta2p1_v*", muonTwoP4);
            eventWeight *= (1 - (1 - trigEff1.first)*(1 - trigEff2.first))/(1 - (1 - trigEff1.second)*(1 - trigEff2.second));
        }*/
    }


    ///////////////////
    // Fill jet info //
    ///////////////////

    if (bjets.size() > 0) {
        bjetP4.SetPtEtaPhiM(bjets[0]->pt, bjets[0]->eta, bjets[0]->phi, bjets[0]->mass);
        bjetD0     = bjets[0]->d0;
        bjetTag    = bjets[0]->csv;
        //bjetPUID   = bjets[0]->mva;

        if (isRealData) {
            bjetFlavor = 0.;
        } else {
            bjetFlavor = bjets[0]->mcFlavor;
        }

    } else {
        bjetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        bjetD0     = 0.;
        bjetTag    = 0.;
        bjetPUID   = 0.;
        bjetFlavor = 0.;
    }

    if (fwdjets.size() > 0) {
        jetP4.SetPtEtaPhiM(fwdjets[0]->pt, fwdjets[0]->eta, fwdjets[0]->phi, fwdjets[0]->mass);
        jetD0     = fwdjets[0]->d0;
        jetTag    = 0.;
        //jetPUID   = fwdjets[0]->mva;

        if (isRealData) {
            jetFlavor = 0.;
        } else {
            jetFlavor = fwdjets[0]->mcFlavor;
        }
    } else if (jets.size() > 0) {
        jetP4.SetPtEtaPhiM(jets[0]->pt, jets[0]->eta, jets[0]->phi, jets[0]->mass);
        jetD0  = jets[0]->d0;
        jetTag = jets[0]->csv;

        //jetTag = jets[0]->csv;
        //jetPUID   = fwdjets[0]->mva;

        if (isRealData) {
            jetFlavor = 0.;
        } else {
            jetFlavor = jets[0]->mcFlavor;
        }
    } else {
        jetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        jetD0     = 0.;
        jetTag    = 0.;
        jetFlavor = 0.;
        jetPUID   = 0.;
    } 

    if (genjets.size() > 0) {
        genBJetP4.SetPtEtaPhiM(genjets[0]->genpt, genjets[0]->geneta, genjets[0]->genphi, genjets[0]->genm);
        genBJetTag = genjets[0]->csv;
    } else {
        genBJetP4.SetPtEtaPhiM(0., 0., 0., 0.);
        genBJetTag = 0;
    }

    outTree->Fill();
    this->passedEvents++;

    if (!isRealData & isAccepted)
        this->passIdEvents++;
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
    std::cout << "  output   : " << params->get_output_filename("demoFile") << std::endl;
    std::cout << "           : Processed " << this->fileCount << " files with " << this->unskimmedEventCount << " unskimmed events." << std::endl;
    std::cout << "           : Selected " << this->passedEvents << " / " << this->totalEvents << " events." << std::endl;

    if (fInfo->runNum == 1)
    {
        std::cout << "           : Accepted " << this->genAcceptedEvents << " / " << this->genEvents << " gen-level events." << std::endl;
        std::cout << "           : Found " << this->passIdEvents << " / " << this->genAcceptedEvents << " events passing ID requirements." << std::endl << std::endl;
        std::cout << "  Acceptance: " << (float)this->genAcceptedEvents / (float)this->genEvents << std::endl;
        std::cout << "  Efficiency: " << (float)this->passIdEvents / (float)this->genAcceptedEvents << std::endl;
    }

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
