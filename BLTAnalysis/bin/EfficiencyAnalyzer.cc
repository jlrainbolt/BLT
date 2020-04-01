#include "EfficiencyAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;


bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 


EfficiencyAnalyzer::EfficiencyAnalyzer() : BLTSelector()
{

}

EfficiencyAnalyzer::~EfficiencyAnalyzer()
{

}

void EfficiencyAnalyzer::Begin(TTree *tree)
{


    //
    //  SETUP
    //

    // Parse command line option
    std::string tmp_option = GetOption();
    std::vector<std::string> options;
    std::regex re_whitespace("(\\s+)");  // split by white space
    std::copy(std::sregex_token_iterator(tmp_option.begin(), tmp_option.end(), re_whitespace, -1),
              std::sregex_token_iterator(), std::back_inserter(options));

    const std::string cmssw_base = getenv("CMSSW_BASE");
    const std::string data_dir = "/src/BLT/BLTAnalysis/data/";


    // Set the parameters
    params.reset(new Parameters());
    params->setup(options);

    // Particle selector, cuts
    cuts.reset(new Cuts());
    particleSelector.reset(new ParticleSelector(*params, *cuts));

    // Weight utilities
    weights.reset(new WeightUtils(params->period, params->selection, false));



    //
    //  OUTPUT TREE
    //

    string outFileName = params->get_output_filename("output");
    string outTreeName = params->get_output_treename("tree");

    outFile = new TFile(outFileName.c_str(),"RECREATE");
    outFile->cd();
    outTree = new TTree(outTreeName.c_str(), "bltTree");


    // Branches
    outTree->Branch(    "runNumber",                &runNumber);
    outTree->Branch(    "evtNumber",                &evtNumber);
    outTree->Branch(    "lumiSection",              &lumiSection);
    outTree->Branch(    "isMatched",                &isMatched);
    outTree->Branch(    "decayChannel",             &decayChannel);

    outTree->Branch(    "genWeight",                &genWeight);
    outTree->Branch(    "PUWeight",                 &PUWeight);
    outTree->Branch(    "nPU",                      &nPU);
    outTree->Branch(    "nPV",                      &nPV);

    outTree->Branch(    "nLooseMuons",              &nLooseMuons);
    outTree->Branch(    "nLooseElectrons",          &nLooseElectrons);
    outTree->Branch(    "nLooseLeptons",            &nLooseLeptons);
    outTree->Branch(    "nTightMuons",              &nTightMuons);
    outTree->Branch(    "nTightElectrons",          &nTightElectrons);
    outTree->Branch(    "nTightLeptons",            &nTightLeptons);
    outTree->Branch(    "nDressedMuons",            &nDressedMuons);
    outTree->Branch(    "nDressedElectrons",        &nDressedElectrons);
    outTree->Branch(    "nDressedLeptons",          &nDressedLeptons);

    outTree->Branch(    "muonP4",                   &muonP4_,               32000,      1);
    outTree->Branch(    "muonQ",                    &muonQ);
    outTree->Branch(    "muonIsolation",            &muonIsolation);
    outTree->Branch(    "muonIsTight",              &muonIsTight);
    outTree->Branch(    "muonIsLoose",              &muonIsLoose);
    outTree->Branch(    "muonIsIsolated",           &muonIsIsolated);
    outTree->Branch(    "muonIsMatched",            &muonIsMatched);
    outTree->Branch(    "muonMatchIdx",             &muonMatchIdx);

    outTree->Branch(    "electronP4",               &electronP4_,           32000,      1);
    outTree->Branch(    "electronQ",                &electronQ);
    outTree->Branch(    "electronIsolation",        &electronIsolation);
    outTree->Branch(    "electronIsTight",          &electronIsTight);
    outTree->Branch(    "electronIsLoose",          &electronIsLoose);
    outTree->Branch(    "electronIsV2Iso",          &electronIsV2Iso);
    outTree->Branch(    "electronIsMatched",        &electronIsMatched);
    outTree->Branch(    "electronMatchIdx",         &electronMatchIdx);

    outTree->Branch(    "dressedMuonP4",            &dressedMuonP4_,        32000,      1);
    outTree->Branch(    "dressedMuonQ",             &dressedMuonQ);
    outTree->Branch(    "dressedMuonIsMatched",     &dressedMuonIsMatched);
    outTree->Branch(    "dressedMuonMatchIdx",      &dressedMuonMatchIdx);

    outTree->Branch(    "dressedElectronP4",        &dressedElectronP4_,    32000,      1);
    outTree->Branch(    "dressedElectronQ",         &dressedElectronQ);
    outTree->Branch(    "dressedElectronIsMatched", &dressedElectronIsMatched);
    outTree->Branch(    "dressedElectronMatchIdx",  &dressedElectronMatchIdx);



    //
    //  HISTOGRAMS
    //

    // Total event counter
    string outHistName = params->get_output_treename("TotalEvents");
    hTotalEvents = new TH1D(outHistName.c_str(), "TotalEvents", 10, 0.5, 10.5);

    // Acceptance counters
    outHistName = params->get_output_treename("SelectedEvents");
    hSelectedEvents = new TH1D(outHistName.c_str(), "SelectedEvents", 10, 0.5, 10.5);

    outHistName = params->get_output_treename("MatchedEvents");
    hMatchedEvents = new TH1D(outHistName.c_str(), "MatchedEvents", 10, 0.5, 10.5);


    ReportPostBegin();
}


Bool_t EfficiencyAnalyzer::Process(Long64_t entry)
{
    if (params->selection != "dressed")
        return kTRUE;



    //
    //  CLEAR CONTAINERS
    //
    
    nLooseMuons     = 0;            nLooseElectrons = 0;            nLooseLeptons   = 0;
    nTightMuons     = 0;            nTightElectrons = 0;            nTightLeptons   = 0;
    genWeight       = 1;            PUWeight        = 1;            nPU             = 0;
    isMatched       = kFALSE;

    muonP4_->Delete();              muonIsMatched.clear();
    muonQ.clear();                  muonMatchIdx.clear();           muonIsolation.clear();
    muonIsTight.clear();            muonIsLoose.clear();            muonIsIsolated.clear();

    electronP4_->Delete();          electronIsMatched.clear();
    electronQ.clear();              electronMatchIdx.clear();       electronIsolation.clear();
    electronIsTight.clear();        electronIsLoose.clear();        electronIsV2Iso.clear();

    nDressedMuons = 0;              nDressedElectrons = 0;          nDressedLeptons = 0;
    dressedMuonP4_->Delete();       dressedMuonQ.clear();           dressedMuonMatchIdx.clear();
    dressedElectronP4_->Delete();   dressedElectronQ.clear();       dressedElectronMatchIdx.clear();
    dressedMuonIsMatched.clear();   dressedElectronIsMatched.clear();



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

    TString dataSetGroup = params->datasetgroup,    sampleName = params->dataset;;
    const bool isSignal     = sampleName.EqualTo("ZZTo4L") || sampleName.EqualTo("ZZTo4L_aMC") || sampleName.EqualTo("ZZTo4L_M-1");
    const bool isDrellYan   = sampleName.Contains("DYJetsToLL_M-50");

    // Reject (accidental?) data events
    if (isData)
        return kTRUE;



    //
    //  EVENT INFO
    //

    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;
    nPV         = fPVArr->GetEntries();

    // Gen weight
    genWeight = fGenEvtInfo->weight > 0 ? 1 : -1; 
    if (genWeight < 0)
        hTotalEvents->Fill(10);

    // Pileup weight
    nPU = fInfo->nPUmean;
    PUWeight = weights->GetPUWeight(nPU);

    float weight = genWeight * PUWeight;

    hTotalEvents->Fill(2);





    ////
    ////
    ////    GEN OBJECTS
    ////
    ////


    vector<TLorentzVector> dressedMuonP4, dressedElectronP4;
    vector<TGenParticle*> genMuons, genElectrons;

    for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
    {
        TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
        int motherIndex = particle->parent;

        if (motherIndex < 0)            // seg faults are bad!
            continue;

        TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(motherIndex);

        // Look for taus from a Z
        if ((abs(particle->pdgId) == 15) && (mother->pdgId == 23))
            return kTRUE;

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
        if (params->selection == "dressed")
        {
            TLorentzVector lepP4;
            copy_p4(particle, lepP4);

            for (int j = 0; j < fGenParticleArr->GetEntries(); j++)
            {
                TGenParticle* gamma = (TGenParticle*) fGenParticleArr->At(j);

                if ((gamma->pdgId != 22) || (gamma->status != 1))
                    continue;

                TLorentzVector gammaP4;
                copy_p4(gamma, 0, gammaP4);

                if (lepP4.DeltaR(gammaP4) < DR_DRESS)
                    lepP4 = lepP4 + gammaP4;
            }

            particle->pt = lepP4.Pt();
            particle->eta = lepP4.Eta();
            particle->phi = lepP4.Phi();
            particle->mass = lepP4.M();
        }

        if      (abs(particle->pdgId) == 13)
        {
            if (particle->pt < MU_PT_MIN)
                continue;
            if (fabs(particle->eta) > MU_ETA_MAX)
                continue;

            genMuons.push_back(particle);
            nDressedMuons++;
        }
        else if (abs(particle->pdgId) == 11)
        {
            if (particle->pt < EL_PT_MIN)
                continue;
            if (fabs(particle->eta) > EL_ETA_MAX)
                continue;

            genElectrons.push_back(particle);
            nDressedElectrons++;
        }
    } // END particle loop

    nDressedLeptons = nDressedMuons + nDressedElectrons;


    // Require correct number of leptons

    if (isSignal && (nDressedLeptons != 4))
        return kTRUE;
    if (isDrellYan && (nDressedLeptons != 2))
        return kTRUE;
    hTotalEvents->Fill(2);


    // Fill containers
    sort(genMuons.begin(), genMuons.end(), sort_by_higher_pt<TGenParticle>);
    sort(genElectrons.begin(), genElectrons.end(), sort_by_higher_pt<TGenParticle>);

    TLorentzVector genLeptonsP4;
    for (unsigned i = 0; i < nDressedMuons; i++)
    {
        TLorentzVector *p4_ = (TLorentzVector*) dressedMuonP4_->ConstructedAt(i);
        copy_p4(genMuons[i], MUON_MASS, p4_);

        TLorentzVector p4;
        copy_p4(genMuons[i], MUON_MASS, p4);
        dressedMuonP4.push_back(p4);
        genLeptonsP4 = genLeptonsP4 + p4;

        int charge = -1 * copysign(1, genMuons[i]->pdgId);
        dressedMuonQ.push_back(charge);

        dressedMuonMatchIdx.push_back(-1);
        dressedMuonIsMatched.push_back(kFALSE);
    }

    for (unsigned i = 0; i < nDressedElectrons; i++)
    {
        TLorentzVector *p4_ = (TLorentzVector*) dressedElectronP4_->ConstructedAt(i);
        copy_p4(genElectrons[i], ELE_MASS, p4_);

        TLorentzVector p4;
        copy_p4(genElectrons[i], ELE_MASS, p4);
        dressedElectronP4.push_back(p4);
        genLeptonsP4 = genLeptonsP4 + p4;

        int charge = -1 * copysign(1, genElectrons[i]->pdgId);
        dressedElectronQ.push_back(charge);

        dressedElectronMatchIdx.push_back(-1);
        dressedElectronIsMatched.push_back(kFALSE);
    }



    //
    //  SELECTION
    //


    // Total mass

    if (genLeptonsP4.M() < M_MIN || genLeptonsP4.M() > M_MAX)
        return kTRUE;
    hTotalEvents->Fill(3);



    // Dilepton mass

    for (unsigned j = 1; j < nDressedMuons; j++)
    {
        for (unsigned i = 0; i < j; i++)
        {
            if (dressedMuonQ[i] != dressedMuonQ[j])
            {
                TLorentzVector dimuonP4 = dressedMuonP4[i] + dressedMuonP4[j];

                if (dimuonP4.M() < MLL_MIN)
                    return kTRUE;
            }
        }
    }
    for (unsigned j = 1; j < nDressedElectrons; j++)
    {
        for (unsigned i = 0; i < j; i++)
        {
            if (dressedElectronQ[i] != dressedElectronQ[j])
            {
                TLorentzVector dielecP4 = dressedElectronP4[i] + dressedElectronP4[j];

                if (dielecP4.M() < MLL_MIN)
                    return kTRUE;
            }
        }
    }
    hTotalEvents->Fill(4);



    // Categorize

    TLorentzVector dressedMuonsP4, dressedElecsP4;

    for (unsigned i = 0; i < nDressedMuons; i++)
        dressedMuonsP4 = dressedMuonsP4 + dressedMuonP4[i];
    for (unsigned i = 0; i < nDressedElectrons; i++)
        dressedElecsP4 = dressedElecsP4 + dressedElectronP4[i];

    unsigned C = 0;                             // Index

    if      (nDressedMuons == 2 && nDressedElectrons == 0)    // mumu = 3
        C = 3;

    else if (nDressedMuons == 0 && nDressedElectrons == 2)    // ee   = 4
        C = 4;

    else if (nDressedMuons == 4 && nDressedElectrons == 0)    // 4m   = 6
        C = 6;

    else if (nDressedMuons == 2 && nDressedElectrons == 2)    // 2m2e = 7
        C = 7;

    else if (nDressedMuons == 0 && nDressedElectrons == 4)    // 4e   = 9
        C = 9;

    else
        return kTRUE;

    unsigned D = (C < 5) ? 2 : 5;
    decayChannel = C;
    hTotalEvents->Fill(5);



    // "Fiducial" acceptance

    vector<TLorentzVector> sorted_leps = dressedMuonP4;
    sorted_leps.insert(sorted_leps.end(), dressedElectronP4.begin(), dressedElectronP4.end());
    sort(sorted_leps.begin(), sorted_leps.end(), P4SortCondition);

    bool isFiducial = kTRUE;

    if (sorted_leps[0].Pt() < PT1_MIN)
        isFiducial = kFALSE;
    if (sorted_leps[1].Pt() < PT2_MIN)
        isFiducial = kFALSE;

    if (!isFiducial)
        return kTRUE;

    hTotalEvents->Fill(6);
    hSelectedEvents->Fill(1, weight);
    hSelectedEvents->Fill(C, weight);
    hSelectedEvents->Fill(D, weight);





    ////
    ////
    ////    RECO OBJECTS
    ////
    ////


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

        // Store only loose muons passing cuts
        if (muon->pt < MU_PT_MIN)
            continue;
        if (fabs(muon->eta) > MU_ETA_MAX)
            continue;

        if (particleSelector->PassMuonID(muon, cuts->looseHZZMuonID))
            muons.push_back(muon);
    }
    sort(muons.begin(), muons.end(), sort_by_higher_pt<TMuon>);


    // Electrons
    vector<TElectron*> electrons;
    for (int i = 0; i < fElectronArr->GetEntries(); i++)
    {
        TElectron* electron = (TElectron*) fElectronArr->At(i);
        assert(electron);

        // Store only loose electrons passing cuts
        if (electron->pt < EL_PT_MIN)
            continue;
        if (fabs(electron->eta) > EL_ETA_MAX)
            continue;

        if (particleSelector->PassElectronID(electron, cuts->looseHZZElectronID))
            electrons.push_back(electron);
    }
    sort(electrons.begin(), electrons.end(), sort_by_higher_pt<TElectron>);

    vector<TLorentzVector> muonP4, electronP4;



    //
    //  MUONS
    //

    for (unsigned i = 0; i < muons.size(); i++)
    {
        TMuon* muon = muons[i];

        TLorentzVector *p4_ = (TLorentzVector*)muonP4_->ConstructedAt(i);
        copy_p4(muon, MUON_MASS, p4_);

        TLorentzVector p4;
        copy_p4(muon, MUON_MASS, p4);
        muonP4.push_back(p4);

        muonQ.push_back(muon->q);

        muonIsolation.push_back(particleSelector->GetMuonIso(muon));

        muonIsTight.push_back(particleSelector->PassMuonID(muon, cuts->tightHZZMuonID));
        muonIsLoose.push_back(particleSelector->PassMuonID(muon, cuts->looseHZZMuonID));
        muonIsIsolated.push_back(particleSelector->PassMuonIso(muon, cuts->wpHZZMuonIso));

        muonMatchIdx.push_back(-1);
        muonIsMatched.push_back(kFALSE);

        if (muonIsLoose.back())
            nLooseMuons++;
        if (muonIsTight.back())
            nTightMuons++;
    }



    //  ELECTRONS
    //

    for (unsigned i = 0; i < electrons.size(); i++)
    {
        TElectron* electron = electrons[i];

        TLorentzVector *p4_ = (TLorentzVector*)electronP4_->ConstructedAt(i);
        copy_p4(electron, ELE_MASS, p4_);

        TLorentzVector p4;
        copy_p4(electron, ELE_MASS, p4);
        electronP4.push_back(p4);

        electronQ.push_back(electron->q);

        electronIsolation.push_back(particleSelector->GetElectronIso(electron));

        electronIsTight.push_back(particleSelector->PassElectronID(electron, cuts->tightHZZElectronID));  
        electronIsLoose.push_back(particleSelector->PassElectronID(electron, cuts->looseHZZElectronID));  
        electronIsV2Iso.push_back(electron->pass2017isoV2wpHZZ);

        electronMatchIdx.push_back(-1);
        electronIsMatched.push_back(kFALSE);

        if (electronIsLoose.back())
            nLooseElectrons++;
        if (electronIsTight.back())
            nTightElectrons++;
    }

    nLooseLeptons = nLooseMuons + nLooseElectrons;
    nTightLeptons = nTightMuons + nTightElectrons;



    //
    //  MATCHING
    //


    // Muons

    if ((muonP4.size() > 0) && (dressedMuonP4.size() > 0))
    {
        for (unsigned i = 0; i < muonP4.size(); i++)
        {
            // Find DeltaR between reco lep and each gen lep
            vector<Float_t> deltaR(dressedMuonP4.size());

            for (unsigned j = 0; j < dressedMuonP4.size(); j++)
            {
                if (dressedMuonIsMatched[j])
                    deltaR[j] = 999999;
                else
                    deltaR[j] = muonP4[i].DeltaR(dressedMuonP4[j]);
            }

            // Find index of minimum DeltaR
            unsigned m = min_element(deltaR.begin(), deltaR.end()) - deltaR.begin();
            if (deltaR[m] > MATCH_DR_MAX)
                break;

            muonMatchIdx[i] = m;
            muonIsMatched[i] = muonIsTight[i];
            dressedMuonMatchIdx[m] = i;
            dressedMuonIsMatched[m] = muonIsTight[i];
        }
    }


    // Electrons

    if ((electronP4.size() > 0) && (dressedElectronP4.size() > 0))
    {
        for (unsigned i = 0; i < electronP4.size(); i++)
        {
            // Find DeltaR between reco lep and each gen lep
            vector<Float_t> deltaR(dressedElectronP4.size());

            for (unsigned j = 0; j < dressedElectronP4.size(); j++)
            {
                if (dressedElectronIsMatched[j])
                    deltaR[j] = 999999;
                else
                    deltaR[j] = electronP4[i].DeltaR(dressedElectronP4[j]);
            }

            // Find index of minimum DeltaR
            unsigned m = min_element(deltaR.begin(), deltaR.end()) - deltaR.begin();
            if (deltaR[m] > MATCH_DR_MAX)
                break;

            electronMatchIdx[i] = m;
            electronIsMatched[i] = electronIsTight[i];
            dressedElectronMatchIdx[m] = i;
            dressedElectronIsMatched[m] = electronIsTight[i];
        }
    }


    // Look for unmatched leptons

    bool foundUnmatched = kFALSE;
    for (unsigned i = 0; i < nDressedMuons; i++)
    {
        if (!dressedMuonIsMatched[i])
            foundUnmatched = kTRUE;
    }
    for (unsigned i = 0; i < nDressedElectrons; i++)
    {
        if (!dressedElectronIsMatched[i])
            foundUnmatched = kTRUE;
    }
    if (!foundUnmatched)
        isMatched = kTRUE;


    if (isMatched)
    {
        hTotalEvents->Fill(6);
        hMatchedEvents->Fill(1, weight);
        hMatchedEvents->Fill(C, weight);
        hMatchedEvents->Fill(D, weight);
    }


    // Fill tree

    outTree->Fill();
    this->passedEvents++;
    return kTRUE;
}

void EfficiencyAnalyzer::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void EfficiencyAnalyzer::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void EfficiencyAnalyzer::ReportPostTerminate()
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
    std::unique_ptr<EfficiencyAnalyzer> selector(new EfficiencyAnalyzer());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
