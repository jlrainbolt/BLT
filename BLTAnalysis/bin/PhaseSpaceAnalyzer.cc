#include "PhaseSpaceAnalyzer.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;

bool P4SortCondition(TLorentzVector p1, TLorentzVector p2) {return (p1.Pt() > p2.Pt());} 


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

    // Weight utilities
    weights.reset(new WeightUtils(params->period, params->selection, false));

    const bool isSignal     = params->datasetgroup == "zz_4l";

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
    if (isSignal)
        outTree->Branch(    "sampleName",               &sampleName);
    outTree->Branch(    "runNumber",                &runNumber);
    outTree->Branch(    "evtNumber",                &evtNumber);
    outTree->Branch(    "lumiSection",              &lumiSection);
    outTree->Branch(    "genWeight",                &genWeight);
    outTree->Branch(    "decayChannel",             &decayChannel);
    outTree->Branch(    "isFiducial",               &isFiducial);

//  outTree->Branch(    "nomWeight",                &nomWeight);
//  outTree->Branch(    "qcdID",                    &qcdID);
//  outTree->Branch(    "qcdWeight",                &qcdWeight);
//  outTree->Branch(    "pdfID",                    &pdfID);
//  outTree->Branch(    "pdfWeight",                &pdfWeight);


    // Counters
    outTree->Branch(    "nMuons",                   &nMuons);
    outTree->Branch(    "nElectrons",               &nElectrons);
    outTree->Branch(    "nLeptons",                 &nLeptons);
    outTree->Branch(    "nZs",                      &nZs);


    // Leptons
    outTree->Branch(    "muonP4",                   &muonP4_,       32000,      1);
    outTree->Branch(    "muonQ",                    &muonQ);
    outTree->Branch(    "muonMother",               &muonMother);
    outTree->Branch(    "muonZIndex",               &muonZIndex);

    outTree->Branch(    "electronP4",               &electronP4_,   32000,      1);
    outTree->Branch(    "electronQ",                &electronQ);
    outTree->Branch(    "electronMother",           &electronMother);
    outTree->Branch(    "electronZIndex",           &electronZIndex);

    outTree->Branch(    "leptonsP4",                &leptonsP4);


    // Z bosons
    outTree->Branch(    "zP4",                      &zP4_,          32000,      1);
    outTree->Branch(    "zStatus",                  &zStatus);
    outTree->Branch(    "zIndex",                   &zIndex);
    outTree->Branch(    "zsP4",                     &zsP4);



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
//  fLHEWeightArr   = 0;

    fChain->SetBranchAddress(   "Info",         &fInfo,             &b_Info);
    fChain->SetBranchAddress(   "GenEvtInfo",   &fGenEvtInfo,       &b_GenEvtInfo);
    fChain->SetBranchAddress(   "GenParticle",  &fGenParticleArr,   &b_GenParticleArr);
//  fChain->SetBranchAddress(   "LHEWeight",    &fLHEWeightArr,     &b_LHEWeightArr);
}

Bool_t PhaseSpaceAnalyzer::Process(Long64_t entry)
{

    //
    //  CLEAR CONTAINERS
    //

    isFiducial = kFALSE;     // qcdID.clear();              qcdWeight.clear();
    decayChannel = 0;        // pdfID.clear();              pdfWeight.clear();
    genWeight = 1;              nomWeight = 1;

    nMuons = 0;                 nElectrons = 0;             nLeptons = 0;               nZs = 0;

    muonP4_->Delete();          muonQ.clear();              muonMother.clear();
    electronP4_->Delete();      electronQ.clear();          electronMother.clear();
    muonZIndex.clear();         electronZIndex.clear();     leptonsP4.Clear();

    zP4_->Delete();             zStatus.clear();            zIndex.clear();             zsP4.Clear();



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

    const bool isData       = (fInfo->runNum != 1);
    const bool isSignal     = params->datasetgroup == "zz_4l";
    const bool isDrellYan   = params->datasetgroup == "zjets_m-50";

    weights->SetSampleName(params->dataset);
    sampleName = params->dataset;

    // Reject (accidental?) data events
    if (isData)
        return kTRUE;



    //
    //  EVENT INFO
    //

    if (isSignal)
        genWeight = weights->GetSampleWeight();

    runNumber   = fInfo->runNum;
    evtNumber   = fInfo->evtNum;
    lumiSection = fInfo->lumiSec;
/*
    for (int i = 0; i < fLHEWeightArr->GetEntries(); i++)
    {
        TLHEWeight* lhe = (TLHEWeight*) fLHEWeightArr->At(i);

        // QCD scales
        if (i < 9)
        {
            qcdID.push_back(lhe->id);
            qcdWeight.push_back(lhe->weight / nomWeight);
        }

        // PDF sets
        else
        {
            pdfID.push_back(lhe->id);
            pdfWeight.push_back(lhe->weight / nomWeight);
        }
    }
*/




    ////
    ////
    ////    PARTICLE LOOP
    ////
    ////


    hPhaseSpaceEvents->Fill(1, genWeight);

    vector<TLorentzVector> muonP4, elecP4;
    leptonsP4 = TLorentzVector();

    for (int i = 0; i < fGenParticleArr->GetEntries(); i++)
    {
        TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(i);
        int motherIndex = particle->parent;

        if (motherIndex < 0)            // seg faults are bad!
            continue;

        TGenParticle* mother = (TGenParticle*) fGenParticleArr->At(motherIndex);
        int motherID = mother->pdgId;

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

        int charge = -1 * copysign(1, particle->pdgId);

        if      (abs(particle->pdgId) == 13)
        {
            TLorentzVector *p4_ = (TLorentzVector*) muonP4_->ConstructedAt(nMuons);
            copy_p4(particle, MUON_MASS, p4_);

            TLorentzVector p4;
            copy_p4(particle, MUON_MASS, p4);

            leptonsP4 = leptonsP4 + p4;
            muonP4.push_back(p4);
            muonQ.push_back(charge);
            muonMother.push_back(motherID);
            muonZIndex.push_back(motherIndex);
            nMuons++;
        }
        else if (abs(particle->pdgId) == 11)
        {
            TLorentzVector *p4_ = (TLorentzVector*) electronP4_->ConstructedAt(nElectrons);
            copy_p4(particle, ELE_MASS, p4_);

            TLorentzVector p4;
            copy_p4(particle, ELE_MASS, p4);

            leptonsP4 = leptonsP4 + p4;
            elecP4.push_back(p4);
            electronQ.push_back(charge);
            electronMother.push_back(motherID);
            electronZIndex.push_back(motherIndex);
            nElectrons++;
        }
    } // END particle loop

    nLeptons = nMuons + nElectrons;



    //
    //  Z COUNTING
    //

    zsP4 = TLorentzVector();

    for (unsigned i = 0; i < nMuons; i++)
    {
        // If this muon's mother Z has not already been found
        if (find(zIndex.begin(), zIndex.end(), muonZIndex[i]) == zIndex.end())
        {
            zIndex.push_back(muonZIndex[i]);
            nZs++;
        }
    }
    for (unsigned i = 0; i < nElectrons; i++)
    {
        // If this electron's mother Z has not already been found
        if (find(zIndex.begin(), zIndex.end(), electronZIndex[i]) == zIndex.end())
        {
            zIndex.push_back(electronZIndex[i]);
            nZs++;
        }
    }

    // Get the (unique) Zs from GenParticleArray
    for (unsigned i = 0; i < nZs; i++)
    {
        TGenParticle* particle = (TGenParticle*) fGenParticleArr->At(zIndex[i]);

        TLorentzVector *p4_ = (TLorentzVector*) zP4_->ConstructedAt(i);
        copy_p4(particle, p4_);

        TLorentzVector p4;
        copy_p4(particle, p4);
        zsP4 = zsP4 + p4;
        zStatus.push_back(particle->status);
    }






    ////
    ////
    ////    SELECTION
    ////
    ////


    if (params->selection != "dressed")
        return kTRUE;


    // Require correct number of leptons

    if (isSignal && (nLeptons != 4))
        return kTRUE;
    if (isDrellYan && (nLeptons != 2))
        return kTRUE;
    hTotalEvents->Fill(2);



    // Total mass

    if (leptonsP4.M() < M_MIN || leptonsP4.M() > M_MAX)
        return kTRUE;
    hTotalEvents->Fill(3);



    // Dilepton mass

    for (unsigned j = 1; j < nMuons; j++)
    {
        for (unsigned i = 0; i < j; i++)
        {
            if (muonQ[i] != muonQ[j])
            {
                TLorentzVector dimuonP4 = muonP4[i] + muonP4[j];

                if (dimuonP4.M() < MLL_MIN)
                    return kTRUE;
            }
        }
    }
    for (unsigned j = 1; j < nElectrons; j++)
    {
        for (unsigned i = 0; i < j; i++)
        {
            if (electronQ[i] != electronQ[j])
            {
                TLorentzVector dielecP4 = elecP4[i] + elecP4[j];

                if (dielecP4.M() < MLL_MIN)
                    return kTRUE;
            }
        }
    }
    hTotalEvents->Fill(4);



    // Categorize

    TLorentzVector muonsP4, elecsP4;

    for (unsigned i = 0; i < nMuons; i++)
        muonsP4 = muonsP4 + muonP4[i];
    for (unsigned i = 0; i < nElectrons; i++)
        elecsP4 = elecsP4 + elecP4[i];

    hPhaseSpaceEvents->Fill(1, genWeight);
    unsigned C = 0;                             // Index

    if      (nMuons == 2 && nElectrons == 0)    // mumu = 3
        C = 3;

    else if (nMuons == 0 && nElectrons == 2)    // ee   = 4
        C = 4;

    else if (nMuons == 4 && nElectrons == 0)    // 4m   = 6
        C = 6;

    else if (nMuons == 2 && nElectrons == 2     // 2m2e = 7
            && muonsP4.M() > elecsP4.M())
        C = 7;

    else if (nMuons == 2 && nElectrons == 2     // 2e2m = 8
            && muonsP4.M() < elecsP4.M())
        C = 8;

    else if (nMuons == 0 && nElectrons == 4)    // 4e   = 9
        C = 9;

    unsigned D = (C < 6) ? 2 : 5;
    hPhaseSpaceEvents->Fill(C, genWeight);
    hPhaseSpaceEvents->Fill(D, genWeight);
    decayChannel = C;



    // Fiducial acceptance

    std::vector<TLorentzVector> sorted_leps = muonP4;
    sorted_leps.insert(sorted_leps.end(), elecP4.begin(), elecP4.end());
    sort(sorted_leps.begin(), sorted_leps.end(), P4SortCondition);

    isFiducial = kTRUE;

    for (unsigned i = 0; i < sorted_leps.size(); i++)
    {
        if (fabs(sorted_leps[i].Eta()) > ETA_MAX)
            isFiducial = kFALSE;
    }

    if (sorted_leps[0].Pt() < PT1_MIN)
        isFiducial = kFALSE;
    if (sorted_leps[1].Pt() < PT2_MIN)
        isFiducial = kFALSE;

    if (nLeptons > 2)
    {
        if (sorted_leps[2].Pt() < PT_MIN)
            isFiducial = kFALSE;
        if (sorted_leps[3].Pt() < PT_MIN)
            isFiducial = kFALSE;
    }



    // Fill tree

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
