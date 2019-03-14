#include "SelectedEventRetriever.h"

//
// See header file for class documentation
//

using namespace baconhep;
using namespace std;


SelectedEventRetriever::SelectedEventRetriever() : BLTSelector()
{

}

SelectedEventRetriever::~SelectedEventRetriever()
{

}

void SelectedEventRetriever::Begin(TTree *tree)
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


    //
    //  SELECTED EVENTS
    //

    // Open appropriate "selected" tree
    TString tempFileName = "root://cmseos.fnal.gov//store/user/jrainbol/Selected/" + params->period
                + "/background_" + params->datasetgroup + ".root";
//              + "/selected_" + params->datasetgroup + ".root";
    TString tempTreeName = params->selection + "_" + params->datasetgroup;

    TFile *tempFile = TFile::Open(tempFileName);
    TTreeReader tempTree(tempTreeName, tempFile);


    // Get run, lumi, event numbers
    TTreeReaderValue    <Int_t>     runNum_     (tempTree,  "runNum");
    TTreeReaderValue    <Int_t>     evtNum_     (tempTree,  "evtNum");
    TTreeReaderValue    <Int_t>     lumiSec_    (tempTree,  "lumiSec");

    nEvents = tempTree.GetEntries(kTRUE);

    cout << endl << "SELECTED EVENTS" << endl;
    cout << "Run, lumi, event:" << endl;
    while (tempTree.Next())
    {
        cout << *runNum_ << ", " << *lumiSec_ << ", " << *evtNum_ << endl;
        runNumber.push_back(*runNum_);
        lumiSection.push_back(*lumiSec_);
        evtNumber.push_back(*evtNum_);
    }
    cout << endl;
    tempFile->Close();


    //
    //  OUTPUT TREE
    //

    string outFileName = params->get_output_filename("output");

    outFile = new TFile(outFileName.c_str(), "RECREATE");
    outFile->cd();
    outTree = tree->CloneTree(0);


    ReportPostBegin();
}


Bool_t SelectedEventRetriever::Process(Long64_t entry)
{


    //
    //  START
    //

    GetEntry(entry, 1);  // load all branches
    this->totalEvents++;

    if (entry%10000==0)  
        std::cout << "... Processing event " << entry << " Run: " << fInfo->runNum 
                  << " Lumi: " << fInfo->lumiSec << " Event: " << fInfo->evtNum 
                  << std::endl;



    //
    //  SELECTION
    //

    // Check event against selection
    bool isMatch = kFALSE;

    for (unsigned i = 0; i < nEvents; i++)
    {
        if (
                (fInfo->runNum == runNumber[i])
            &&  (fInfo->evtNum == evtNumber[i])
            &&  (fInfo->lumiSec == lumiSection[i])
           )
        {
            isMatch = kTRUE;
            break;
        }
    }


    // Fill output tree with matched events
    if (isMatch)
    {
        outTree->Fill();
        this->passedEvents++;
    }
    return kTRUE;
}

void SelectedEventRetriever::Terminate()
{
    outFile->Write();
    outFile->Close();

    ReportPostTerminate();
}

void SelectedEventRetriever::ReportPostBegin()
{
    std::cout << "  ==== Begin Job =============================================" << std::endl;
    std::cout << *params << std::endl;
    std::cout << "  ============================================================" << std::endl;
}

void SelectedEventRetriever::ReportPostTerminate()
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
    std::unique_ptr<SelectedEventRetriever> selector(new SelectedEventRetriever());

    try {
        selector->MakeMeSandwich(argc, argv);  //<===the real main function is here

    } catch (const std::exception& e) {
        std::cerr << "An exception is caught: " << e.what() << std::endl;
        throw;

        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
