#include "TFile.h"
#include "TString.h"
#include "TH1.h"



using namespace std;



void makePUWeights()
{

    // File info
    TString path = "../data/";
    TString outFileName = "pu_weights_2012.root";
    TString dataFileName = "puProfile_Data_8TeV.root",  mcFileName = "puProfile_Summer12_53X.root";

    TString outHistName = "weights",    inHistName = "pileup";


    // Get histograms from files
    TH1 *dataHist, *mcHist;

    TFile *dataFile = new TFile(path + dataFileName);
    dataFile->GetObject(inHistName, dataHist);
    dataHist->SetDirectory(0);
    dataFile->Close();

    TFile *mcFile = new TFile(path + mcFileName);
    mcFile->GetObject(inHistName, mcHist);
    mcHist->SetDirectory(0);
    mcFile->Close();


    // Create output histogram
    TH1 *outHist = (TH1*) dataHist->Clone();
    outHist->Reset();

    outHist->Divide(dataHist, mcHist, 1./dataHist->Integral(), 1./mcHist->Integral());
    outHist->Sumw2(kFALSE);
    outHist->SetName(outHistName);


    // Save
    TFile *outFile = new TFile(path + outFileName, "RECREATE");
    outFile->cd();
    outHist->Write();
    outFile->Close();
}
