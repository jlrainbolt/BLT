#include "TFile.h"
#include "TString.h"
#include "TH1.h"



using namespace std;



void calcPUScenario2017()
{
    // PU Scenario from
    //      https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py#L13
    // as recommended by
    //      https://twiki.cern.ch/twiki/bin/view/CMS/HiggsZZ4l2018#Pileup_Reweighting


    
    // Here we go, yo

    TString path = "../data/pileup/";
    TString outFileName = "pileup_2017_69200_100bins.root";
    TString dataFileName = "DataPileupHistogram2017_69200_100bins.root", dataHistName = "pileup";
    TString mcHistName = "RunIIFall17", sfHistName = "pileup_sf";



    // So what's the what's the what's the scenario?

    const Int_t NTIMES = 99;

    Double_t x[NTIMES] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
		15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 
		27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 
		39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 
		51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 
		63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 
		75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 
		87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98
    };
    
    Double_t w[NTIMES] = {
        3.39597497605e-05,
        6.63688402133e-06,
        1.39533611284e-05,
        3.64963078209e-05,
        6.00872171664e-05,
        9.33932578027e-05,
        0.000120591524486,
        0.000128694546198,
        0.000361697233219,
        0.000361796847553,
        0.000702474896113,
        0.00133766053707,
        0.00237817050805,
        0.00389825605651,
        0.00594546732588,
        0.00856825906255,
        0.0116627396044,
        0.0148793350787,
        0.0179897368379,
        0.0208723871946,
        0.0232564170641,
        0.0249826433945,
        0.0262245860346,
        0.0272704617569,
        0.0283301107549,
        0.0294006137386,
        0.0303026836965,
        0.0309692426278,
        0.0308818046328,
        0.0310566806228,
        0.0309692426278,
        0.0310566806228,
        0.0310566806228,
        0.0310566806228,
        0.0307696426944,
        0.0300103336052,
        0.0288355370103,
        0.0273233309106,
        0.0264343533951,
        0.0255453758796,
        0.0235877272306,
        0.0215627588047,
        0.0195825559393,
        0.0177296309658,
        0.0160560731931,
        0.0146022004183,
        0.0134080690078,
        0.0129586991411,
        0.0125093292745,
        0.0124360740539,
        0.0123547104433,
        0.0123953922486,
        0.0124360740539,
        0.0124360740539,
        0.0123547104433,
        0.0124360740539,
        0.0123387597772,
        0.0122414455005,
        0.011705203844,
        0.0108187105305,
        0.00963985508986,
        0.00827210065136,
        0.00683770076341,
        0.00545237697118,
        0.00420456901556,
        0.00367513566191,
        0.00314570230825,
        0.0022917978982,
        0.00163221454973,
        0.00114065309494,
        0.000784838366118,
        0.000533204105387,
        0.000358474034915,
        0.000238881117601,
        0.0001984254989,
        0.000157969880198,
        0.00010375646169,
        6.77366175538e-05,
        4.39850477645e-05,
        2.84298066026e-05,
        1.83041729561e-05,
        1.17473542058e-05,
        7.51982735129e-06,
        6.16160108867e-06,
        4.80337482605e-06,
        3.06235473369e-06,
        1.94863396999e-06,
        1.23726800704e-06,
        7.83538083774e-07,
        4.94602064224e-07,
        3.10989480331e-07,
        1.94628487765e-07,
        1.57888581037e-07,
        1.2114867431e-07,
        7.49518929908e-08,
        4.6060444984e-08,
        2.81008884326e-08,
        1.70121486128e-08,
        1.02159894812e-08
    };



    // Yes yes y'all

    Int_t maxPileupBin = 100, numPileupBins = 100;      // Binning from HZZ twiki
    TH1D *mcHist = new TH1D(mcHistName, mcHistName, numPileupBins, 0, maxPileupBin);
    mcHist->Sumw2(kFALSE);
    mcHist->FillN(NTIMES, x, w);



    // Wow how now, wow how now brown cow?

    TH1 *dataHist;

    TFile *dataFile = new TFile(path + dataFileName);
    dataFile->GetObject(dataHistName, dataHist);
    dataHist->SetDirectory(0);
    dataFile->Close();



    // Yo Mr. Busta Rhymes, tell him what I did

    TH1D *sfHist = new TH1D(*mcHist);
    sfHist->Divide(dataHist, mcHist, 1./dataHist->Integral(), 1./mcHist->Integral());
    sfHist->Sumw2(kFALSE);
    sfHist->SetName(sfHistName);
    sfHist->SetTitle(sfHistName);



    // RAARR RAARR (like a dungeon dragon)

    sfHist->Draw("*H");
    cout << "Mean SF: " << sfHist->GetMean(2) << endl;
//  for (int i = 0; i <= numPileupBins+1; i++)
//      cout << sfHist->GetBinContent(i) << endl;



    // Bust it out before the Busta bust another rhyme

//  TFile *outFile = new TFile(path + outFileName, "RECREATE");
//  sfHist->Write();
//  dataHist->Write();
//  mcHist->Write();
//  outFile->Close();



//  // Observe the vibe and check out the scenario!!

//  cout << "Wrote output to " << path + outFileName << endl;
}
