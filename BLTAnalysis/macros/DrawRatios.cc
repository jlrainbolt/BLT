// STL
#include <iostream>
#include <vector>
#include <tuple>

// ROOT
#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRatioPlot.h"

using namespace std;
#include "/uscms/home/jrainbol/analysis/include/PlotUtils.hh"


/*
**  DrawRatios
**
**  Draws ratio plots from tag-and-probe studies
*/ 

void DrawRatios()
{

    //
    //  INPUT FILES
    //

    const TString tag = "iso35";
    const TString suffix = "muon_2017",    YEAR_STR = "2017";
    const TString EOS_PATH = "root://cmseos.fnal.gov//store/user/jrainbol/TagAndProbe/";

    TString inName  = "mumu_" + tag + "_" + suffix + ".root";
    TFile   *inFile = TFile::Open(EOS_PATH + YEAR_STR + "/" + inName);
    cout << endl << endl << "Opened " << inName << endl << endl;

    // Get directory keys
    TString dirName = "hists_" + suffix;
    vector<TString> hnames;
    TKey *histKey;
    TDirectory *dir = (TDirectory*) inFile->Get(dirName);
    TIter next(dir->GetListOfKeys());
    while ((histKey = (TKey*) next()))
    {
        TString hname = histKey->GetName();
        hnames.push_back(hname);
    }
    cout << "Got directory keys from " << inName << endl;
    hnames.erase(hnames.begin(), hnames.begin() + 2);

    const unsigned H = hnames.size();

    TString ssName = "mumuSS_" + tag + "_" + suffix + ".root";
    TFile   *ssFile = TFile::Open(EOS_PATH + YEAR_STR + "/" + ssName);
    cout << endl << endl << "Opened " << ssName << endl;



    //
    //  GET HISTOGRAMS
    //

    // Get histograms
    TH1 *mumu[H], *mumuSS[H];

    for (unsigned h = 0; h < H; h++)
    {
        inFile->GetObject(dirName + "/" + hnames[h], mumu[h]);
        mumu[h]->Rebin(2);
        mumu[h]->SetDirectory(0);

        ssFile->GetObject(dirName + "/" + hnames[h], mumuSS[h]);
        mumuSS[h]->Rebin(2);
        mumuSS[h]->SetDirectory(0);
    }
    cout << "Got histograms" << endl;

    inFile->Close();
    ssFile->Close();



    //
    //  OUTPUT FILE
    //

    TString outName = "ratios_" + tag + "_" + suffix + ".root";
    TFile *outFile  = new TFile(outName, "RECREATE");

    // Draw
    TCanvas *c[H];
    for (unsigned h = 0; h < H; h++)
    {
/*
        TString title = data[h]->GetYaxis()->GetTitle();
        TString width; width.Form("%g", data[h]->GetBinWidth(1));
        title.Prepend("Events/" + width + "");
        title.ReplaceAll("\\mbox{", " ");
        title.ReplaceAll("}", "");

        TString xtitle = data[h]->GetXaxis()->GetTitle();
        xtitle.ReplaceAll("\\ ", " ");
        xtitle.ReplaceAll("\\ell", "l");
        xtitle.ReplaceAll("\\cos", "cos ");
        xtitle.ReplaceAll("\\mbox{", "");
        xtitle.ReplaceAll("})", ")");
        xtitle.ReplaceAll("}_", "_");
        xtitle.ReplaceAll("  ", " ");

        data[h]->SetXTitle(xtitle);
        bkg[h]->SetXTitle(xtitle);
        res[h]->SetXTitle(xtitle);
        gen[h]->SetXTitle(xtitle);
        reco[h]->SetXTitle(xtitle);
        sel[h]->SetXTitle(xtitle);
        ps[h]->SetXTitle(xtitle);
*/
        

        // Reco vs. gen

        c[h] = new TCanvas(hnames[h], "", lCanvasSize, lCanvasSize);
//      Facelift(c);

        mumu[h]->SetStats(0);
        mumu[h]->SetLineColor(kBlack);
        mumu[h]->SetLineWidth(2);
        mumu[h]->SetMarkerStyle(20);
        mumu[h]->SetMarkerSize(2);

//      mumuSS[h]->SetMaximum(1.25 * gen[h]->GetMaximum());
        mumuSS[h]->SetStats(0);
        mumuSS[h]->SetLineColor(kBlue);
        mumuSS[h]->SetMarkerColor(kBlue);
        mumuSS[h]->SetLineWidth(2);
        mumuSS[h]->SetLineWidth(2);
        mumuSS[h]->SetMarkerStyle(22);
        mumuSS[h]->SetMarkerSize(2);

        TRatioPlot *r = new TRatioPlot(mumu[h], mumuSS[h], "diff");
        r->SetH1DrawOpt("E");
        r->SetH2DrawOpt("E");
        r->SetSeparationMargin(0.01);

        float LeftPosition = 0.5,       LeftMargin = 2. * lCanvasMargin - lLegendMargin;
        float RightPosition = 1,        RightMargin = -lLegendMargin;
        float TopPosition = 1,          TopMargin = -lLegendMargin;
        float BottomPosition = TopPosition - 0.085 * 3;
        float BottomMargin = 2. * lCanvasMargin - lLegendMargin;

        TLegend *l = new TLegend(LeftPosition + LeftMargin, BottomPosition - TopMargin,
                TopPosition + TopMargin, TopPosition + TopMargin);
        l->AddEntry(mumu[h], "Opp. sign #mu#mu", "LP");
        l->AddEntry(mumuSS[h], "Same sign #mu#mu", "LP");
//      Facelift(l);

        c[h]->cd();
        r->Draw();
        r->GetUpperPad()->cd();
        l->Draw();

        r->GetLowerRefYaxis()->SetTitle("OS - SS");
/*
        Facelift(r_reco_gen->GetUpperRefYaxis());
        r_reco_gen->GetUpperRefYaxis()->SetTitle(title);
        r_reco_gen->GetUpperRefYaxis()->SetTitleOffset(lTitleOffsetY);
        r_reco_gen->GetUpperPad()->Modified();

        Facelift(r_reco_gen->GetLowerRefXaxis());
        Facelift(r_reco_gen->GetLowerRefYaxis());
        r_reco_gen->GetLowerRefGraph()->SetMinimum(0.8);
        r_reco_gen->GetLowerRefGraph()->SetMaximum(1.2);
        r_reco_gen->SetLowBottomMargin(3 * lCanvasMargin);
        r_reco_gen->SetLeftMargin(1.2 * lCanvasMargin);
        r_reco_gen->GetLowerPad()->Modified();

        r_reco_gen->GetUpperPad()->cd();
        l_reco_gen->Draw();
        c_reco_gen->SaveAs(".pdf");
*/
        c[h]->Write();
    }

    outFile->Close();
    cout << "Wrote ratios to " << outName << endl;
}
