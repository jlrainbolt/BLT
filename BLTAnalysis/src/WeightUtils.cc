#include "BLT/BLTAnalysis/interface/WeightUtils.h"


WeightUtils::WeightUtils(string dataPeriod, string selection, bool isRealData)
{
    _dataPeriod = dataPeriod;
    _selection  = selection;
    _isRealData = isRealData;
    std::string fileName;

    const std::string cmssw_base = getenv("CMSSW_BASE");
    const std::string data_dir = cmssw_base + "/src/BLT/BLTAnalysis/data/";



    //
    //  PILEUP
    //

    fileName = data_dir + "pu_weights_" + _dataPeriod + ".root";
    TFile* puFile = new TFile(fileName.c_str(), "OPEN");
    _puReweight = (TH1*)puFile->Get("weights");



    //
    //  TRIGGERS    FIXME
    //

    // Double muon

    // Leg 1
    fileName = data_dir + "/muon_trigger/sf_Mu17Leg_Eta0to09.root";
    TFile* f_DoubleMuTrigSF_leg1_0 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_data");
    _eff_doubleMu_leg1_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_0->Get("eff_mc"); 

    fileName = data_dir + "/muon_trigger/sf_Mu17Leg_Eta09to12.root";
    TFile* f_DoubleMuTrigSF_leg1_1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_data");
    _eff_doubleMu_leg1_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_1->Get("eff_mc"); 

    fileName = data_dir + "/muon_trigger/sf_Mu17Leg_Eta12to21.root";
    TFile* f_DoubleMuTrigSF_leg1_2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_data");
    _eff_doubleMu_leg1_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_2->Get("eff_mc"); 

    fileName = data_dir + "/muon_trigger/sf_Mu17Leg_Eta21to24.root";
    TFile* f_DoubleMuTrigSF_leg1_3 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg1_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_data");
    _eff_doubleMu_leg1_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg1_3->Get("eff_mc"); 

    // Leg 2
    fileName = data_dir + "/muon_trigger/sf_Mu8Leg_Eta0to09.root";
    TFile* f_DoubleMuTrigSF_leg2_0 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[0] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_data");
    _eff_doubleMu_leg2_MC[0]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_0->Get("eff_mc"); 

    fileName = data_dir + "/muon_trigger/sf_Mu8Leg_Eta09to12.root";
    TFile* f_DoubleMuTrigSF_leg2_1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[1] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_data");
    _eff_doubleMu_leg2_MC[1]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_1->Get("eff_mc"); 

    fileName = data_dir + "/muon_trigger/sf_Mu8Leg_Eta12to21.root";
    TFile* f_DoubleMuTrigSF_leg2_2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[2] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_data");
    _eff_doubleMu_leg2_MC[2]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_2->Get("eff_mc"); 

    fileName = data_dir + "/muon_trigger/sf_Mu8Leg_Eta21to24.root";
    TFile* f_DoubleMuTrigSF_leg2_3 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleMu_leg2_DATA[3] = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_data");
    _eff_doubleMu_leg2_MC[3]   = (TGraphAsymmErrors*)f_DoubleMuTrigSF_leg2_3->Get("eff_mc"); 


    // Double electron

    // Leg 1
    fileName = data_dir + "/electron_trigger/SFs_Leg1_Ele23_HZZSelection_Tag35.root";
    TFile* f_elTrigSF_leg1 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleEle_leg1_DATA = (TH2*)f_elTrigSF_leg1->Get("EGamma_EffData2D");
    _eff_doubleEle_leg1_MC   = (TH2*)f_elTrigSF_leg1->Get("EGamma_EffMC2D"); 

    // Leg 2
    fileName = data_dir + "/electron_trigger/SFs_Leg2_Ele12_HZZSelection_Tag35.root";
    TFile* f_elTrigSF_leg2 = new TFile(fileName.c_str(), "OPEN");
    _eff_doubleEle_leg2_DATA = (TH2*)f_elTrigSF_leg2->Get("EGamma_EffData2D");
    _eff_doubleEle_leg2_MC   = (TH2*)f_elTrigSF_leg2->Get("EGamma_EffMC2D");



    //
    //  ID & RECO
    //

    // Muon HZZ ID
    fileName = data_dir + "muon_id_sf_" + _dataPeriod + ".root";
    TFile* f_hzz_muIdSF = new TFile(fileName.c_str(), "OPEN");
    _hzz_muIdSF = (TH2*) f_hzz_muIdSF->Get("FINAL");
    _hzz_muIdErr = (TH2*) f_hzz_muIdSF->Get("ERROR");


    // Electron

    // HZZ ID
    fileName = data_dir + "electron_id_sf_" + _dataPeriod + ".root";
    TFile* f_hzz_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    _hzz_eleIdSF = (TH2*) f_hzz_eleIdSF->Get("EGamma_SF2D");

    fileName = data_dir + "electron_id_sf_" + _dataPeriod + "_gap.root";
    TFile* f_hzz_eleIdSF_gap = new TFile(fileName.c_str(), "OPEN"); 
    _hzz_eleIdSF_gap = (TH2*) f_hzz_eleIdSF_gap->Get("EGamma_SF2D");


    // Reco
    fileName = data_dir + "electron_reco_sf_" + _dataPeriod + ".root";
    TFile* f_eleRecoSF = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_RECO = (TH2*) f_eleRecoSF->Get("EGamma_SF2D");

    fileName = data_dir + "electron_reco_sf_" + _dataPeriod + "_lowEt.root";
    TFile* f_eleRecoSF_lowEt = new TFile(fileName.c_str(), "OPEN"); 
    _eleSF_RECO_lowEt = (TH2*) f_eleRecoSF_lowEt->Get("EGamma_SF2D");
}




//
//  SETTERS
//

void WeightUtils::SetDataBit(bool isRealData)
{
    _isRealData = isRealData;
}


void WeightUtils::SetDataPeriod(string dataPeriod)
{
    _dataPeriod = dataPeriod;
}


void WeightUtils::SetSelection(string selection)
{
    _selection = selection;
}


void WeightUtils::SetSampleName(string sampleName)
{
    _sampleName = sampleName;
}



//
//  PILEUP
//

float WeightUtils::GetPUWeight(float nPU) const
{
    if (_isRealData)
        return 1;

    if (nPU < 0)
        return 0;

    return _puReweight->GetBinContent(_puReweight->FindBin(nPU)); 
}



//
//  MC WEIGHT
//

float WeightUtils::GetSampleWeight() const
{
    if (_isRealData)
        return 1;

    float lumi = 1, xsec = 1, ngen = 1000, eff = 1;

    if      (_dataPeriod == "2016")
        lumi = 36.42;
    else if (_dataPeriod == "2017")
        lumi = 41.37;
    else if (_dataPeriod == "2018")
        lumi = 58.83;

    if ((_selection == "ee") && (_dataPeriod == "2017"))
        eff = 0.991;

    if      (_sampleName == "DYJetsToLL_M-50")
    {
        xsec = 5765.4;

        if (_dataPeriod == "2017")
            ngen = 80924255;
    }
    else if (_sampleName == "TTJets")
    {
        xsec = 831.76;

        if (_dataPeriod == "2017")
            ngen = 15173839;
    }
    else if (_sampleName == "TTTo2L2Nu")
    {
        xsec = 87.31;

        if (_dataPeriod == "2017")
            ngen = 65899840;
    }
    xsec *= 1000.;

    return xsec * lumi * eff / ngen;
}




//
//  TRIGGERS
//

// THIS IS FOR 2016!
std::pair<float, float> WeightUtils::GetDoubleMuonTriggerEff(const baconhep::TMuon* muon, const int leg) const
{
    float effData = 1, effMC = 1;

    if (_isRealData)
        return std::make_pair(effData, effMC);

    float binningEta[] = {0, 0.9, 1.2, 2.1, 2.4};
    unsigned etaBin = 0;
    for (int i = 0; i < 4; i++)
    {
        if (fabs(muon->eta) > binningEta[i] && fabs(muon->eta) <= binningEta[i+1])
        {
            etaBin = i;
            break;
        }
    }

    if      (leg == 1)
    {
        effData = _eff_doubleMu_leg1_DATA[etaBin]->Eval(muon->pt);
        effMC   = _eff_doubleMu_leg1_MC[etaBin]->Eval(muon->pt);
    }
    else if (leg == 2)
    {
        effData = _eff_doubleMu_leg2_DATA[etaBin]->Eval(muon->pt);
        effMC   = _eff_doubleMu_leg2_MC[etaBin]->Eval(muon->pt);
    }

    return std::make_pair(effData, effMC);
}


// THIS IS FOR 2016!
std::pair<float, float> WeightUtils::GetDoubleElectronTriggerEff(const baconhep::TElectron* electron, const int leg) const
{
    float effData = 1, effMC = 1;

    if (_isRealData)
        return std::make_pair(effData, effMC);

    if      (leg == 1)
    {
        effData = _eff_doubleEle_leg1_DATA->GetBinContent(_eff_doubleEle_leg1_DATA->FindBin(electron->scEta, electron->pt));
        effMC   = _eff_doubleEle_leg1_MC->GetBinContent(_eff_doubleEle_leg1_MC->FindBin(electron->scEta, electron->pt));
    }
    else if (leg == 2)
    {
        effData = _eff_doubleEle_leg2_DATA->GetBinContent(_eff_doubleEle_leg2_DATA->FindBin(electron->scEta, electron->pt));
        effMC   = _eff_doubleEle_leg2_MC->GetBinContent(_eff_doubleEle_leg2_MC->FindBin(electron->scEta, electron->pt));
    }

    return std::make_pair(effData, effMC);
}



//
//  ID & RECO
//

float WeightUtils::GetHZZMuonIDSF(const baconhep::TMuon* muon) const
{
    if (_isRealData)
        return 1;

    TH2 *hist = _hzz_muIdSF;

    float   maxPt = hist->GetYaxis()->GetXmax();
    int     bin = hist->FindBin(muon->eta, min((double) muon->pt, 0.99 * maxPt));

    if (hist->IsBinUnderflow(bin) || hist->IsBinOverflow(bin))
        return 1;
    else
        return hist->GetBinContent(bin);
}


float WeightUtils::GetHZZElectronIDSF(const baconhep::TElectron* electron) const 
{
    if (_isRealData)
        return 1;

    TH2 *hist;
    if (electron->fiducialBits & kIsGap)
        hist = _hzz_eleIdSF_gap;
    else
        hist = _hzz_eleIdSF;

    float   maxPt = hist->GetYaxis()->GetXmax();
    int     bin = hist->FindBin(electron->scEta, min((double) electron->pt, 0.99 * maxPt));

    if (hist->IsBinUnderflow(bin) || hist->IsBinOverflow(bin))
        return 1;
    else
        return hist->GetBinContent(bin);
}

float WeightUtils::GetElectronRecoSF(const baconhep::TElectron* electron) const 
{
    if (_isRealData)
        return 1;

    TH2 *hist;
    if (electron->pt < _eleSF_RECO->GetYaxis()->GetXmin())
        hist = _eleSF_RECO_lowEt;
    else
        hist = _eleSF_RECO;

    float   maxPt = hist->GetYaxis()->GetXmax();
    int     bin = hist->FindBin(electron->scEta, min((double) electron->pt, 0.99 * maxPt));

    if (hist->IsBinUnderflow(bin) || hist->IsBinOverflow(bin))
        return 1;
    else
        return hist->GetBinContent(bin);
}
