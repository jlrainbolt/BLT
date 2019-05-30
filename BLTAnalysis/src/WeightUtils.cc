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
    puFile->GetObject("weights", _puReweight);
    _puReweight->SetDirectory(0);
    puFile->Close();



    //
    //  ID & RECO
    //

    // Muon HZZ ID
    fileName = data_dir + "muon_id_sf_" + _dataPeriod + ".root";
    TFile* f_hzz_muIdSF = new TFile(fileName.c_str(), "OPEN");
    f_hzz_muIdSF->GetObject("FINAL", _hzz_muIdSF);
    _hzz_muIdSF->SetDirectory(0);
    f_hzz_muIdSF->Close();


    // Electron

    // HZZ ID
    fileName = data_dir + "electron_id_sf_" + _dataPeriod + ".root";
    TFile* f_hzz_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    f_hzz_eleIdSF->GetObject("EGamma_SF2D", _hzz_eleIdSF);
    _hzz_eleIdSF->SetDirectory(0);
    f_hzz_eleIdSF->Close();

    fileName = data_dir + "electron_id_sf_" + _dataPeriod + "_gap.root";
    TFile* f_hzz_eleIdSF_gap = new TFile(fileName.c_str(), "OPEN"); 
    f_hzz_eleIdSF_gap->GetObject("EGamma_SF2D", _hzz_eleIdSF_gap);
    _hzz_eleIdSF_gap->SetDirectory(0);
    f_hzz_eleIdSF_gap->Close();


    // Reco
    fileName = data_dir + "electron_reco_sf_" + _dataPeriod + ".root";
    TFile* f_eleRecoSF = new TFile(fileName.c_str(), "OPEN"); 
    f_eleRecoSF->GetObject("EGamma_SF2D", _eleSF_RECO);
    _eleSF_RECO->SetDirectory(0);
    f_eleRecoSF->Close();

    if (_dataPeriod != "2018")
    {
        fileName = data_dir + "electron_reco_sf_" + _dataPeriod + "_lowEt.root";
        TFile* f_eleRecoSF_lowEt = new TFile(fileName.c_str(), "OPEN"); 
        f_eleRecoSF_lowEt->GetObject("EGamma_SF2D", _eleSF_RECO_lowEt);
        _eleSF_RECO_lowEt->SetDirectory(0);
        f_eleRecoSF_lowEt->Close();
    }
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
//  ID & RECO
//

float WeightUtils::GetHZZMuonIDSF(const baconhep::TMuon* muon) const
{
    if (_isRealData)
        return 1;

    TH2 *hist = _hzz_muIdSF;

    float   minPt = hist->GetYaxis()->GetXmin();
    float   maxPt = hist->GetYaxis()->GetXmax();

    double  pt = muon->pt;
    pt = min(pt, 0.99 * maxPt);
    pt = max(pt, 1.01 * minPt);

    int bin = hist->FindBin(muon->eta, pt);

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

    float   minPt = hist->GetYaxis()->GetXmin();
    float   maxPt = hist->GetYaxis()->GetXmax();

    double  pt = electron->pt;
    pt = min(pt, 0.99 * maxPt);
    pt = max(pt, 1.01 * minPt);

    int bin = hist->FindBin(electron->scEta, pt);

    if (hist->IsBinUnderflow(bin) || hist->IsBinOverflow(bin))
        return 1;
    else
        return hist->GetBinContent(bin);
}

float WeightUtils::GetElectronRecoSF(const baconhep::TElectron* electron) const 
{
    if (_isRealData)
        return 1;

    TH2 *hist = _eleSF_RECO;

    float   minPt = hist->GetYaxis()->GetXmin();
    float   maxPt = hist->GetYaxis()->GetXmax();

    double  pt = electron->pt;

    if ((pt < minPt) && (_dataPeriod != "2018"))
    {
        hist = _eleSF_RECO_lowEt;
        minPt = hist->GetYaxis()->GetXmin();
        maxPt = hist->GetYaxis()->GetXmax();
    }

    pt = min(pt, 0.99 * maxPt);
    pt = max(pt, 1.01 * minPt);

    int bin = hist->FindBin(electron->scEta, pt);

    if (hist->IsBinUnderflow(bin) || hist->IsBinOverflow(bin))
        return 1;
    else
        return hist->GetBinContent(bin);
}
