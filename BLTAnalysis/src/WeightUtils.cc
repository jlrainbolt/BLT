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
    fileName = data_dir + "MuonScaleFactors_2011_2012.root";
    TFile* f_hzz_muIdSF = new TFile(fileName.c_str(), "OPEN");
    f_hzz_muIdSF->GetObject("TH2D_ALL_2012", _hzz_muIdSF);
    _hzz_muIdSF->SetDirectory(0);
    f_hzz_muIdSF->Close();


    // Electron HZZ ID
    fileName = data_dir + "CombinedMethod_ScaleFactors_RecoIdIsoSip.root";
    TFile* f_hzz_eleIdSF = new TFile(fileName.c_str(), "OPEN"); 
    f_hzz_eleIdSF->GetObject("h_electronScaleFactor_RecoIdIsoSip", _hzz_eleIdSF);
    _hzz_eleIdSF->SetDirectory(0);
    f_hzz_eleIdSF->Close();
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


void WeightUtils::SetSampleName(string sampleName)
{
    _sampleName = sampleName;
}


void WeightUtils::SetSelection(string selection)
{
    _selection = selection;
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
//  SAMPLE WEIGHT
//

float WeightUtils::GetSampleWeight() const
{
    if (_isRealData)
        return 1;

    float xsec = 1, ngen = 1;

    if      ((_sampleName == "ZZTo4e") || (_sampleName == "ZZTo4mu") || (_sampleName == "ZZTo4tau"))
        xsec = XSEC_4x;
    else if ((_sampleName == "ZZTo2e2mu") || (_sampleName == "ZZTo2e2tau") || (_sampleName == "ZZTo2mu2tau"))
        xsec = XSEC_2x2y;
    else
        return 1;

    if      (_sampleName == "ZZTo4e")
        ngen = NGEN_4e;
    else if (_sampleName == "ZZTo4mu")
        ngen = NGEN_4mu;
    else if (_sampleName == "ZZTo4tau")
        ngen = NGEN_4tau;
    else if (_sampleName == "ZZTo2e2mu")
        ngen = NGEN_2e2mu;
    else if (_sampleName == "ZZTo2e2tau")
        ngen = NGEN_2e2tau;
    else if (_sampleName == "ZZTo2mu2tau")
        ngen = NGEN_2mu2tau;

    xsec /= XSEC_4L;
    ngen /= NGEN_4L;

    return xsec / ngen;
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

    TH2 *hist = _hzz_eleIdSF;

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
