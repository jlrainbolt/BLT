BLT
===

BLT is a framework for analyzing bacon ntuples originally developed by @jiafulow.

This branch is designed for the BF(Z -> 4l) measurement over the 8 TeV (2012) MiniAOD.

Setup
=====

This branch uses CMSSW_10_2_13 on cmslpc-sl6:

### Build and source environment

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
```

### Checkout dependencies

Depends on NWUHEP/BaconAna tag 04:

```
git clone -b 04 git@github.com:jlrainbolt/BaconAna
```

### Checkout and compile BLT code

```
git clone -b lacey_8TeV git@github.com:jlrainbolt/BLT.git
scram b -j 12
```

## Running the analyzer

7 input arguments are mandatory: [input file] [no of events] [dataset] [datasetgroup] [selection] [period] [jobid]

You can use:
* [selection] = ["loose" (>= 2 SFOS loose leptons), "tight" (>= 2 SFOS tight leptons)]
* [period] = ["2012"]

```
cd BLT/BLTAnalysis/test
MultileptonAnalyzer input/DYJetsToLL_M-50.txt 100000 DYJetsToLL_M-50 zjets_m-50 loose 2012 0
```

## Running a BLT analyzer with condor

Must be on cmslpc

```
cd BLT/BLTAnalysis/test
./batch_cfg.py
```
