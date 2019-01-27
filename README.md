BLT
===

BLT is a framework for analyzing bacon ntuples originally developed by @jiafulow.

This branch is designed for the BF(Z -> 4l) measurement over the full 2017 dataset (Moriond 2018).

Setup
=====

This branch uses CMSSW_9_4_12 on cmslpc-sl6:

### Build and source environment

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
cmsrel CMSSW_9_4_12
cd CMSSW_9_4_12/src
cmsenv
git cms-init
```

### Checkout dependencies

Depends on jlrainbolt/BaconAna 2016legacy branch:

```
git clone -b 2016legacy git@github.com:jlrainbolt/BaconAna
```

### Checkout and compile BLT code

```
git clone -b trimmed_legacy git@github.com:jlrainbolt/BLT.git
scram b -j 12
```

## Running the analyzer

7 input arguments are mandatory: [input file] [no of events] [dataset] [datasetgroup] [selection] [period] [jobid]
[input file]

Currently [selection] = "emu"

```
cd BLT/BLTAnalysis/test
MultileptonAnalyzer input/DYJetsToLL_M-50.txt 1000 DYJetsToLL_M-50 DYJetsToLL emu 2016 0
```

## Running a BLT analyzer with condor

Must be on cmslpc

```
cd BLT/BLTAnalysis/test
./batch_lpc_cfg.py
```
