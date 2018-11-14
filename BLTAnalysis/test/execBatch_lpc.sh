#!/bin/sh

echo "Job submitted on host `hostname` on `date`"
echo ">>> arguments: $@"

### Required parameters #####
DATANAME=$1
COUNT=$2
SUFFIX=$3
SELECTION=$4
PERIOD=$5

### Transfer files, prepare directory ###
TOPDIR=$PWD

# lpc
export SCRAM_ARCH=slc6_amd64_gcc530
export CMSSW_VERSION=CMSSW_8_0_20
source /cvmfs/cms.cern.ch/cmsset_default.sh

# temporary fix
tar -xzf source.tar.gz
mv $CMSSW_VERSION tmp
scram project CMSSW $CMSSW_VERSION
cp -r tmp/src/* $CMSSW_VERSION/src
cd $CMSSW_VERSION/src
eval `scram runtime -sh`

# this used to work, now it don't
#tar -xzf source.tar.gz
#cd $CMSSW_VERSION/src/
#scramv1 b ProjectRename

cmsenv
scramv1 b -j8 #ProjectRename
cd BLT/BLTAnalysis/test
cp $TOPDIR/input_${DATANAME}_${COUNT}.txt input.txt

echo $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
echo $PATH
pwd
cat input.txt

### Run the analyzer
MultileptonAnalyzer input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT
#FakeSelector input.txt -1 $DATANAME $SUFFIX $SELECTION $PERIOD $COUNT

### Copy output and cleanup ###
cp output_${DATANAME}_${COUNT}.root ${_CONDOR_SCRATCH_DIR}
