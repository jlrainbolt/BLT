#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg         = bm.JobConfig
selection   = 'loose'
period      = '2012'
executable  = 'execBatch.sh'
location    = 'lpc'

data_samples    = ['double_mu', 'double_el', 'single_mu', 'single_el']
#data_samples    = []
mc_samples      = ['signal', 'zjets', 'ttbar', 'diboson', 'higgs']
#mc_samples      = []



''' 
Set job configurations.  
'''

path = '/eos/uscms/store/group/lpcbacon/04'


# DATA #
data_dict = {}

data_dict['double_mu'] = \
[
    cfg(data_name = 'DoubleMuParked_Run2012A',
        path     = '{0}/DoubleMuParked_2012A-22Jan2013'.format(path),
        nJobs    = 1,
        suffix   = 'muon_2012'
        ),
    cfg(data_name = 'DoubleMuParked_Run2012B',
        path     = '{0}/DoubleMuParked_2012B-22Jan2013'.format(path),
        nJobs    = 8,
        suffix   = 'muon_2012'
        ),
    cfg(data_name = 'DoubleMuParked_Run2012C',
        path     = '{0}/DoubleMuParked_2012C-22Jan2013'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2012'
        ),
    cfg(data_name = 'DoubleMuParked_Run2012D',
        path     = '{0}/DoubleMuParked_2012D-22Jan2013'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2012'
        ),
]

data_dict['single_mu'] = \
[
    cfg(data_name = 'SingleMu_Run2012A',
        path     = '{0}/SingleMu_2012A-22Jan2013'.format(path),
        nJobs    = 1,
        suffix   = 'muon_2012'
        ),
    cfg(data_name = 'SingleMu_Run2012B',
        path     = '{0}/SingleMu_2012B-22Jan2013'.format(path),
        nJobs    = 7,
        suffix   = 'muon_2012'
        ),
    cfg(data_name = 'SingleMu_Run2012C',
        path     = '{0}/SingleMu_2012C-22Jan2013'.format(path),
        nJobs    = 11,
        suffix   = 'muon_2012'
        ),
    cfg(data_name = 'SingleMu_Run2012D',
        path     = '{0}/SingleMu_2012D-22Jan2013'.format(path),
        nJobs    = 13,
        suffix   = 'muon_2012'
        ),
]

data_dict['double_el'] = \
[
    cfg(data_name = 'DoubleElectron_Run2012A',
        path     = '{0}/2012_Data_Multicrab_DoubleElectron_Run2012A-22Jan2013-v1/161006_183851/0000'.format(path),
        nJobs    = 4,
        suffix   = 'electron_2012'
        ),
    cfg(data_name = 'DoubleElectron_Run2012B',
        path     = '{0}/2012_Data_Multicrab_DoubleElectron_Run2012B-22Jan2013-v1/161006_184140/0000'.format(path),
        nJobs    = 8,
        suffix   = 'electron_2012'
        ),
    cfg(data_name = 'DoubleElectron_Run2012C',
        path     = '{0}/2012_Data_Multicrab_DoubleElectron_Run2012C-22Jan2013-v1/161006_184501/0000'.format(path),
        nJobs    = 12,
        suffix   = 'electron_2012'
        ),
    cfg(data_name = 'DoubleElectron_Run2012D',
        path     = '{0}/2012_Data_Multicrab_DoubleElectron_Run2012D-22Jan2013-v1/161006_184923/0000'.format(path),
        nJobs    = 12,
        suffix   = 'electron_2012'
        ),
]

data_dict['single_mu'] = \
[
    cfg(data_name = 'SingleElectron_Run2012A',
        path     = '{0}/SingleElectron_2012A-22Jan2013'.format(path),
        nJobs    = 3,
        suffix   = 'electron_2012'
        ),
    cfg(data_name = 'SingleElectron_Run2012B',
        path     = '{0}/SingleElectron_2012B-22Jan2013'.format(path),
        nJobs    = 14,
        suffix   = 'electron_2012'
        ),
    cfg(data_name = 'SingleElectron_Run2012C',
        path     = '{0}/SingleElectron_2012C-22Jan2013'.format(path),
        nJobs    = 20,
        suffix   = 'electron_2012'
        ),
    cfg(data_name = 'SingleElectron_Run2012D',
        path     = '{0}/SingleElectron_2012D-22Jan2013'.format(path),
        nJobs    = 20,
        suffix   = 'electron_2012'
        ),
]


# MONTE CARLO #
mc_dict = {}

mc_dict['signal'] = \
[
    cfg(data_name = 'ZZTo4mu',
        path     = '{0}/Summer12_ZZTo4mu'.format(path),
        nJobs    = 3,
        suffix   = 'zz_4l'
        ),
    cfg(data_name = 'ZZTo2e2mu',
        path     = '{0}/Summer12_ZZTo2e2mu'.format(path),
        nJobs    = 2,
        suffix   = 'zz_4l'
        ),
    cfg(data_name = 'ZZTo4e',
        path     = '{0}/Summer12_ZZTo4e'.format(path),
        nJobs    = 2,
        suffix   = 'zz_4l'
        ),
    cfg(data_name = 'ZZTo4tau',
        path     = '{0}/Summer12_ZZTo4tau'.format(path),
        nJobs    = 1,
        suffix   = 'zz_4l'
        ),
    cfg(data_name = 'ZZTo2e2tau',
        path     = '{0}/Summer12_ZZTo2e2tau'.format(path),
        nJobs    = 1,
        suffix   = 'zz_4l'
        ),
    cfg(data_name = 'ZZTo2mu2tau',
        path     = '{0}/Summer12_ZZTo2mu2tau'.format(path),
        nJobs    = 1,
        suffix   = 'zz_4l'
        ),
]

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/Summer12_DYJetsToLL_M-50_TuneZ2Star'.format(path),
        nJobs    = 8,
        suffix   = 'zjets_m-50'
        ),
#   cfg(data_name = 'DYJetsToLL_M-10to50',
#       path     = '{0}/Summer12_DYJetsToLL_M-10To50_TuneZ2Star'.format(path),
#       nJobs    = 9,
#       suffix   = 'zjets_m-10'
#       ),
]

mc_dict['ttbar'] = \
[
    cfg(data_name = 'TTJets',
        path     = '{0}/Summer12_TTJets_FullLeptMGDecays'.format(path),
        nJobs    = 6,
        suffix   = 'ttbar'
        ),
    cfg(data_name = 'TTZJets',
        path     = '{0}/Summer12_TTZJets'.format(path),
        nJobs    = 1,
        suffix   = 'ttz_2l2nu'
        ),
]

mc_dict['diboson'] = \
[
    cfg(data_name = 'WW',
        path     = '{0}/Summer12_WW_TuneZ2star'.format(path),
        nJobs    = 3,
        suffix   = 'ww_2l2nu'
        ),
    cfg(data_name = 'WZJetsTo2L2Q',
        path     = '{0}/Summer12_WZJetsTo2L2Q_TuneZ2star'.format(path),
        nJobs    = 1,
        suffix   = 'wz_2l2q'
        ),
    cfg(data_name = 'WZTo3LNu',
        path     = '{0}/Summer12_WZJetsTo3LNu_TuneZ2'.format(path),
        nJobs    = 1,
        suffix   = 'wz_3lnu'
        ),
    cfg(data_name = 'ZZJetsTo2L2Nu',
        path     = '{0}/Summer12_ZZJetsTo2L2Nu_TuneZ2star'.format(path),
        nJobs    = 1,
        suffix   = 'zz_2l2nu'
        ),
    cfg(data_name = 'ZZJetsTo2L2Q',
        path     = '{0}/Summer12_ZZJetsTo2L2Q_TuneZ2star'.format(path),
        nJobs    = 1,
        suffix   = 'zz_2l2q'
        ),
    cfg(data_name = 'ZZ',
        path     = '{0}/Summer12_ZZ_TuneZ2star'.format(path),
        nJobs    = 3,
        suffix   = 'zz'
        ),
]

mc_dict['higgs'] = \
[
#   cfg (data_name = 'GluGluHToZZTo4L',
#       path     = '{0}/ggH_zz_4l'.format(path),
#       nJobs    = 3,
#       suffix   = 'ggH_zz_4l'
#      ),
#   cfg(data_name = 'VBF_HToZZTo4L',
#       path     = '{0}/vbfH_zz_4l'.format(path),
#       nJobs    = 2,
#       suffix   = 'vbfH_zz_4l'
#      ),
]



batch_list = []
batch_list += sum([data_dict[n] for n in data_samples], []) 
batch_list += sum([mc_dict[n] for n in mc_samples], []) 

batch = bm.BatchMaster(config_list  = batch_list, 
                      stage_dir     = 'batch',
                      selection     = selection,
                      period        = period,
                      executable    = executable,
                      location      = 'lpc'
                     )
batch.submit_to_batch()
