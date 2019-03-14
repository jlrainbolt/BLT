#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg         = bm.JobConfig
selection   = 'loose'
period      = '2016'
executable  = 'execBatch.sh'
location    = 'lpc'

data_samples    = ['double_mu', 'double_eg']
mc_samples      = ['signal', 'zjets', 'ttbar', 'diboson', 'higgs']



''' 
Set job configurations.  
'''

# DATA #
data_dict = {}

path = '/eos/uscms/store/user/jrainbol/Bacon/2016'
data_dict['double_mu'] = \
[
    cfg(data_name = 'DoubleMuon_Run2016B_v1',
        path     = '{0}/DoubleMuon_Run2016B-17Jul2018_ver1-v1'.format(path),
        nJobs    = 1,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016B_v2',
        path     = '{0}/DoubleMuon_Run2016B-17Jul2018_ver2-v1'.format(path),
        nJobs    = 15,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016C_v1',
        path     = '{0}/DoubleMuon_Run2016C-17Jul2018-v1'.format(path),
        nJobs    = 5,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016D_v1',
        path     = '{0}/DoubleMuon_Run2016D-17Jul2018-v1'.format(path),
        nJobs    = 5,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016E_v1',
        path     = '{0}/DoubleMuon_Run2016E-17Jul2018-v1'.format(path),
        nJobs    = 5,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016F_v1',
        path     = '{0}/DoubleMuon_Run2016F-17Jul2018-v1'.format(path),
        nJobs    = 5,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016G_v1',
        path     = '{0}/DoubleMuon_Run2016G-17Jul2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016H_v1',
        path     = '{0}/DoubleMuon_Run2016H-17Jul2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2016'
        ),
]

path = '/eos/uscms/store/user/jrainbol/BaconProd/eg_2016'
data_dict['double_eg'] = \
[
    cfg(data_name = 'DoubleEG_Run2016B_v1',
        path     = '{0}/RunB_ver1'.format(path),
        nJobs    = 1,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016B_v2',
        path     = '{0}/RunB_ver2'.format(path),
        nJobs    = 25,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016B_v2_nfl',
        path     = '{0}/RunB_ver2_nfl'.format(path),
        nJobs    = 10,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016C_v1',
        path     = '{0}/RunC'.format(path),
        nJobs    = 10,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016D_v1',
        path     = '{0}/RunD'.format(path),
        nJobs    = 10,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016E_v1',
        path     = '{0}/RunE'.format(path),
        nJobs    = 10,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016F_v1',
        path     = '{0}/RunF'.format(path),
        nJobs    = 10,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016G_v1',
        path     = '{0}/RunG'.format(path),
        nJobs    = 19,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016G_v1_nfl',
        path     = '{0}/RunG_nfl'.format(path),
        nJobs    = 1,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016H_v1',
        path     = '{0}/RunH'.format(path),
        nJobs    = 20,
        suffix   = 'electron_2016'
        ),
#   cfg(data_name = 'DoubleEG_Run2016H_v1_nfl',
#       path     = '{0}/RunH_nfl'.format(path),
#       nJobs    = 0,
#       suffix   = 'electron_2016'
#       ),
]


# MONTE CARLO #
mc_dict = {}

path = '/eos/uscms/store/user/jrainbol/Bacon/2016'
mc_dict['signal'] = \
[
    cfg(data_name = 'ZZTo4L',
        path     = '{0}/Summer16_ZZTo4L_powheg'.format(path),
        nJobs    = 7,
        suffix   = 'zz_4l'
        ),
]

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/Summer16_DYJetsToLL_M-50_amcatnlo'.format(path),
        nJobs    = 100,
        suffix   = 'zjets_m-50'
        ),
]

path = '/eos/uscms/store/user/jrainbol/BaconProd'
mc_dict['ttbar'] = \
[
    cfg(data_name = 'TTJets',
        path     = '{0}/ttbar'.format(path),
        nJobs    = 15,
        suffix   = 'ttbar'
        ),
]

mc_dict['diboson'] = \
[
    cfg(data_name = 'WWTo2L2Nu',
        path     = '{0}/ww_2l2nu'.format(path),
        nJobs    = 5,
        suffix   = 'ww_2l2nu'
        ),
    cfg(data_name = 'WZTo2L2Q',
        path     = '{0}/wz_2l2q'.format(path),
        nJobs    = 50,
        suffix   = 'wz_2l2q'
        ),
    cfg(data_name = 'WZTo3LNu',
        path     = '{0}/wz_3lnu'.format(path),
        nJobs    = 30,
        suffix   = 'wz_3lnu'
        ),
    cfg(data_name = 'ZZTo2L2Nu',
        path     = '{0}/zz_2l2nu'.format(path),
        nJobs    = 100,
        suffix   = 'zz_2l2nu'
        ),
    cfg(data_name = 'ZZTo2L2Q',
        path     = '{0}/zz_2l2q'.format(path),
        nJobs    = 5,
        suffix   = 'zz_2l2q'
        ),
]

mc_dict['higgs'] = \
[
    cfg (data_name = 'GluGluHToZZTo4L',
        path     = '{0}/ggH_zz_4l'.format(path),
        nJobs    = 3,
        suffix   = 'ggH_zz_4l'
       ),
    cfg(data_name = 'VBF_HToZZTo4L',
        path     = '{0}/vbfH_zz_4l'.format(path),
        nJobs    = 2,
        suffix   = 'vbfH_zz_4l'
       ),
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