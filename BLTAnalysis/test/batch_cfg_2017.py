#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg         = bm.JobConfig
selection   = 'loose'
period      = '2017'
executable  = 'execBatch.sh'
location    = 'lpc'

data_samples    = ['double_mu', 'double_eg']
mc_samples      = ['zjets', 'ttbar', 'diboson', 'higgs']



''' 
Set job configurations.  
'''

# DATA #
data_dict = {}

path = '/eos/uscms/store/user/jbueghly/jbueghly_2017_data/'
data_dict['double_mu'] = \
[
    cfg(data_name = 'muon_2017B_v1',
        path     = '{0}/DoubleMuon_Run2017B-31Mar2018-v1'.format(path),
        nJobs    = 52,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017C_v1',
        path     = '{0}/DoubleMuon_Run2017C-31Mar2018-v1'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017D_v1',
        path     = '{0}/DoubleMuon_Run2017D-31Mar2018-v1'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017E_v1',
        path     = '{0}/DoubleMuon_Run2017E-31Mar2018-v1'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017F_v1',
        path     = '{0}/DoubleMuon_Run2017F-31Mar2018-v1'.format(path),
        nJobs    = 51,
        suffix   = 'muon_2017'
        ),
]

data_dict['double_eg'] = \
[
    cfg(data_name = 'electron_2017B_v1',
        path     = '{0}/DoubleEG_Run2017B-31Mar2018-v1'.format(path),
        nJobs    = 51,
        suffix   = 'electron_2017'
        ),
    cfg(data_name = 'electron_2017C_v1',
        path     = '{0}/DoubleEG_Run2017C-31Mar2018-v1'.format(path),
        nJobs    = 51,
        suffix   = 'electron_2017'
        ),
    cfg(data_name = 'electron_2017D_v1',
        path     = '{0}/DoubleEG_Run2017D-31Mar2018-v1'.format(path),
        nJobs    = 51,
        suffix   = 'electron_2017'
        ),
    cfg(data_name = 'electron_2017E_v1',
        path     = '{0}/DoubleEG_Run2017E-31Mar2018-v1'.format(path),
        nJobs    = 51,
        suffix   = 'electron_2017'
        ),
    cfg(data_name = 'electron_2017F_v1',
        path     = '{0}/DoubleEG_Run2017F-31Mar2018-v1'.format(path),
        nJobs    = 51,
        suffix   = 'electron_2017'
        ),
]


# MONTE CARLO #
mc_dict = {}

path = '/eos/uscms/store/group/lpcbacon/2017/'
mc_dict['signal'] = \
[
    cfg(data_name = 'ZZTo4L',
        path     = '{0}/Fall17_ZZTo4L_powheg'.format(path),
        nJobs    = 20,
        suffix   = 'zz_4l'
        ),
]

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/Fall17_DYJetsToLL_M-50_amcatnlo'.format(path),
        nJobs    = 50,
        suffix   = 'zjets_m-50'
        ),
]

mc_dict['ttbar'] = \
[
    cfg(data_name = 'TTJets',
        path     = '{0}/Fall17_TTJets_amcatnlo'.format(path),
        nJobs    = 50,
        suffix   = 'ttbar'
        ),
]

mc_dict['diboson'] = \
[
    cfg(data_name = 'WWTo2L2Nu',
        path     = '{0}/Fall17_WWTo2L2Nu_powheg'.format(path),
        nJobs    = 5,
        suffix   = 'ww_2l2nu'
        ),
    cfg(data_name = 'WZTo2L2Q',
        path     = '{0}/Fall17_WZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 50,
        suffix   = 'wz_2l2q'
        ),
    cfg(data_name = 'WZTo3LNu',
        path     = '{0}/Fall17_WZTo3LNu_amcatnlo'.format(path),
        nJobs    = 30,
        suffix   = 'wz_3lnu'
        ),
#   cfg(data_name = 'ZZTo2L2Nu',
#       path     = '{0}/Fall17_ZZTo2L2Nu_powheg'.format(path),
#       nJobs    = 10,
#       suffix   = 'zz_2l2nu'
#       ),
    cfg(data_name = 'ZZTo2L2Q',
        path     = '{0}/Fall17_ZZTo2L2Q_amcatnlo'.format(path),
        nJobs    = 50,
        suffix   = 'zz_2l2q'
        ),
    cfg(data_name = 'ZZTo4L',
        path     = '{0}/Fall17_ZZTo4L_powheg'.format(path),
        nJobs    = 20,
        suffix   = 'zz_4l'
        ),
]

mc_dict['higgs'] = \
[
    cfg (data_name = 'GluGluHToZZTo4L',
        path     = '{0}/Fall17_GluGluHToZZTo4L_powheg'.format(path),
        nJobs    = 3,
        suffix   = 'ggH_zz_4l'
       ),
    cfg(data_name = 'VBF_HToZZTo4L',
        path     = '{0}/Fall17_VBF_HToZZTo4L_powheg'.format(path),
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

