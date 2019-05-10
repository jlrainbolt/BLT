#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg         = bm.JobConfig
selection   = 'mumuIso'
period      = '2017'
executable  = 'execBatchTnP.sh'
location    = 'lpc'

#data_samples    = ['single_mu', 'single_eg']
data_samples    = ['single_mu']
mc_samples      = ['zjets', 'ttbar']
#mc_samples      = ['ttbar']



''' 
Set job configurations.  
'''

# DATA #
data_dict = {}

path = '/eos/uscms/store/group/lpcbacon/jlr'
data_dict['single_mu'] = \
[
    cfg(data_name = 'muon_2017B_v1',
        path     = '{0}/SingleMuon_Run2017B-31Mar2018-v1'.format(path),
        nJobs    = 60,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017C_v1',
        path     = '{0}/SingleMuon_Run2017C-31Mar2018-v1'.format(path),
        nJobs    = 120,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017D_v1',
        path     = '{0}/SingleMuon_Run2017D-31Mar2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017E_v1',
        path     = '{0}/SingleMuon_Run2017E-31Mar2018-v1'.format(path),
        nJobs    = 130,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017F_v1',
        path     = '{0}/SingleMuon_Run2017F-31Mar2018-v1'.format(path),
        nJobs    = 200,
        suffix   = 'muon_2017'
        ),
]

data_dict['double_eg'] = \
[
    cfg(data_name = 'electron_2017B_v1',
        path     = '{0}/DoubleEG_Run2017B-31Mar2018-v1'.format(path),
        nJobs    = 15,
        suffix   = 'electron_2017'
        ),
    cfg(data_name = 'electron_2017C_v1',
        path     = '{0}/DoubleEG_Run2017C-31Mar2018-v1'.format(path),
        nJobs    = 30,
        suffix   = 'electron_2017'
        ),
    cfg(data_name = 'electron_2017D_v1',
        path     = '{0}/DoubleEG_Run2017D-31Mar2018-v1'.format(path),
        nJobs    = 5,
        suffix   = 'electron_2017'
        ),
    cfg(data_name = 'electron_2017E_v1',
        path     = '{0}/DoubleEG_Run2017E-31Mar2018-v1'.format(path),
        nJobs    = 20,
        suffix   = 'electron_2017'
        ),
    cfg(data_name = 'electron_2017F_v1',
        path     = '{0}/DoubleEG_Run2017F-31Mar2018-v1'.format(path),
        nJobs    = 25,
        suffix   = 'electron_2017'
        ),
]


# MONTE CARLO #
mc_dict = {}

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/DYJetsToLL_M-50_amcatnlo_RunIIFall17_new_pmx_v14_ext1-v1'.format(path),
        nJobs    = 350,
        suffix   = 'zjets_m-50'
        ),
]

mc_dict['ttbar'] = \
[
    cfg(data_name = 'TTJets',
        path     = '{0}/TTJets_amcatnlo_RunIIFall17_new_pmx_v14-v1'.format(path),
        nJobs    = 200,
        suffix   = 'ttbar'
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

