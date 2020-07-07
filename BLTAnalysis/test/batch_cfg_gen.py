#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg         = bm.JobConfig
selection   = 'dressed'
period      = '2012'
executable  = 'execBatchGen.sh'
location    = 'lpc'

mc_samples      = ['signal', 'zjets']



''' 
Set job configurations.  
'''

# MONTE CARLO #
mc_dict = {}

path = '/eos/uscms/store/group/lpcbacon/04'
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
]

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/Summer12_DYJetsToLL_M-50_TuneZ2Star'.format(path),
        nJobs    = 8,
        suffix   = 'zjets_m-50'
        ),
]


batch_list = []
batch_list += sum([mc_dict[n] for n in mc_samples], []) 

batch = bm.BatchMaster(config_list  = batch_list, 
                      stage_dir     = 'batch',
                      selection     = selection,
                      period        = period,
                      executable    = executable,
                      location      = 'lpc'
                     )
batch.submit_to_batch()
