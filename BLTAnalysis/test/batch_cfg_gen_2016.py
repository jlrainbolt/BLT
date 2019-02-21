#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
selection  = 'hard'
period     = '2016'
executable = 'execBatchGen.sh'
location   = 'lpc'

mc_samples  = ['zz_4l', 'zjets']



''' 
Set job configurations.  
'''

# MONTE CARLO #
mc_dict = {}

path = '/eos/uscms/store/user/jrainbol/Bacon/2016'
mc_dict['zz_4l'] = \
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



batch_list = []
batch_list += sum([mc_dict[n] for n in mc_samples], []) 

batch = bm.BatchMaster(config_list = batch_list, 
                      stage_dir   = 'batch',
                      selection  = selection,
                      period     = period,
                      executable = executable,
                      location   = 'lpc'
                     )
batch.submit_to_batch()

