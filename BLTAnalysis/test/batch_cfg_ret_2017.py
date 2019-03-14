#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
selection  = '4l'
period     = '2017'
executable = 'execBatchRetriever.sh'
location   = 'lpc'

mc_samples  = ['zjets']



''' 
Set job configurations.  
'''

# MONTE CARLO #
mc_dict = {}

path = '/eos/uscms/store/group/lpcbacon/jlr'

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/DYJetsToLL_M-50_amcatnlo_RunIIFall17_new_pmx_v14_ext1-v1'.format(path),
        nJobs    = 175,
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

