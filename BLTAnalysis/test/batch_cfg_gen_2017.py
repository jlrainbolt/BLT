#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
selection  = 'hard'
period     = '2017'
executable = 'execBatchGen.sh'
location   = 'lpc'

#mc_samples  = ['zz_4l', 'zjets']
mc_samples  = ['zz_4l']



''' 
Set job configurations.  
'''

# MONTE CARLO #
mc_dict = {}

path = '/eos/uscms/store/group/lpcbacon/jlr'
mc_dict['zz_4l'] = \
[
    cfg(data_name = 'ZZTo4L',
        path     = '{0}/ZZTo4L_powheg_RunIIFall17_new_pmx_v14-v1'.format(path),
        nJobs    = 15,
        suffix   = 'zz_4l'
        ),
]

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/DYJetsToLL_M-50_amcatnlo_RunIIFall17_new_pmx_v14_ext1-v1'.format(path),
        nJobs    = 350,
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

