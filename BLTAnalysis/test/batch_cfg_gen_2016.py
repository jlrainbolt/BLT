#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg        = bm.JobConfig
selection  = 'dressed'
period     = '2016'
executable = 'execBatchGen.sh'
location   = 'lpc'

mc_samples  = ['zz_4l', 'zjets']



''' 
Set job configurations.  
'''

path = '/eos/uscms/store/group/lpcbacon/jlr'

# MONTE CARLO #
mc_dict = {}

mc_dict['zz_4l'] = \
[
    cfg(data_name = 'ZZTo4L',
        path     = '{0}/ZZTo4L_powheg_RunIISummer16_v3-v1'.format(path),
        nJobs    = 15,
        suffix   = 'zz_4l'
        ),
]

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/DYJetsToLL_M-50_amcatnlo_RunIISummer16_v3_ext2-v1'.format(path),
        nJobs    = 240,
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

