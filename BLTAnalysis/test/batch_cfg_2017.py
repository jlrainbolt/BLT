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
#mc_samples      = ['signal', 'zjets', 'ttbar', 'diboson', 'higgs', 'triboson']
#data_samples    = []
mc_samples      = ['zjets', 'signal']



''' 
Set job configurations.  
'''

# DATA #
data_dict = {}

path = '/eos/uscms/store/group/lpcbacon/jlr/'
data_dict['double_mu'] = \
[
    cfg(data_name = 'muon_2017B_v1',
        path     = '{0}/DoubleMuon_Run2017B-31Mar2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017C_v1',
        path     = '{0}/DoubleMuon_Run2017C-31Mar2018-v1'.format(path),
        nJobs    = 30,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017D_v1',
        path     = '{0}/DoubleMuon_Run2017D-31Mar2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017E_v1',
        path     = '{0}/DoubleMuon_Run2017E-31Mar2018-v1'.format(path),
        nJobs    = 30,
        suffix   = 'muon_2017'
        ),
    cfg(data_name = 'muon_2017F_v1',
        path     = '{0}/DoubleMuon_Run2017F-31Mar2018-v1'.format(path),
        nJobs    = 40,
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

mc_dict['signal'] = \
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
#   cfg(data_name = 'DYJetsToLL_M-10to50',
#       path     = '{0}/DYJetsToLL_M-10to50_madgraph_RunIIFall17_v14_ext1-v2'.format(path),
#       nJobs    = 65,
#       suffix   = 'zjets_m-10'
#       ),
]

mc_dict['ttbar'] = \
[
    cfg(data_name = 'TTJets',
        path     = '{0}/TTJets_amcatnlo_RunIIFall17_new_pmx_v14-v1'.format(path),
        nJobs    = 200,
        suffix   = 'ttbar'
        ),
    cfg(data_name = 'TTTo2L2Nu',
        path     = '{0}/TTTo2L2Nu_powheg_RunIIFall17_new_pmx_v14-v1'.format(path),
        nJobs    = 30,
        suffix   = 'tt_2l2nu'
        ),
    cfg(data_name = 'TTZTo2L2Nu',
        path     = '{0}/TTZToLLNuNu_powheg_RunIIFall17_v14-v1'.format(path),
        nJobs    = 30,
        suffix   = 'ttz_2l2nu'
        ),
]

mc_dict['diboson'] = \
[
    cfg(data_name = 'WWTo2L2Nu',
        path     = '{0}/WWTo2L2Nu_powheg_RunIIFall17_v14-v1'.format(path),
        nJobs    = 5,
        suffix   = 'ww_2l2nu'
        ),
    cfg(data_name = 'WZTo2L2Q',
        path     = '{0}/WZTo2L2Q_amcatnlo_RunIIFall17_v14-v1'.format(path),
        nJobs    = 70,
        suffix   = 'wz_2l2q'
        ),
    cfg(data_name = 'WZTo3LNu',
        path     = '{0}/WZTo3LNu_amcatnlo_RunIIFall17_new_pmx_v14-v1'.format(path),
        nJobs    = 25,
        suffix   = 'wz_3lnu'
        ),
    cfg(data_name = 'ZZTo2L2Nu',
        path     = '{0}/ZZTo2L2Nu_powheg_RunIIFall17_v14-v1'.format(path),
        nJobs    = 20,
        suffix   = 'zz_2l2nu'
        ),
    cfg(data_name = 'ZZTo2L2Q',
        path     = '{0}/ZZTo2L2Q_amcatnlo_RunIIFall17_v14-v1'.format(path),
        nJobs    = 75,
        suffix   = 'zz_2l2q'
        ),
]

mc_dict['higgs'] = \
[
    cfg (data_name = 'GluGluHToZZTo4L',
        path     = '{0}/GluGluHToZZTo4L_M125_powheg_RunIIFall17_v14_ext1-v1'.format(path),
        nJobs    = 2,
        suffix   = 'ggH_zz_4l'
       ),
    cfg(data_name = 'VBF_HToZZTo4L',
        path     = '{0}/VBF_HToZZTo4L_M125_powheg_RunIIFall17_v14_ext1-v1'.format(path),
        nJobs    = 2,
        suffix   = 'vbfH_zz_4l'
       ),
]

mc_dict['triboson'] = \
[
    cfg(data_name = 'WWZJetsTo4L2Nu',
        path     = '{0}/WWZJetsTo4L2Nu_amcatnlo_RunIIFall17_v14-v1'.format(path),
        nJobs    = 5,
        suffix   = 'wwz_4l2nu'
        ),
    cfg(data_name = 'WZZJetsTo4L2Nu',
        path     = '{0}/WZZJetsTo4L2Nu_amcatnlo_RunIIFall17_v14-v1'.format(path),
        nJobs    = 5,
        suffix   = 'wzz_4l2nu'
        ),
    cfg(data_name = 'ZZZJetsTo4L2Nu',
        path     = '{0}/ZZZJetsTo4L2Nu_amcatnlo_RunIIFall17_v14-v1'.format(path),
        nJobs    = 5,
        suffix   = 'zzz_4l2nu'
        ),
    cfg(data_name = 'ZZGJetsTo4L2Nu',
        path     = '{0}/ZZGJetsTo4L2Nu_amcatnlo_RunIIFall17_v14-v1'.format(path),
        nJobs    = 5,
        suffix   = 'zzg_4l2nu'
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

