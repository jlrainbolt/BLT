#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg         = bm.JobConfig
selection   = 'tight'
period      = '2018'
executable  = 'execBatch.sh'
location    = 'lpc'

#data_samples    = ['double_mu', 'single_mu', 'egamma']
data_samples    = []
#mc_samples      = ['signal', 'zjets', 'ttbar', 'diboson', 'higgs', 'triboson']
mc_samples      = ['zjets']



''' 
Set job configurations.  
'''

path = '/eos/uscms/store/group/lpcbacon/jlr/'

# DATA #
data_dict = {}

data_dict['double_mu'] = \
[
    cfg(data_name = 'DoubleMuon_Run2018A',
        path     = '{0}/DoubleMuon_Run2018A-17Sep2018-v2'.format(path),
        nJobs    = 40,
        suffix   = 'muon_2018'
        ),
    cfg(data_name = 'DoubleMuon_Run2018B',
        path     = '{0}/DoubleMuon_Run2018B-17Sep2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2018'
        ),
    cfg(data_name = 'DoubleMuon_Run2018C',
        path     = '{0}/DoubleMuon_Run2018C-17Sep2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2018'
        ),
    cfg(data_name = 'DoubleMuon_Run2018D',
        path     = '{0}/DoubleMuon_Run2018D-PromptReco-v2'.format(path),
        nJobs    = 70,
        suffix   = 'muon_2018'
        ),
]

data_dict['single_mu'] = \
[
    cfg(data_name = 'SingleMuon_Run2018A',
        path     = '{0}/SingleMuon_Run2018A-17Sep2018-v2'.format(path),
        nJobs    = 190,
        suffix   = 'muon_2018'
        ),
    cfg(data_name = 'SingleMuon_Run2018B',
        path     = '{0}/SingleMuon_Run2018B-17Sep2018-v1'.format(path),
        nJobs    = 90,
        suffix   = 'muon_2018'
        ),
    cfg(data_name = 'SingleMuon_Run2018C',
        path     = '{0}/SingleMuon_Run2018C-17Sep2018-v1'.format(path),
        nJobs    = 85,
        suffix   = 'muon_2018'
        ),
    cfg(data_name = 'SingleMuon_Run2018D',
        path     = '{0}/SingleMuon_Run2018D-PromptReco-v2'.format(path),
        nJobs    = 385,
        suffix   = 'muon_2018'
        ),
]

data_dict['egamma'] = \
[
    cfg(data_name = 'EGamma_Run2018A',
        path     = '{0}/EGamma_Run2018A-17Sep2018-v2'.format(path),
        nJobs    = 190,
        suffix   = 'electron_2018'
        ),
    cfg(data_name = 'EGamma_Run2018B',
        path     = '{0}/EGamma_Run2018B-17Sep2018-v1'.format(path),
        nJobs    = 80,
        suffix   = 'electron_2018'
        ),
    cfg(data_name = 'EGamma_Run2018C',
        path     = '{0}/EGamma_Run2018C-17Sep2018-v1'.format(path),
        nJobs    = 75,
        suffix   = 'electron_2018'
        ),
    cfg(data_name = 'EGamma_Run2018D',
        path     = '{0}/EGamma_Run2018D-PromptReco-v2'.format(path),
        nJobs    = 330,
        suffix   = 'electron_2018'
        ),
]


# MONTE CARLO #
mc_dict = {}

mc_dict['signal'] = \
[
#   cfg(data_name = 'ZZTo4L',
#       path     = '{0}/ZZTo4L_powheg_RunIIAutumn18_v15_ext1-v2'.format(path),
#       nJobs    = 15,
#       suffix   = 'zz_4l'
#       ),
#   cfg(data_name = 'ZZTo4L_aMC',
#       path     = '{0}/ZZTo4L_amcatnlo_RunIIAutumn18_v15-v1'.format(path),
#       nJobs    = 65,
#       suffix   = 'zz_4l_aMC'
#       ),
    cfg(data_name = 'ZZTo4L_M-1',
        path     = '{0}/ZZTo4L_M-1_powheg_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 125,
        suffix   = 'zz_4l_m-1'
        ),
]

mc_dict['zjets'] = \
[
    cfg(data_name = 'DYJetsToLL_M-50',
        path     = '{0}/DYJetsToLL_M-50_madgraph_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 200,
        suffix   = 'zjets_m-50'
        ),
]

mc_dict['ttbar'] = \
[
    cfg(data_name = 'TTJets',
        path     = '{0}/TTJets_amcatnlo_RunIIAutumn18_v15_ext1-v2'.format(path),
        nJobs    = 200,
        suffix   = 'ttbar'
        ),
    cfg(data_name = 'TTTo2L2Nu',
        path     = '{0}/TTTo2L2Nu_powheg_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 100,
        suffix   = 'tt_2l2nu'
        ),
    cfg(data_name = 'TTZTo2L2Nu',
        path     = '{0}/TTZToLLNuNu_powheg_RunIIAutumn18_v15_ext1-v2'.format(path),
        nJobs    = 25,
        suffix   = 'ttz_2l2nu'
        ),
]

mc_dict['diboson'] = \
[
    cfg(data_name = 'WWTo2L2Nu',
        path     = '{0}/WWTo2L2Nu_powheg_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 15,
        suffix   = 'ww_2l2nu'
        ),
    cfg(data_name = 'WZTo2L2Q',
        path     = '{0}/WZTo2L2Q_amcatnlo_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 70,
        suffix   = 'wz_2l2q'
        ),
    cfg(data_name = 'WZTo3LNu',
        path     = '{0}/WZTo3LNu_amcatnlo_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 20,
        suffix   = 'wz_3lnu'
        ),
    cfg(data_name = 'ZZTo2L2Nu',
        path     = '{0}/ZZTo2L2Nu_powheg_RunIIAutumn18_v15_ext1-v2'.format(path),
        nJobs    = 15,
        suffix   = 'zz_2l2nu'
        ),
    cfg(data_name = 'ZZTo2L2Q',
        path     = '{0}/ZZTo2L2Q_amcatnlo_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 70,
        suffix   = 'zz_2l2q'
        ),
]

mc_dict['higgs'] = \
[
    cfg (data_name = 'GluGluHToZZTo4L',
        path     = '{0}/GluGluHToZZTo4L_M125_powheg_RunIIAutumn18_v15-v2'.format(path),
        nJobs    = 2,
        suffix   = 'ggH_zz_4l'
       ),
    cfg(data_name = 'VBF_HToZZTo4L',
        path     = '{0}/VBF_HToZZTo4L_M125_powheg_RunIIAutumn18_v15-v2'.format(path),
        nJobs    = 2,
        suffix   = 'vbfH_zz_4l'
       ),
]

mc_dict['triboson'] = \
[
    cfg(data_name = 'WWZJetsTo4L2Nu',
        path     = '{0}/WWZJetsTo4L2Nu_amcatnlo_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 5,
        suffix   = 'wwz_4l2nu'
        ),
    cfg(data_name = 'WZZJetsTo4L2Nu',
        path     = '{0}/WZZJetsTo4L2Nu_amcatnlo_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 5,
        suffix   = 'wzz_4l2nu'
        ),
    cfg(data_name = 'ZZZJetsTo4L2Nu',
        path     = '{0}/ZZZJetsTo4L2Nu_amcatnlo_RunIIAutumn18_v15-v1'.format(path),
        nJobs    = 5,
        suffix   = 'zzz_4l2nu'
        ),
    cfg(data_name = 'ZZGJetsTo4L2Nu',
        path     = '{0}/ZZGJetsTo4L2Nu_amcatnlo_RunIIAutumn18_v15-v1'.format(path),
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

