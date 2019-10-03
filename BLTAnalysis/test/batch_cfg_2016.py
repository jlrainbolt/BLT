#! /usr/bin/env python
import BLT.BLTAnalysis.BatchMaster as bm
import sys


''' Specify parameters '''
cfg         = bm.JobConfig
selection   = 'tight'
period      = '2016'
executable  = 'execBatch.sh'
location    = 'lpc'

#data_samples    = ['double_mu', 'double_eg', 'single_mu', 'single_e']
data_samples    = []
#mc_samples      = ['signal', 'zjets', 'ttbar', 'diboson', 'higgs']
mc_samples      = ['zjets']



''' 
Set job configurations.  
'''

path = '/eos/uscms/store/group/lpcbacon/jlr'

# DATA #
data_dict = {}

data_dict['double_mu'] = \
[
    cfg(data_name = 'DoubleMuon_Run2016B_v1',
        path     = '{0}/DoubleMuon_Run2016B-17Jul2018_ver1-v1'.format(path),
        nJobs    = 2,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016B_v2',
        path     = '{0}/DoubleMuon_Run2016B-17Jul2018_ver2-v1'.format(path),
        nJobs    = 25,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016C_v1',
        path     = '{0}/DoubleMuon_Run2016C-17Jul2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016D_v1',
        path     = '{0}/DoubleMuon_Run2016D-17Jul2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016E_v1',
        path     = '{0}/DoubleMuon_Run2016E-17Jul2018-v1'.format(path),
        nJobs    = 15,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016F_v1',
        path     = '{0}/DoubleMuon_Run2016F-17Jul2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016G_v1',
        path     = '{0}/DoubleMuon_Run2016G-17Jul2018-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'DoubleMuon_Run2016H_v1',
        path     = '{0}/DoubleMuon_Run2016H-17Jul2018-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016'
        ),
]

data_dict['single_mu'] = \
[
    cfg(data_name = 'SingleMuon_Run2016B_v1',
        path     = '{0}/SingleMuon_Run2016B-17Jul2018_ver1-v1'.format(path),
        nJobs    = 1,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'SingleMuon_Run2016B_v2',
        path     = '{0}/SingleMuon_Run2016B-17Jul2018_ver2-v1'.format(path),
        nJobs    = 50,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'SingleMuon_Run2016C_v1',
        path     = '{0}/SingleMuon_Run2016C-17Jul2018-v1'.format(path),
        nJobs    = 20,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'SingleMuon_Run2016D_v1',
        path     = '{0}/SingleMuon_Run2016D-17Jul2018-v1'.format(path),
        nJobs    = 40,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'SingleMuon_Run2016E_v1',
        path     = '{0}/SingleMuon_Run2016E-17Jul2018-v1'.format(path),
        nJobs    = 40,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'SingleMuon_Run2016F_v1',
        path     = '{0}/SingleMuon_Run2016F-17Jul2018-v1'.format(path),
        nJobs    = 30,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'SingleMuon_Run2016G_v1',
        path     = '{0}/SingleMuon_Run2016G-17Jul2018-v1'.format(path),
        nJobs    = 80,
        suffix   = 'muon_2016'
        ),
    cfg(data_name = 'SingleMuon_Run2016H_v1',
        path     = '{0}/SingleMuon_Run2016H-17Jul2018-v1'.format(path),
        nJobs    = 100,
        suffix   = 'muon_2016'
        ),
]

data_dict['double_eg'] = \
[
    cfg(data_name = 'DoubleEG_Run2016B_v1',
        path     = '{0}/DoubleEG_Run2016B-17Jul2018_ver1-v1'.format(path),
        nJobs    = 5,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016B_v2',
        path     = '{0}/DoubleEG_Run2016B-17Jul2018_ver2-v1'.format(path),
        nJobs    = 70,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016C_v1',
        path     = '{0}/DoubleEG_Run2016C-17Jul2018-v1'.format(path),
        nJobs    = 20,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016D_v1',
        path     = '{0}/DoubleEG_Run2016D-17Jul2018-v1'.format(path),
        nJobs    = 15,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016E_v1',
        path     = '{0}/DoubleEG_Run2016E-17Jul2018-v1'.format(path),
        nJobs    = 15,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016F_v1',
        path     = '{0}/DoubleEG_Run2016F-17Jul2018-v1'.format(path),
        nJobs    = 10,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016G_v1',
        path     = '{0}/DoubleEG_Run2016G-17Jul2018-v1'.format(path),
        nJobs    = 25,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'DoubleEG_Run2016H_v1',
        path     = '{0}/DoubleEG_Run2016H-17Jul2018-v1'.format(path),
        nJobs    = 25,
        suffix   = 'electron_2016'
        ),
]

data_dict['single_e'] = \
[
    cfg(data_name = 'SingleElectron_Run2016B_v1',
        path     = '{0}/SingleElectron_Run2016B-17Jul2018_ver1-v1'.format(path),
        nJobs    = 1,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'SingleElectron_Run2016B_v2',
        path     = '{0}/SingleElectron_Run2016B-17Jul2018_ver2-v1'.format(path),
        nJobs    = 140,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'SingleElectron_Run2016C_v1',
        path     = '{0}/SingleElectron_Run2016C-17Jul2018-v1'.format(path),
        nJobs    = 60,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'SingleElectron_Run2016D_v1',
        path     = '{0}/SingleElectron_Run2016D-17Jul2018-v1'.format(path),
        nJobs    = 90,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'SingleElectron_Run2016E_v1',
        path     = '{0}/SingleElectron_Run2016E-17Jul2018-v1'.format(path),
        nJobs    = 70,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'SingleElectron_Run2016F_v1',
        path     = '{0}/SingleElectron_Run2016F-17Jul2018-v1'.format(path),
        nJobs    = 50,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'SingleElectron_Run2016G_v1',
        path     = '{0}/SingleElectron_Run2016G-17Jul2018-v1'.format(path),
        nJobs    = 110,
        suffix   = 'electron_2016'
        ),
    cfg(data_name = 'SingleElectron_Run2016H_v1',
        path     = '{0}/SingleElectron_Run2016H-17Jul2018-v1'.format(path),
        nJobs    = 90,
        suffix   = 'electron_2016'
        ),
]



# MONTE CARLO #
mc_dict = {}

mc_dict['signal'] = \
[
#   cfg(data_name = 'ZZTo4L',
#       path     = '{0}/ZZTo4L_powheg_RunIISummer16_v3-v1'.format(path),
#       nJobs    = 15,
#       suffix   = 'zz_4l'
#       ),
#   cfg(data_name = 'ZZTo4L_aMC',
#       path     = '{0}/ZZTo4L_amcatnlo_RunIISummer16_v3_ext1-v2'.format(path),
#       nJobs    = 20,
#       suffix   = 'zz_4l_aMC'
#       ),
    cfg(data_name = 'ZZTo4L_M-1',
        path     = '{0}/ZZTo4L_M-1_powheg_RunIISummer16_v3-v1'.format(path),
        nJobs    = 125,
        suffix   = 'zz_4l_m-1'
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

mc_dict['ttbar'] = \
[
    cfg(data_name = 'TTJets',
        path     = '{0}/TTJets_amcatnlo_RunIISummer16_v3-v2'.format(path),
        nJobs    = 85,
        suffix   = 'ttbar'
        ),
    cfg(data_name = 'TTTo2L2Nu',
        path     = '{0}/TTTo2L2Nu_ttHtranche3_powheg_RunIISummer16_v3-v2'.format(path),
        nJobs    = 140,
        suffix   = 'tt_2l2nu'
        ),
    cfg(data_name = 'TTZTo2L2Nu',
        path     = '{0}/TTZToLLNuNu_amcatnlo_RunIISummer16_v3_ext3-v1'.format(path),
        nJobs    = 20,
        suffix   = 'ttz_2l2nu'
        ),
]

mc_dict['diboson'] = \
[
    cfg(data_name = 'WWTo2L2Nu',
        path     = '{0}/WWTo2L2Nu_powheg_RunIISummer16_v3-v2'.format(path),
        nJobs    = 4,
        suffix   = 'ww_2l2nu'
        ),
    cfg(data_name = 'WZTo2L2Q',
        path     = '{0}/WZTo2L2Q_amcatnlo_RunIISummer16_v3-v2'.format(path),
        nJobs    = 65,
        suffix   = 'wz_2l2q'
        ),
    cfg(data_name = 'WZTo3LNu',
        path     = '{0}/WZTo3LNu_amcatnlo_RunIISummer16_v3-v2'.format(path),
        nJobs    = 25,
        suffix   = 'wz_3lnu'
        ),
    cfg(data_name = 'ZZTo2L2Nu',
        path     = '{0}/ZZTo2L2Nu_powheg_RunIISummer16_v3_ext1-v2'.format(path),
        nJobs    = 100,
        suffix   = 'zz_2l2nu'
        ),
    cfg(data_name = 'ZZTo2L2Q',
        path     = '{0}/ZZTo2L2Q_powheg_RunIISummer16_v3-v2'.format(path),
        nJobs    = 1,
        suffix   = 'zz_2l2q'
        ),
]

mc_dict['higgs'] = \
[
    cfg (data_name = 'GluGluHToZZTo4L',
        path     = '{0}/GluGluHToZZTo4L_M125_powheg_RunIISummer16_v3-v1'.format(path),
        nJobs    = 2,
        suffix   = 'ggH_zz_4l'
       ),
    cfg(data_name = 'VBF_HToZZTo4L',
        path     = '{0}/VBF_HToZZTo4L_M125_powheg_RunIISummer16_v3_ext1-v1'.format(path),
        nJobs    = 1,
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
