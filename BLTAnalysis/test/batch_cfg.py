#! /usr/bin/env python
import BatchMaster as bm
import sys


''' Specify parameters '''
system     = 'lpc'
cfg        = bm.JobConfig
selection  = 'amumu'
period     = '2016'

if system == 'lpc':
    path       = '/tthome/share/bacon/production'
    executable = 'execBatch.sh'
elif system == 'tier3':
    path       = '/tthome/share/bacon/production'
    executable = 'execBatch.sh'

''' 
    Set job configurations.  The order of arguments is: (Dataset, path to data,
    number of jobs, arguments to pass to executable, output directory name)
'''

dataList = []
dataList.extend([
    cfg(dataName = 'muon_2012A',
        path     = '{0}/04/SingleMu_2012A-22Jan2013'.format(path),
        nJobs    = 1,
        args     = 'muon_2012A muon 2012'
       ),
    cfg(dataName = 'muon_2012B',
        path     = '{0}/04/SingleMu_2012B-22Jan2013'.format(path),
        nJobs    = 15,
        args     = 'muon_2012B muon 2012'
       ),
    cfg(dataName = 'muon_2012C',
        path     = '{0}/04/SingleMu_2012C-22Jan2013'.format(path),
        nJobs    = 15,
        args     = 'muon_2012C muon 2012'
       ),
    cfg(dataName = 'muon_2012D',
        path     = '{0}/04/SingleMu_2012D-22Jan2013'.format(path),
        nJobs    = 15,
        args     = 'muon_2012D muon 2012'
       )
    ])

batch = bm.BatchMaster(configList = dataList, 
                      shortQueue = False,
                      stageDir   = 'batch',
                      executable = executable,
                      selection  = '{0}_{1}'.format(selection, period)
                     )
batch.submit_to_batch()

