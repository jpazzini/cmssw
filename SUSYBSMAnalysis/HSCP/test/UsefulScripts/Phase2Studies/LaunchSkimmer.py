#!/usr/bin/env python

import string, os, sys
import SUSYBSMAnalysis.HSCP.LaunchOnCondor as LaunchOnCondor

datasets     = [
	'/MinBias_140PU_TuneCUETP8M1_14TeV-pythia8/PhaseIIFall16DR82-PU140_90X_upgrade2023_realistic_v1-v1/AODSIM',
	'/MinBias_200PU_TuneCUETP8M1_14TeV-pythia8/PhaseIIFall16DR82-PU200_90X_upgrade2023_realistic_v1-v1/AODSIM'
]

outdir       = 'out'
server       = 'root://cms-xrd-global.cern.ch/'

def outDirName (dataset):
    return dataset.split('/')[1]
   
def getDatasetFiles (dataset):
    return os.popen('das_client --limit=0 --query "file dataset=%s"' % dataset).read().split()

def createOutStructure ():
    os.system ('rm -rf %s && mkdir %s' % (outdir, outdir))
    for dataset in datasets:
        os.system ('mkdir -p %s/%s' % (outdir, outDirName (dataset)))

def initProxy():
      print "You are going to run on a sample over grid using either CRAB or the AAA protocol, it is therefore needed to initialize your grid certificate"
      os.system('mkdir -p ~/x509_user_proxy; voms-proxy-init --voms cms -valid 192:00 --out ~/x509_user_proxy/x509_proxy')#all must be done in the same command to avoid environement problems.  Note that the first sourcing is only needed in Louvain


if sys.argv[1] == '1':
   initProxy ()
   createOutStructure()

   JobName = "dEdxSkimmer"
   FarmDirectory = "FARM"
   LaunchOnCondor.SendCluster_Create(FarmDirectory, JobName)
   LaunchOnCondor.Jobs_Queue = '8nh'
   LaunchOnCondor.Jobs_InitCmds   = ['export HOME=%s' % os.environ['HOME'], 'export X509_USER_PROXY=$HOME/x509_user_proxy/x509_proxy']

   for dataset in datasets:
      datasetMark = outDirName (dataset)
      print '===========================\n%s\n' % datasetMark
      Files = getDatasetFiles(dataset)
      for i in range (0, len(Files)):
         os.system ('cp dEdxSkimmer_Template_cfg.py dEdxSkimmer_cff.py')
         f = open ('dEdxSkimmer_cff.py', 'a')
         f.write ('\n')
         f.write ('process.Out.fileName = cms.untracked.string(\'dEdxSkim_%s_%i.root\')\n' % (datasetMark, i))
         f.write ('process.source.fileNames.extend([\'%s/%s\'])\n' % (server,Files[i]))
         f.close()
         LaunchOnCondor.Jobs_FinalCmds = ['cp dEdxSkim*.root %s/%s/' % (outdir, datasetMark)]
         LaunchOnCondor.SendCluster_Push (["CMSSW", "dEdxSkimmer_cff.py"])
         os.system ('rm -f dEdxSkimmer_cff.py')
#   LaunchOnCondor.SendCluster_Submit ()
