import sys, os
import FWCore.ParameterSet.Config as cms

isSignal = False
isBckg = False
isData = True
isSkimmedSample = False
GTAG = '80X_dataRun2_Prompt_v8'
OUTPUTFILE = 'HSCP.root'
LUMITOPROCESS = ''

#debug input files 
#this list is overwritten by CRAB
InputFileList = cms.untracked.vstring(
   'root://eoscms//eos/cms/store/data/Run2015B/DoubleMuon/RECO/PromptReco-v1/000/251/521/00000/5E8B4110-6429-E511-A64C-02163E0127DF.root',
)

#main EDM tuple cfg that depends on the above parameters
execfile( os.path.expandvars('${CMSSW_BASE}/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/HSCParticleProducer_cfg.py') )
