# ================================================
#   MicroAOD and Jet Tree Producer
#   Y. Haddad 01/2015
# ================================================

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from flashgg.MicroAOD.flashggJets_cfi import flashggBTag, maxJetCollections
process = cms.Process("FLASHggEXOValidation")
process.load("FWCore.MessageService.MessageLogger_cfi")

# +++++ the number of processed events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 100 ) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 10 )

# +++++ the source file
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring("/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_111.root")) 

#process.MessageLogger.cerr.threshold = 'ERROR'

process.TFileService = cms.Service("TFileService",fileName  = cms.string("exoValidationTrees_test.root"))

# process.flashggEXOValidationTreeMaker = cms.EDAnalyzer('FlashggEXOValidationTreeMaker',
#                                                        GenParticleTag           = cms.untracked.InputTag('prunedGenParticles'),
#                                                        JetTagDz                 = cms.InputTag("flashggJets"),
#                                                        StringTag		= cms.string("PF"),
#                                                    )
#

#JetCollectionVInputTagPFCHS = cms.VInputTag()
#JetCollectionVInputTagPUPPI = cms.VInputTag()
#for i in range(0,5):
#    JetCollectionVInputTagPFCHS.append(cms.InputTag('flashggPFCHSJets' + str(i)))
    #JetCollectionVInputTagPUPPI.append(cms.InputTag('flashggPUPPIJets' + str(i)))

    
process.flashggEXOValidationTreeMaker = cms.EDAnalyzer('FlashggEXOValidationTreeMaker',
                                                            DiPhotonTag     = cms.InputTag('flashggDiPhotons'),
                                                            )
#

process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('testetstestst.root'),
                               outputCommands = cms.untracked.vstring())

process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.p = cms.Path( 
    process.flashggEXOValidationTreeMaker 
    )
process.e = cms.EndPath(process.out)

# customization for job splitting, lumi weighting, etc.
from diphotons.MetaData.JobConfig import customize
customize.setDefault("maxEvents",-1)
customize(process)

