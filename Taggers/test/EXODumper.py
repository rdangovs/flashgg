#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
from FWCore.ParameterSet.VarParsing import VarParsing
from flashgg.MetaData.samples_utils import SamplesManager

process = cms.Process("EXOTagDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("/store/group/phys_higgs/cmshgg/musella/flashgg/EXOMoriond16/1_2_0-136-ge8a0efc/DoubleEG/EXOMoriond16-1_2_0-136-ge8a0efc-v1-Run2015D-16Dec2015-v2/160211_163340/0000/myMicroAODOutputFile_387.root"))

process.TFileService = cms.Service( "TFileService",
                                    fileName = cms.string("EXOTagsDump.root"),
                                    closeFileFast = cms.untracked.bool(True) )

from flashgg.Taggers.flashggTagOutputCommands_cff import tagDefaultOutputCommand
import flashgg.Taggers.dumperConfigTools as cfgTools
from  flashgg.Taggers.tagsDumpers_cfi import createTagDumper

process.exoTagDumper = createTagDumper("EXOTag")
process.exoTagDumper.dumpTrees = True
process.exoTagDumper.dumpHistos = False
process.exoTagDumper.dumpWorkspace = False

dipho_vars = [
        "eventID := exoTag







