#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("Bfinder")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

## Source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
    'file:///home/ltsai/ReceivedFile/DATA/8028_2016RunG_AOD_07Aug17.root'
    )
        )


from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring('HLT_Dimuon16_Jpsi_v*'))

## Geometry and Detector Conditions (needed for a few patTuple production steps)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016LegacyRepro_v4', '')

# preselect pat muon and pat tracks
process.load("BPHAnalysis.PreselectFilter.FilterConf9_cfi")
# remove MC dependence
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, names=['All'], outputModules=[] )

### ExcitedBc analyizer
process.analysis = cms.EDAnalyzer('Bfinder',
        Bchannel                = cms.vint32(
            1,#RECONSTRUCTION: J/psi + p+  K-
            1,#RECONSTRUCTION: J/psi + p-  K+
            1,#RECONSTRUCTION: J/psi + lambda (p+, pi-) 
            1,),#RECONSTRUCTION: J/psi + lambda (p-, pi+) 
        MuonLabel   = cms.string('selectedMuons'),         #selectedPatMuons
        TrackLabel  = cms.string('selectedTracks'),    #selectedPat
        BSLabel     = cms.string("offlineBeamSpot::RECO"),
        PVLabel     = cms.string("offlinePrimaryVertices::RECO"),
        GenLabel    = cms.string('genParticles'),
        tkPtCut = cms.double(0.3),
        jpsiPtCut = cms.double(16.0),
        bPtCut = cms.double(16.0),
        RunOnMC = cms.bool(False),
        doTkPreCut = cms.bool(True),
        doMuPreCut = cms.bool(True)
        )

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('bfinder.root')
        )


process.p = cms.Path(
      process.hltHighLevel
    * process.myMuonsSequence
    * process.myTrackSequence
    * process.analysis
)
