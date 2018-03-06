import FWCore.ParameterSet.Config as cms

process = cms.Process("bphAnalysis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
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


process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
#'file:///home/ltsai/Data/CMSSW_8_0_21/LbToJPsiLam0/08A76328-3362-E711-90FD-0025905A611C.root'
'file:///home/ltsai/Data/8028_2016RunG_AOD_07Aug17.root'
))
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016LegacyRepro_v4', '')

HLTName='HLT_DoubleMu4_JpsiTrk_Displaced_v*'
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(HLTName))

#process.load("BPHAnalysis.SpecificDecay.testProducer_cfi")
process.load("BPHAnalysis.PreselectFilter.FilterConf_cfi")

# remove MC dependence
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, names=['All'], outputModules=[] )

from BPHAnalysis.SpecificDecay.newSelectForWrite_cfi import recoSelect
process.lbWriteSpecificDecay = cms.EDProducer('lbWriteSpecificDecay',

# the label used in calling data
    bsPointLabel = cms.string('offlineBeamSpot::RECO'),
    pVertexLabel = cms.string('offlinePrimaryVertices::RECO'),
    gpCandsLabel = cms.string('selectedTracks'),
    patMuonLabel = cms.string('selectedMuons'),
    dedxHrmLabel = cms.string('dedxHarmonic2::RECO'),

# the label of output product
    oniaName      = cms.string('oniaFitted'),
    TkTkName      = cms.string('TkTkFitted'),
    LbToTkTkName  = cms.string('LbToTkTkFitted'),
    #Lam0Name      = cms.string('Lam0Fitted'),
    #LbToLam0Name  = cms.string('LbToLam0Fitted'),
    writeVertex   = cms.bool( True ),
    writeMomentum = cms.bool( True ),
    recoSelect    = cms.VPSet(recoSelect)
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('recoWriteSpecificDecay.root'),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_lbWriteSpecificDecay_*_bphAnalysis",
        "keep *_offlineBeamSpot_*_RECO",
        "keep *_offlinePrimaryVertices_*_RECO",
        "keep *_genParticles__HLT"
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('myfilterpath') )
)
process.myfilterpath = cms.Path(
      process.hltHighLevel
    * process.selectedTracks
    * process.selectedMuons
    * process.lbWriteSpecificDecay
)


process.e = cms.EndPath(process.out)
