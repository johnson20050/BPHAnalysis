import FWCore.ParameterSet.Config as cms

process = cms.Process("bphAnalysis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

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
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
'file:/home/ltsai/Data/miniAOD_LambdaBToLambdaMumuToPPiMuMu_EtaPtFilter_13Tev-pythia6-evtgen.root'
))

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

from BPHAnalysis.SpecificDecay.recoSelectForWrite_cfi import recoSelect

process.bphWriteSpecificDecay = cms.EDProducer('BPHWriteSpecificDecay',
    pVertexLabel = cms.string('offlineSlimmedPrimaryVertices::PAT'),
    patMuonLabel = cms.string('slimmedMuons::PAT'),
    pcCandsLabel = cms.string('packedPFCandidates::PAT'),
    oniaName = cms.string('oniaFitted'),
    sdName   = cms.string('kx0Cand'),
    ssName   = cms.string('phiCand'),
    buName   = cms.string('buFitted'),
    bdName   = cms.string('bdFitted'),
    bsName   = cms.string('bsFitted'),
    writeVertex   = cms.bool( True ),
    writeMomentum = cms.bool( True ),
    recoSelect = cms.VPSet(recoSelect)
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('reco.root'),
    outputCommands = cms.untracked.vstring(
      "keep *",
      "keep *_writeBPHSpecificDecay_*_*",
      "drop *_patSelectedTracks_*_*",
      "drop *_CandidateSelectedTracks_*_*",
      "drop *_TriggerResults_*_bphAnalysis",
      "drop *_random*_*_bphAnalysis"
    ),
)

process.p = cms.Path(
    process.bphWriteSpecificDecay
)

process.e = cms.EndPath(process.out)

