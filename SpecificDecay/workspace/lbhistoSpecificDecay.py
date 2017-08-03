import FWCore.ParameterSet.Config as cms

process = cms.Process("bckAnalysis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
    'file:reco.root'
))

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.TFileService = cms.Service('TFileService',
  fileName = cms.string('his.root'),
  closeFileFast = cms.untracked.bool(True)
)

process.lbHistoSpecificDecay = cms.EDAnalyzer('lbHistoSpecificDecay',
    oniaCandsLabel =         cms.string('lbWriteSpecificDecay:oniaCand:bphAnalysis'),
    LbToJPsiLam0CandsLabel = cms.string('lbWriteSpecificDecay:LbToLam0Fitted:bphAnalysis'),
    LbToJPsiTkTkCandsLabel = cms.string('lbWriteSpecificDecay:LbToTkTkFitted:bphAnalysis'),
)

process.p = cms.Path(
    process.lbHistoSpecificDecay
)


