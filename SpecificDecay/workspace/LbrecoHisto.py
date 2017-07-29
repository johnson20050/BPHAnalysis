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
    #'file:/home/ltsai/ReceivedFile/DATA/testFile.root'
    #'file:/home/ltsai/Work/LbFrame/TEST/CMSSW_8_0_13_patch1/src/BPHAnalysis/SpecificDecay/workspace/reco_MC.root'
    'file:reco_MC.root'
))
#from BPHAnalysis.SpecificDecay.fileLinkv4 import readFiles
#process.source = cms.Source("PoolSource", fileNames = readFiles)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v2', '')             # 8_0_10  DATA
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')             # 8_0_10  DATA

process.bphHistoSpecificDecay = cms.EDAnalyzer('LbHistoSpecificDecay',
        oniaCandsLabel = cms.string("lbWriteSpecificDecay:oniaFitted:bphAnalysis"      ),
        Lam0CandsLabel = cms.string("lbWriteSpecificDecay:Lam0Cand:bphAnalysis"      ),
     LbToLam0CandsLabel = cms.string("lbWriteSpecificDecay:LbToLam0Fitted:bphAnalysis"),
     LbToTkTkCandsLabel = cms.string("lbWriteSpecificDecay:LbToTkTkFitted:bphAnalysis"),
    #    _lamCandsLabel = cms.string('oniaV0Tracks:Kshort:RECO'),          
    #    _lamCandsLabel = cms.string('onia2MuMuPAT::RECO'),          
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("histo.root"),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.bphHistoSpecificDecay
)



















