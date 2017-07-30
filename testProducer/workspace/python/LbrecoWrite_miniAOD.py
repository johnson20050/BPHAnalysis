import FWCore.ParameterSet.Config as cms

process = cms.Process("bphAnalysis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

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


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
    'file:///home/ltsai/Data/miniAOD_LambdaBToLambdaMumuToPPiMuMu_EtaPtFilter_13Tev-pythia6-evtgen.root'    # MC
    #'file:///home/ltsai/Data/miniAOD_data/6E3E2983-743B-E611-B6FF-02163E01417D.root'
    #'file:///home/ltsai/Data/miniAOD_Charmonium.root'                                                       # Data
))




from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')             # MC
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_ICHEP16_repro_v0', '')             # 8_0_10  DATA

from BPHAnalysis.testProducer.LbrecoSelectForWrite_cfi import recoSelect
process.TFileService = cms.Service("TFileService",
    #fileName = cms.string('result_MCminiAOD.root')
    fileName = cms.string('result_DATAminiAOD.root')
)


process.lbSpecificDecay = cms.EDAnalyzer('testProducer',

# the label used in calling data
   #pVertexLabel = cms.string('offlineSlimmedPrimaryVertices::RECO'),
   #patMuonLabel = cms.string('slimmedMuons::RECO'),
   #pcCandsLabel = cms.string('packedPFCandidates::RECO'),
    pVertexLabel = cms.string('offlineSlimmedPrimaryVertices::PAT'),
    patMuonLabel = cms.string('slimmedMuons::PAT'),
    pcCandsLabel = cms.string('packedPFCandidates::PAT'),

    oniaName      = cms.string('oniaFitted'),
# the label of output product
    Lam0Name      = cms.string('Lam0Cand'),
    LbToLam0Name  = cms.string('LbToLam0Fitted'),
    LbToTkTkName  = cms.string('LbToTkTkFitted'),
    writeVertex   = cms.bool( True ),
    writeMomentum = cms.bool( True ),
    recoSelect    = cms.VPSet(recoSelect)
)


process.p = cms.Path(
    process.lbSpecificDecay
)


