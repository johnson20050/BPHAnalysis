import FWCore.ParameterSet.Config as cms

process = cms.Process("bphAnalysis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

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

process.CandidateSelectedTracks = cms.EDProducer( "ConcreteChargedCandidateProducer",
                src=cms.InputTag("oniaSelectedTracks::RECO"),   # BPHSkim
                #src=cms.InputTag("generalTracks::RECO"),       # AOD
                particleType=cms.string('pi+')
)

from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.patSelectedTracks = patGenericParticles.clone(src=cms.InputTag("CandidateSelectedTracks"))

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
    #'file:/home/ltsai/ReceivedFile/BPHanalysis/BPHanalysisTest.root'                   # 8_0_8 BPHSkim
     'file:/home/ltsai/ReceivedFile/DATA/BPHSkim_Charmonium.root'                       # 8_0_8 BPHSkim Charmonium
))

#from BPHAnalysis.SpecificDecay.fileLinkv3 import readFiles
#process.source = cms.Source("PoolSource", fileNames = readFiles)



from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v2', '')             # 8_0_10  DATA

from BPHAnalysis.testAnalyzer.LbrecoSelectForWrite_cfi import recoSelect

process.lbSpecificDecay = cms.EDAnalyzer('LbSpecificDecay',

# the label used in calling data
    pVertexLabel = cms.string('offlinePrimaryVertices::RECO'),
    #patMuonLabel = cms.string('selectedPatMuons'),
    gpCandsLabel = cms.string('patSelectedTracks'),
    #pfCandsLabel = cms.string('particleFlow'),
    ccCandsLabel = cms.string('onia2MuMuPAT::RECO'),

    oniaName      = cms.string('oniaFitted'),
# the label of output product
    Lam0Name      = cms.string('Lam0Cand'),
    LbToLam0Name  = cms.string('LbToLam0Fitted'),
    LbToTkTkName  = cms.string('LbToTkTkFitted'),
    writeVertex   = cms.bool( True ),
    writeMomentum = cms.bool( True ),
    recoSelect    = cms.VPSet(recoSelect)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('result.root')
)

process.p = cms.Path(
    process.CandidateSelectedTracks *
    process.lbSpecificDecay
)


