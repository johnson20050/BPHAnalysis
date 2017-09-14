import FWCore.ParameterSet.Config as cms

process = cms.Process("bphAnalysis")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
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
                src=cms.InputTag("oniaSelectedTracks::RECO"), # BPHSkim
                #src=cms.InputTag("generalTracks::RECO"),     # AOD
                particleType=cms.string('pi+')
)

from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.patSelectedTracks = patGenericParticles.clone(src=cms.InputTag("CandidateSelectedTracks"))

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
'file:///afs/cern.ch/user/l/ltsai/ReceivedFile/DATA/CMSSW_8_0_28/2016RunB_Charmonium/14C79216-F03D-E711-AE2F-D067E5F91695.root'
))
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016LegacyRepro_v3', '')

#from BPHAnalysis.SpecificDecay.LbrecoSelectForWrite_cfi import recoSelect
from BPHAnalysis.SpecificDecay.newSelectForWrite_cfi import recoSelect


HLTName='HLT_DoubleMu4_JpsiTrk_Displaced_v*'
from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltHighLevel= hltHighLevel.clone(HLTPaths = cms.vstring(HLTName))

process.lbWriteSpecificDecay = cms.EDProducer('lbWriteSpecificDecay',

# the label used in calling data
    pVertexLabel = cms.string('offlinePrimaryVertices::RECO'),
    #patMuonLabel = cms.string('selectedPatMuons'), # AOD
    ccCandsLabel = cms.string('onia2MuMuPAT::RECO'),# BPHSkim
    gpCandsLabel = cms.string('patSelectedTracks'), # BPHSkim && AOD

# the label of output product
    oniaName      = cms.string('oniaFitted'),
    Lam0Name      = cms.string('Lam0Cand'),
    TkTkName      = cms.string('TkTkFitted'),
    LbToLam0Name  = cms.string('LbToLam0Fitted'),
    LbToTkTkName  = cms.string('LbToTkTkFitted'),
    PhiName   = cms.string('phiCand'),
    BsName   = cms.string('bsFitted'),
    writeVertex   = cms.bool( True ),
    writeMomentum = cms.bool( True ),
    recoSelect    = cms.VPSet(recoSelect)
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('recoWriteSpecificDecay.root'),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_lbWriteSpecificDecay_*_bphAnalysis"
        "keep *_offlineBeamSpot_*_RECO"
        "keep *_offlinePrimaryVertices_*_RECO"
    )
)

process.p = cms.Path(
    process.hltHighLevel *
    process.CandidateSelectedTracks *
    process.lbWriteSpecificDecay
)

process.e = cms.EndPath(process.out)

