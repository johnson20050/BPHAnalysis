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

process.CandidateSelectedTracks = cms.EDProducer( "ConcreteChargedCandidateProducer",
                #src=cms.InputTag("oniaSelectedTracks::RECO"),
                src=cms.InputTag("generalTracks::RECO"),
                particleType=cms.string('pi+')
)

from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.patSelectedTracks = patGenericParticles.clone(src=cms.InputTag("CandidateSelectedTracks"))

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
    '',
))

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

from BPHAnalysis.SpecificDecay.LbrecoSelectForWrite_cfi import recoSelect

# create pat::Muons {{{
process.options.allowUnscheduled = cms.untracked.bool(True)
process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff')
process.load('PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff')

process.selectedPatMuons.cut = cms.string('muonID(\"TMOneStationTight\")'
    ' && abs(innerTrack.dxy) < 0.3'
    ' && abs(innerTrack.dz)  < 20.'
    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
    ' && innerTrack.quality(\"highPurity\")'
)

#make patTracks
from PhysicsTools.PatAlgos.tools.trackTools import makeTrackCandidates
makeTrackCandidates(process,
    label        = 'TrackCands',                  # output collection
    tracks       = cms.InputTag('generalTracks'), # input track collection
    particleType = 'pi+',                         # particle type (for assigning a mass)
    preselection = 'pt > 0.7',                    # preselection cut on candidates
    selection    = 'pt > 0.7',                    # selection on PAT Layer 1 objects
    isolation    = {},                            # isolations to use (set to {} for None)
    isoDeposits  = [],
    mcAs         = None                           # replicate MC match as the one used for Muons
)
process.patTrackCands.embedTrack = True
# create pat::Muons end }}}



process.lbWriteSpecificDecay = cms.EDProducer('lbWriteSpecificDecay',

# the label used in calling data
    pVertexLabel = cms.string('offlinePrimaryVertices::RECO'),
    patMuonLabel = cms.string('selectedPatMuons'),
    #gpCandsLabel = cms.string('patSelectedTracks'),
    #pfCandsLabel = cms.string('particleFlow'),
    #ccCandsLabel = cms.string('onia2MuMuPAT::RECO'),

# the label of output product
    oniaName      = cms.string('oniaFitted'),
    Lam0Name      = cms.string('Lam0Cand'),
    LamxName      = cms.string('LamxCand'),
    TkTkName      = cms.string('TkTkFitted'),
    LbToLam0Name  = cms.string('LbToLam0Fitted'),
    LbToTkTkName  = cms.string('LbToTkTkFitted'),
    writeVertex   = cms.bool( True ),
    writeMomentum = cms.bool( True ),
    recoSelect    = cms.VPSet(recoSelect)
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('recoWriteSpecificDecay.root'),
    outputCommands = cms.untracked.vstring(
      "keep *",
      "keep *_lbWriteSpecificDecay_*_bphAnalysis",
      "drop *_patSelectedTracks_*_*",
      "drop *_CandidateSelectedTracks_*_*",
      "drop *_TriggerResults_*_bphAnalysis",
      "drop *_random*_*_bphAnalysis"
    ),
)

process.p = cms.Path(
    process.lbWriteSpecificDecay
)

process.e = cms.EndPath(process.out)

