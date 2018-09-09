import FWCore.ParameterSet.Config as cms


#--------------------------------------------------------------------------------
# - Muon  Setting
#--------------------------------------------------------------------------------
# prepare pat::Muon from AOD
from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff import *
from PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff import *


selectedPatMuons.cut = cms.string('muonID(\"TMOneStationTight\")'
    ' && abs(innerTrack.dxy) < 0.3'
    ' && abs(innerTrack.dz)  < 20.'
    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
    ' && innerTrack.quality(\"highPurity\")'
)

#myMuonsSequence = cms.EDFilter(
selectedMuons = cms.EDFilter(
    "MuonProducer",
    muonsrc = cms.InputTag("selectedPatMuons"),
    )


#--------------------------------------------------------------------------------
# - Track Setting
#--------------------------------------------------------------------------------
# prepare GenericParticle from AOD
CandidateSelectedTracks = cms.EDProducer( "ConcreteChargedCandidateProducer",
                #src=cms.InputTag("oniaSelectedTracks::RECO"),
                src=cms.InputTag("generalTracks::RECO"),
                particleType=cms.string('pi+')
)
from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
patSelectedTracks = patGenericParticles.clone(src=cms.InputTag("CandidateSelectedTracks"))


# preselect track
#myTrackSequence = cms.EDFilter(
selectedTracks = cms.EDFilter(
    "TrackProducer",
    tracksrc = cms.InputTag("patSelectedTracks"),
    muonsrc = cms.InputTag("selectedPatMuons")
    )


