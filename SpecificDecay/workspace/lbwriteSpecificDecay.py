import FWCore.ParameterSet.Config as cms

process = cms.Process("bphAnalysis")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

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
                src=cms.InputTag("oniaSelectedTracks::RECO"),
                particleType=cms.string('pi+')
)

from PhysicsTools.PatAlgos.producersLayer1.genericParticleProducer_cfi import patGenericParticles
process.patSelectedTracks = patGenericParticles.clone(src=cms.InputTag("CandidateSelectedTracks"))

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
'file:///home/ltsai/Data/8_0_20/001E3B7A-E784-E611-B6C8-008CFA05E8EC.root',
'file:///home/ltsai/Data/8_0_20/02860492-C385-E611-8224-008CFAF0842A.root',
'file:///home/ltsai/Data/8_0_20/0431281D-2B86-E611-8C8F-02163E012ED5.root',
'file:///home/ltsai/Data/8_0_20/0476138F-AC8A-E611-B54C-02163E017657.root',
'file:///home/ltsai/Data/8_0_20/068A3275-6886-E611-B337-02163E011591.root',
'file:///home/ltsai/Data/8_0_20/08A02647-BD85-E611-B736-1CC1DE19274E.root',
'file:///home/ltsai/Data/8_0_20/08C1A145-B985-E611-886E-FA163E6D8C1A.root',
'file:///home/ltsai/Data/8_0_20/0A7D3A63-DB85-E611-823A-FA163E9639D6.root',
'file:///home/ltsai/Data/8_0_20/0C02D086-0986-E611-8BFD-FA163E121CF1.root',
'file:///home/ltsai/Data/8_0_20/0C03F3CA-FC84-E611-9B6F-00259073E544.root',
'file:///home/ltsai/Data/8_0_20/0C8D4BF8-2687-E611-BF99-1CB72C1B2D88.root',
'file:///home/ltsai/Data/8_0_20/0CDCAE2C-8286-E611-B5B5-001E6750489D.root',
'file:///home/ltsai/Data/8_0_20/0E0BEEF0-D085-E611-BF45-02163E015F72.root',
'file:///home/ltsai/Data/8_0_20/10726B43-E885-E611-A059-FA163E35DE22.root',
'file:///home/ltsai/Data/8_0_20/10C1295A-4686-E611-81D2-FA163EB07D21.root',
'file:///home/ltsai/Data/8_0_20/10E1F5AD-7486-E611-BF5E-001E67C9AF38.root',
'file:///home/ltsai/Data/8_0_20/140B4175-7685-E611-A67E-001E67F8FA3D.root',
'file:///home/ltsai/Data/8_0_20/14676347-0286-E611-ADEF-FA163E25D48A.root',
'file:///home/ltsai/Data/8_0_20/14752985-5986-E611-B6AD-FA163EB7F610.root',
'file:///home/ltsai/Data/8_0_20/14AE588F-E285-E611-B6D5-FA163E4A6D71.root',
'file:///home/ltsai/Data/8_0_20/16018E40-9C86-E611-A135-842B2B7680DF.root',
'file:///home/ltsai/Data/8_0_20/1625CB6B-CF85-E611-885F-FA163EA6F910.root',
'file:///home/ltsai/Data/8_0_20/163FFEAE-9386-E611-86C7-002590E7E010.root',
'file:///home/ltsai/Data/8_0_20/16516944-B985-E611-9C36-FA163E8CCFC6.root',
'file:///home/ltsai/Data/8_0_20/1675C2A6-2585-E611-925A-001E6757E03C.root',
'file:///home/ltsai/Data/8_0_20/588FCC5E-AC8A-E611-9FF6-FA163EA08A8F.root',
'file:///home/ltsai/Data/8_0_20/5AAF92D3-F88B-E611-9D0F-02163E01315D.root',
'file:///home/ltsai/Data/8_0_20/68BCDE2F-2C8B-E611-AD55-0CC47A4D76A0.root',
'file:///home/ltsai/Data/8_0_20/88E3FF56-AC8A-E611-B4A5-FA163E3B08F9.root',
'file:///home/ltsai/Data/8_0_20/90754EA5-D38B-E611-AE3A-FA163EF1361B.root',
'file:///home/ltsai/Data/8_0_20/9CF27DC0-BF8A-E611-99C0-FA163E0ED7B3.root',
'file:///home/ltsai/Data/8_0_20/A614E7DE-938A-E611-9222-FA163EEDA31A.root',
'file:///home/ltsai/Data/8_0_20/B0B246A2-878A-E611-8325-FA163E8DCC10.root',
'file:///home/ltsai/Data/8_0_20/D21D2032-2C8B-E611-BE4D-0025905A60AA.root',
'file:///home/ltsai/Data/8_0_20/D24D02DC-D38B-E611-93FE-02163E01765D.root',
'file:///home/ltsai/Data/8_0_20/DAF64466-AC8A-E611-A804-FA163E51EF69.root',
'file:///home/ltsai/Data/8_0_20/E058024B-AC8A-E611-AD9D-FA163E32AAAE.root',
'file:///home/ltsai/Data/8_0_20/F8974060-AC8A-E611-86E3-FA163EAA66FB.root',
'file:///home/ltsai/Data/8_0_20/FCF396F4-938A-E611-AFB4-FA163E65B77A.root',
))

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

from BPHAnalysis.SpecificDecay.LbrecoSelectForWrite_cfi import recoSelect

process.lbWriteSpecificDecay = cms.EDProducer('lbWriteSpecificDecay',
    pVertexLabel = cms.string('offlinePrimaryVertices::RECO'),
    gpCandsLabel = cms.string('patSelectedTracks'),
    ccCandsLabel = cms.string('onia2MuMuPAT::RECO'),
    Lam0Name      = cms.string('Lam0Cand'),
    TkTkName      = cms.string('TkTkCand'),
    LbToLam0Name  = cms.string('LbToLam0Fitted'),
    LbToTkTkName  = cms.string('LbToTkTkFitted'),
    writeVertex   = cms.bool( True ),
    writeMomentum = cms.bool( True ),
    recoSelect = cms.VPSet(recoSelect)
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('reco.root'),
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
    process.CandidateSelectedTracks *
    process.lbWriteSpecificDecay
)

process.e = cms.EndPath(process.out)

