import FWCore.ParameterSet.Config as cms

### This is a simpler version. You should edit the fileNames and apply cmsRun.
### Ref: http://cmslxr.fnal.gov/lxr/source/PhysicsTools/HepMCCandAlgos/plugins/ParticleListDrawer.cc

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    'file:///home/ltsai/Data/miniAOD_LambdaBToLambdaMumuToPPiMuMu_EtaPtFilter_13Tev-pythia6-evtgen.root'

    )
    #eventsToProcess = cms.untracked.VEventRange('1:8-1:6','2:100-3:max'),
    #eventsToSkip = cms.untracked.VEventRange('1:1-1:6','2:100-3:max'),
    #lumisToProcess = cms.untracked.VLuminosityBlockRange('1:1-1:6','2:100-3:max'),
    #lumisToSkip = cms.untracked.VLuminosityBlockRange('1:1-1:6','2:100-3:max'),
)

process.printGenParticle = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(-1)
)

##
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("genParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(False),
    printIndex = cms.untracked.bool(False),
    status = cms.untracked.vint32(3)
)

##
process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src = cms.InputTag("genParticles"),
    printP4 = cms.untracked.bool(True),
    printPtEtaPhi = cms.untracked.bool(True),
    printVertex = cms.untracked.bool(True),
    status = cms.untracked.vint32(3)
)

### Optional: check all charged particles
#process.genParticlesClone = cms.EDFilter("CandViewShallowCloneProducer",
#    src = cms.InputTag("genParticles"),
#    cut = cms.string("!charge = 0 & status=1")
#)

#process.chargedHistos= cms.EDAnalyzer("CandViewHistoAnalyzer",
#    src = cms.InputTag("genParticlesClone"),
#    histograms = cms.VPSet(
#        cms.PSet(
#            min = cms.untracked.double(0.0),
#            max = cms.untracked.double(20.0),
#            nbins = cms.untracked.int32(50),
#            name = cms.untracked.string("charged particle pT"),
#            description = cms.untracked.string("pT [GeV/c]"),
#            plotquantity = cms.untracked.string("pt")
#        )
#    )
#)
#
#process.TFileService = cms.Service(
#        "TFileService",
#        fileName = cms.string("charged.root")
#)

#process.p = cms.Path(process.printGenParticle)
process.p = cms.Path(process.printDecay*process.printGenParticle)
#process.p = cms.Path(process.genParticlesClone*process.printGenParticle*process.chargedHistos)
