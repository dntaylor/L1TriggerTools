import FWCore.ParameterSet.Config as cms

l1Analyzer = cms.EDAnalyzer("L1Analyzer",
    ecalDigis = cms.InputTag("ecalDigis","EcalTriggerPrimitives"),
    hcalDigis = cms.InputTag("hcalDigis"),
    stage2Layer1Digis = cms.InputTag("stage2Layer1Digis"),
    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    electronIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
    electronPairs = cms.InputTag("TODO"),
    genParticles = cms.InputTag("genParticles"),
    storeEcal = cms.bool(True),
    storeHcal = cms.bool(True),
    storeStage2Layer1 = cms.bool(True),
    LSB = cms.double(0.5),
    isMC = cms.bool(False),
)
