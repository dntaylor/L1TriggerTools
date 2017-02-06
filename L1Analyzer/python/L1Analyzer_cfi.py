import FWCore.ParameterSet.Config as cms

l1Analyzer = cms.EDAnalyzer('L1Analyzer',
    ecalDigis = cms.InputTag("ecalDigis","EcalTriggerPrimitives"),
    hcalDigis = cms.InputTag("hcalDigis"),
    stage2Layer1Digis = cms.InputTag("stage2Layer1Digis"),
)
