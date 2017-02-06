import FWCore.ParameterSet.Config as cms

######################
### Option parsing ###
######################
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'l1tTree.root'
options.inputFiles = '/store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-DIGI-RAW/81X_upgrade2017_realistic_v22-v1/10000/06AEB644-DFA6-E611-88B9-0CC47A4C8EBA.root'
options.maxEvents = -1
options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Events to skip")
options.register('reportEvery', 100, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Report every")

options.parseArguments()

##########################
### process definition ###
##########################
from Configuration.StandardSequences.Eras import eras
process = cms.Process("L1Analyzer",eras.Run2_2017)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    skipEvents = cms.untracked.uint32(options.skipEvents),
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(options.outputFile),
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

################
### emulator ###
################
process.load('L1Trigger.L1TCaloLayer1.simCaloStage2Layer1Digis_cfi')
process.simCaloStage2Layer1Digis.useECALLUT = cms.bool(True)
process.simCaloStage2Layer1Digis.useHCALLUT = cms.bool(True)
process.simCaloStage2Layer1Digis.useHFLUT = cms.bool(True)
process.simCaloStage2Layer1Digis.useLSB = cms.bool(True)
process.simCaloStage2Layer1Digis.verbose = cms.bool(True)
#process.simCaloStage2Layer1Digis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
#process.simCaloStage2Layer1Digis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")
process.simCaloStage2Layer1Digis.ecalToken = cms.InputTag("simEcalTriggerPrimitiveDigis")
process.simCaloStage2Layer1Digis.hcalToken = cms.InputTag("simHcalTriggerPrimitiveDigis")

#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_2017_v1_0_inconsistent_cfi.py')

############################
### ntuple customization ###
############################
process.load('L1TriggerTools.L1Analyzer.L1Analyzer_cfi')
process.l1Analyzer.ecalDigis = cms.InputTag("simEcalTriggerPrimitiveDigis")
process.l1Analyzer.hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
process.l1Analyzer.stage2Layer1Digis = cms.InputTag("simCaloStage2Layer1Digis")

process.path = cms.Path(
#    process.ecalDigis
#    +process.hcalDigis
    process.simCaloStage2Layer1Digis
#    +process.simCaloStage2Digis
    +process.l1Analyzer
)
