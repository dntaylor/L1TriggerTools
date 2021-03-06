import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.StandardSequences.Eras import eras
process = cms.Process("L1TCaloLayer1LUTWriter",eras.Run2_2016)

options = VarParsing()
options.register('runNumber', 270217, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Run to analyze')
options.register('outputFile', 'luts.xml', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Output XML File')
options.parseArguments()

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2017_realistic_v20', '')

# temporary
#process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2017_realistic_v15', '')
#
#process.es_ascii = cms.ESSource("HcalTextCalibrations",
#    input = cms.VPSet(
#        cms.PSet(
#            object = cms.string('LutMetadata'),
#            #file = cms.FileInPath('DumpCondLutMetadata_Run999999.NewRcalib.txt')
#            file = cms.FileInPath('L1TriggerTools/L1Analyzer/data/DumpCondLutMetadata_Run999999.NewRcalib.txt')
#        )
#    )
#)
#process.es_prefer = cms.ESPrefer('HcalTextCalibrations','es_ascii')

process.source = cms.Source('EmptySource',
    #firstRun = cms.untracked.uint32(options.runNumber)
)

# Writes LUT for the only event to be processed - ignores data itself.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

#process.load('L1Trigger.L1TCaloLayer1Spy.l1tCaloLayer1LUTWriter_cfi')
#process.l1tCaloLayer1LUTWriter.fileName = options.outputFile
#
#process.p = cms.Path(process.l1tCaloLayer1LUTWriter)

process.load('L1TriggerTools.L1TCaloLayer1Spy.l1tCaloLayer1LUTWriter_cfi')
process.l1tCaloLayer1LUTWriterNew.fileName = options.outputFile

process.p = cms.Path(process.l1tCaloLayer1LUTWriterNew)

process.schedule = cms.Schedule(process.p)

# HF 1x1 TPs need special ES source as of 2016/03/21
#from L1Trigger.Configuration.customiseReEmul import L1TEventSetupForHF1x1TPs
#process = L1TEventSetupForHF1x1TPs(process)

# To get L1 CaloParams, until in GT
#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_2017_v2_1_cfi')
process.load('L1Trigger.L1TCalorimeter.caloStage2Params_2017_v2_1_inconsistent_cfi')

nEt = len(process.caloStage2Params.layer1HCalScaleETBins)
for ieta in range(17,28):
    for ipt in range(nEt):
         pass
         #process.caloStage2Params.layer1HCalScaleFactors[nEt*28+ipt*28+ieta] *= 0.  # zeroes
         #process.caloStage2Params.layer1HCalScaleFactors[nEt*28+ipt*28+ieta] = 1.   # ones
         #process.caloStage2Params.layer1HCalScaleFactors[nEt*28+ipt*28+ieta] *= 2.  # double
         #process.caloStage2Params.layer1HCalScaleFactors[nEt*28+ipt*28+ieta] *= 0.5 # halve

# To get CaloTPGTranscoder
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
