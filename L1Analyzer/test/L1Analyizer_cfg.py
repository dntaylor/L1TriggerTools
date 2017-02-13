import FWCore.ParameterSet.Config as cms

######################
### Option parsing ###
######################
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'l1tTree.root'
#options.inputFiles = '/store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-DIGI-RAW/81X_upgrade2017_realistic_v22-v1/10000/06AEB644-DFA6-E611-88B9-0CC47A4C8EBA.root'
options.inputFiles = '/store/data/Run2016E/SingleElectron/MINIAOD/23Sep2016-v1/100000/00827A71-F98C-E611-9639-0CC47A4D7650.root'
options.maxEvents = -1
options.register('isMC', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Simulation")
options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Events to skip")
options.register('reportEvery', 100, VarParsing.multiplicity.singleton, VarParsing.varType.int, "Report every")
options.register('lumimask', '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt', VarParsing.multiplicity.singleton, VarParsing.varType.string, "Lumimask for data")

options.parseArguments()

##########################
### process definition ###
##########################
from Configuration.StandardSequences.Eras import eras
process = cms.Process("L1Analyzer",eras.Run2_2016)

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    skipEvents = cms.untracked.uint32(options.skipEvents),
)

def getSecondaryFiles(primaryFileList) :
    secondaryFiles = []
    for primaryFile in primaryFileList :
        query = 'parent file=%s' % primaryFile
        for entry in subprocess.Popen("das_client --query='parent file={0}' --limit=0".format(primaryFile), shell=True, stdout=subprocess.PIPE).communicate()[0].splitlines():
            secondaryFiles.append(entry)
    return secondaryFiles

if not options.isMC:
    import subprocess
    process.source.secondaryFileNames = cms.untracked.vstring(getSecondaryFiles(process.source.fileNames))

# output files
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(options.outputFile),
)

# globaltag
envvar = 'mcgt' if options.isMC else 'datagt'
from Configuration.AlCa.GlobalTag import GlobalTag
GT = {'mcgt': 'auto:run2_mc', 'datagt': 'auto:run2_data'}
process.GlobalTag = GlobalTag(process.GlobalTag, GT[envvar], '')

# if data, apply lumimask
if not options.isMC:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(
        filename = options.lumimask
    ).getVLuminosityBlockRange()


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

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')      # CaloTPGTranscoder
process.load('L1Trigger.L1TCalorimeter.caloStage2Params_2016_v3_3_cfi')  # L1 CaloParams
process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

from L1Trigger.Configuration.customiseReEmul import L1TEventSetupForHF1x1TPs,L1TReEmulFromRAW,L1TReEmulFromRAWsimTP
process = L1TEventSetupForHF1x1TPs(process)

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')
process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
process.simHcalTriggerPrimitiveDigis.inputLabel = cms.VInputTag(
    cms.InputTag('hcalDigis'),
    cms.InputTag('hcalDigis')
    )
process.L1TReEmul = cms.Sequence(process.simHcalTriggerPrimitiveDigis * process.simEcalTriggerPrimitiveDigis * process.simCaloStage2Layer1Digis)

##############################
### electron customization ###
##############################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.tagElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("(abs(-log(tan(superCluster.position.theta/2)))<=2.1) && !(1.4442<=abs(-log(tan(superClusterPosition.theta/2)))<=1.566) && pt >= 30.0"),
    #filter = cms.bool(True)
)

process.probeElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string(""),
)

process.electronPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagElectrons@+ probeElectrons@-"),
    cut   = cms.string("76 < mass < 106")
)

process.electronTnP = cms.Sequence(process.tagElectrons * process.probeElectrons * process.electronPairs)

############################
### ntuple customization ###
############################
process.load('L1TriggerTools.L1Analyzer.L1Analyzer_cfi')
process.l1Analyzer.ecalDigis = cms.InputTag("simEcalTriggerPrimitiveDigis")
process.l1Analyzer.hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
process.l1Analyzer.stage2Layer1Digis = cms.InputTag("simCaloStage2Layer1Digis")
process.l1Analyzer.electronIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose")
process.l1Analyzer.electronPairs = cms.InputTag("electronPairs")
process.l1Analyzer.storeEcal = cms.bool(False)
process.l1Analyzer.storeHcal = cms.bool(False)
process.l1Analyzer.storeStage2Layer1 = cms.bool(False)

##################
### final path ###
##################
process.path = cms.Path(
    process.RawToDigi
    *process.L1TReEmul
    *process.egmGsfElectronIDSequence
    *process.electronTnP
    *process.l1Analyzer
)
