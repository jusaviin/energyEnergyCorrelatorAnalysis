from WMCore.Configuration import Configuration
config = Configuration()

isMC='0'
jobTag='PbPb2018skim_HardProbes_MiniAODv1-v1_jet80or100Trigger_2023-01-31'
inputList='PbPb2018Data_miniAOD.txt'
outputFile='HiForestMiniAOD_PbPb2018skim.root'
fileLocation='2'  # Locations: 0 = Purdue, 1 = CERN, 2 = Vanderbilt, 3 = Search with xrootd

config.section_("General")
config.General.requestName = jobTag
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRunSkimmer.sh'
config.JobType.scriptArgs = ['ismc='+isMC,'output='+outputFile,'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','skimmer.tar.gz']
config.JobType.outputFiles = [outputFile]
config.JobType.maxJobRuntimeMin = 120
config.JobType.maxMemoryMB = 1000

config.section_("Data")
config.Data.userInputFiles = open(inputList).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'dataSkim'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = 'T2_CH_CERN'

#"really" force crab to only run at whitelisted sites
#config.section_("Debug")
#config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

