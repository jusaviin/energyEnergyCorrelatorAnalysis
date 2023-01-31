from WMCore.Configuration import Configuration
config = Configuration()

isMC='1'
jobTag='skimjob_PbPb2018MC_test_2020-03-20'
inputList='PbPbMC2018_commonForest.txt'
outputFile='HiForestAOD_PbPbMC2018skim.root'
fileLocation='1'  # Locations: 0 = Purdue, 1 = CERN, 2 = Search with xrootd

config.section_("General")
config.General.requestName = jobTag
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
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
config.Data.outputPrimaryDataset = 'dijetSkim'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_*']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
#config.section_("Debug")
#config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

