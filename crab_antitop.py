from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

import datetime

time = datetime.datetime.today()
timestamp = time.strftime('%d%B')


config.General.requestName = 'MC_generation_antitop'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'TOP_gen_antitop_cfg.py'
# config.JobType.inputFiles = ['GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh', 'st_tch_top_slc7_amd64_gcc700_CMSSW_10_2_17_tarball.tar.xz']
# config.JobType.inputFiles = ['GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh', 'st_tch_top_slc7_amd64_gcc700_CMSSW_10_2_17_tarball.tar.xz']
config.JobType.inputFiles = ['st_tch_antitop_slc7_amd64_gcc700_CMSSW_10_2_17_tarball.tar.xz']

config.Data.outputPrimaryDataset = 'ST_t-channel_antitop_4f_elmuDecays_13TeV-amcatnloFXFX-pythia8'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 10000
NJOBS = 100  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'GT_stEFT_' + timestamp

config.User.voGroup = 'dcms'
config.Site.storageSite = 'T1_DE_KIT_Disk'
