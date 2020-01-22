from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from WMCore.Configuration import Configuration
from CRABAPI.RawCommand import crabCommand

import datetime

time = datetime.datetime.today()
timestamp = time.strftime('%d%B')


def submit(quark, dataset, suffix=''):

    if len(suffix) > 0:
        suffix = '_' + suffix

    config = Configuration()
    config.section_('General')
    config.General.requestName = 'MC_nano_'+quark + suffix
    config.General.workArea = 'crab_projects'
    config.General.transferOutputs = True
    config.General.transferLogs = True 


    config.section_('JobType')
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'cmssw_mc_gen/test/GenAnalyzer_cfg.py'

    config.section_('Data')
    config.Data.inputDataset = dataset

    config.Data.inputDBS = 'phys03'
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 10
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.outputDatasetTag = 'GT_stEFT_nano_' + timestamp + suffix
    config.Data.publication = True

    config.section_('User')
    config.User.voGroup = 'dcms'

    config.section_('Site')
    config.Site.storageSite = 'T1_DE_KIT_Disk'

    crabCommand('submit', config=config)

if __name__ == '__main__':

    import sys
    quark = sys.argv[1]
    dataset = sys.argv[2]
    try:
        suffix = sys.argv[3]
    except:
        suffix = ''
    submit(quark, dataset, suffix)
