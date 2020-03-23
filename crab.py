from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABAPI.RawCommand import crabCommand
from WMCore.Configuration import Configuration

import datetime

time = datetime.datetime.today()
timestamp = time.strftime('%d%B')

def submit(quark, suffix=None):

    suffix = '' if suffix is None else '_' + suffix
    tar_balls = {
            'top' : 'st_tch_top_halfchain_slc7_amd64_gcc700_CMSSW_10_2_17_tarball.tar.xz'
    }

    # tar_ball = tar_balls[quark+suffix]
    tar_ball = 'st_tch_top_halfchain_ctw_5_slc7_amd64_gcc700_CMSSW_10_2_17_tarball.tar.xz'

    config = Configuration()

    config.section_('General')
    config.General.requestName = 'MC_generation_top_halfchain_ctw_5'
    config.General.workArea = 'crab_projects_2'
    config.General.transferOutputs = True
    config.General.transferLogs = True

    config.section_('JobType')
    config.JobType.pluginName = 'PrivateMC'
    config.JobType.psetName = 'TOP_gen_halfchain_ctw_5_crab_cfg.py'
    # config.JobType.inputFiles = ['GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh', 'st_tch_top_slc7_amd64_gcc700_CMSSW_10_2_17_tarball.tar.xz']
    # config.JobType.inputFiles = ['GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh', 'st_tch_top_slc7_amd64_gcc700_CMSSW_10_2_17_tarball.tar.xz']
    config.JobType.inputFiles = [tar_ball]

    config.section_('Data')
    config.Data.outputPrimaryDataset = 'ST_t-channel_4f_elmuDecays_13TeV-amcatnloFXFX-pythia8'
    config.Data.splitting = 'EventBased'
    config.Data.unitsPerJob = 10000
    NJOBS = 10  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
    config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.publication = True
    config.Data.outputDatasetTag = 'GT_stEFT_' + timestamp + suffix

    config.section_('User')
    config.User.voGroup = 'dcms'

    config.section_('Site')
    config.Site.storageSite = 'T1_DE_KIT_Disk'

    crabCommand('submit', config=config)

if __name__ == '__main__':

    import sys
    quark = sys.argv[1]

    try:
        suffix = sys.argv[2]
    except:
        suffix = None

    print(quark, suffix)
    submit(quark, suffix)
