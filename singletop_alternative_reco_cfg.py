import FWCore.ParameterSet.Config as cms

# process = cms.Process("GenAnalyzer")
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
process = cms.Process('NANO', Run2_2018)

# process.load('PhysicsTools.NanoAOD.nano_cff')
# process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load('Configuration.StandardSequences.EndOfProcess_cff')
# process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

import sys

infile = sys.argv[2]

# file_name = infile.rsplit('/')[-1].replace('.txt', '.root')
# 
# rfiles = []
# with open(infile) as f:
#     for line in f:
#         rfiles.append('file:' + line.rstrip())

# file_path = '/ceph/dseith/mc_gen/'#'/ceph/dseith/mc_gen/more_couplings/bms1/'
# file_path = ''
rfiles = infile
print(sys.argv)
file_name = sys.argv[3]

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
                    # 'file:' + file_path + 'TOP_gen_top.root'
                    # 'file:' + file_path + 'TOP_gen_31.root'
                    # 'file:'+infile,
                    rfiles
                    # 'file:/ceph/dseith/pdf/files/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5down_PSweights_13TeV-powheg-pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14-v1_numEvent1.root'

                )
                            )
process.options = cms.untracked.PSet(

)

from  PhysicsTools.NanoAOD.common_cff import *

process.nanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

process.genWeightsTable = cms.EDProducer('GenWeightsTableProducer',
        genEvent = cms.InputTag('generator'),
        # genLumiInfoHeader = cms.InputTag("generator"),
    lheInfo = cms.VInputTag(cms.InputTag("externalLHEProducer"), cms.InputTag("source")),
    preferredPDFs = cms.VPSet( # see https://lhapdf.hepforge.org/pdfsets.html
        cms.PSet( name = cms.string("NNPDF31_nnlo_as_0118_nf_4_mc_hessian"), lhaid = cms.uint32(306046) ),
        cms.PSet( name = cms.string("NNPDF31_nnlo_hessian_pdfas"), lhaid = cms.uint32(306000) ),
        # cms.PSet( name = cms.string("PDF4LHC15_nnlo_30_pdfas"), lhaid = cms.uint32(91400) ),
        # cms.PSet( name = cms.string("NNPDF30_nlo_as_0118"), lhaid = cms.uint32(260000) ), # for some 92X samples. Note that the nominal weight, 260000, is not included in the LHE ...
        # cms.PSet( name = cms.string("NNPDF30_lo_as_0130"), lhaid = cms.uint32(262000) ), # some MLM 80X samples have only this (e.g. /store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root )
        # cms.PSet( name = cms.string("NNPDF30_nlo_nf_4_pdfas"), lhaid = cms.uint32(292000) ), # some FXFX 80X samples have only this (e.g. WWTo1L1Nu2Q, WWTo4Q)
        # cms.PSet( name = cms.string("NNPDF30_nlo_nf_5_pdfas"), lhaid = cms.uint32(292200) ), # some FXFX 80X samples have only this (e.g. DYJetsToLL_Pt, WJetsToLNu_Pt, DYJetsToNuNu_Pt)
    ),
    namedWeightIDs = cms.vstring(),
    namedWeightLabels = cms.vstring(),
    lheWeightPrecision = cms.int32(14),
    maxPdfWeights = cms.uint32(150), 
    debug = cms.untracked.bool(False),
)

process.lheInfoTable = cms.EDProducer("LHETablesProducer",
    lheInfo = cms.VInputTag(cms.InputTag("externalLHEProducer"), cms.InputTag("source")),
    precision = cms.int32(14),
    storeLHEParticles = cms.bool(True) 
)


process.finalGenParticles = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("prunedGenParticles"),
    select = cms.vstring(
    "drop *",
        "keep++ abs(pdgId) == 15 & (pt > 15 ||  isPromptDecayed() )",#  keep full tau decay chain for some taus
    #"drop status==1 & pt < 1", #drop soft stable particle in tau decay
        "keep+ abs(pdgId) == 15 ",  #  keep first gen decay product for all tau
        "+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())", # keep gamma above 10 GeV (or all prompt) and its first parent
    "+keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15", #keep leptons, with at most one mother back in the history
    "drop abs(pdgId)= 2212 && abs(pz) > 1000", #drop LHC protons accidentally added by previous keeps
        "keep (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)", #keep all B and C hadrons
        "keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16",   # keep neutrinos
    "keep status == 3 || (status > 20 && status < 30)", #keep matrix element summary
        "keep isHardProcess() ||  fromHardProcessDecayed()  || fromHardProcessFinalState() || (statusFlags().fromHardProcess() && statusFlags().isLastCopy())",  #keep event summary based on status flags
    "keep  (status > 70 && status < 80 && pt > 15) ", # keep high pt partons right before hadronization
        "keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 37 ",   # keep VIP(articles)s
        #"keep abs(pdgId) == 310 && abs(eta) < 2.5 && pt > 1 ",                                                     # keep K0
        "keep (1000001 <= abs(pdgId) <= 1000039 ) || ( 2000001 <= abs(pdgId) <= 2000015)", #keep SUSY fiction particles
        
   )
)



##################### Tables for final output and docs ##########################
process.genParticleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("finalGenParticles"),
    cut = cms.string(""), #we should not filter after pruning
    name= cms.string("GenPart"),
    doc = cms.string("interesting gen particles "),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the taus
    variables = cms.PSet(
         pt  = Var("pt",  float, precision=8),
         phi = Var("phi", float,precision=8),
         eta  = Var("eta",  float,precision=8),
         mass = Var("?mass>10 || (pdgId==22 && mass > 1) || abs(pdgId)==24 || pdgId==23 || abs(pdgId)>1000000?mass:0", float,precision="?(abs(pdgId)==6 && statusFlags().isLastCopy())?20:8",doc="Mass stored for all particles with mass > 10 GeV and photons with mass > 1 GeV, plus W/Z and BSM particles. For other particles you can lookup from PDGID"),
         pdgId  = Var("pdgId", int, doc="PDG id"),
         status  = Var("status", int, doc="Particle status. 1=stable"),
         genPartIdxMother = Var("?numberOfMothers>0?motherRef(0).key():-1", int, doc="index of the mother particle"),
         statusFlags = (Var(
            "statusFlags().isLastCopyBeforeFSR()                  * 16384 +"
            "statusFlags().isLastCopy()                           * 8192  +"
            "statusFlags().isFirstCopy()                          * 4096  +"
            "statusFlags().fromHardProcessBeforeFSR()             * 2048  +"
            "statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +"
            "statusFlags().isHardProcessTauDecayProduct()         * 512   +"
            "statusFlags().fromHardProcess()                      * 256   +"
            "statusFlags().isHardProcess()                        * 128   +"
            "statusFlags().isDirectHadronDecayProduct()           * 64    +"
            "statusFlags().isDirectPromptTauDecayProduct()        * 32    +"
            "statusFlags().isDirectTauDecayProduct()              * 16    +"
            "statusFlags().isPromptTauDecayProduct()              * 8     +"
            "statusFlags().isTauDecayProduct()                    * 4     +"
            "statusFlags().isDecayedLeptonHadron()                * 2     +"
            "statusFlags().isPrompt()                             * 1      ",
            int, doc=("gen status flags stored bitwise, bits are: "
                "0 : isPrompt, "
                "1 : isDecayedLeptonHadron, "
                "2 : isTauDecayProduct, "
                "3 : isPromptTauDecayProduct, "
                "4 : isDirectTauDecayProduct, "
                "5 : isDirectPromptTauDecayProduct, "
                "6 : isDirectHadronDecayProduct, "
                "7 : isHardProcess, "
                "8 : fromHardProcess, "
                "9 : isHardProcessTauDecayProduct, "
                "10 : isDirectHardProcessTauDecayProduct, "
                "11 : fromHardProcessBeforeFSR, "
                "12 : isFirstCopy, "
                "13 : isLastCopy, "
                "14 : isLastCopyBeforeFSR, ")
            )),
    )
)


process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.NANOEDMAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:/ceph/dseith/singletop/eft/'+file_name),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

from particle_level_cfg import genp
genp(process)


print process.NANOAODSIMEventContent.outputCommands

process.nanoSequence = cms.Sequence(process.nanoMetadata * process.genWeightsTable * process.lheInfoTable * process.finalGenParticles * process.genParticleTable * process.mergedGenParticles * process.genParticles2HepMC * process.genParticles2HepMCHiggsVtx * process.particleLevel * process.rivetProducerHTXS * process.rivetLeptonTable * process.rivetMetTable * process.genJetTable * process.metMCTable)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v4', '')


process.nanoAOD_step = cms.Path(process.nanoSequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOEDMAODSIMoutput_step = cms.EndPath(process.NANOEDMAODSIMoutput)
process.schedule = cms.Schedule(process.nanoAOD_step, process.endjob_step, process.NANOEDMAODSIMoutput_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# process.options.numberOfThreads=cms.untracked.uint32(48)
# process.options.numberOfStreams=cms.untracked.uint32(0)
# process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)



# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
