import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *


def genp(process):
    process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
        inputPruned = cms.InputTag("prunedGenParticles"),
        inputPacked = cms.InputTag("packedGenParticles"),
        # inputPruned = cms.InputTag("genParticles"),
        # inputPacked = cms.InputTag("genParticles"),
    )

    process.genParticles2HepMC = cms.EDProducer("GenParticles2HepMCConverter",
        genParticles = cms.InputTag("mergedGenParticles"),
        genEventInfo = cms.InputTag("generator"),
        signalParticlePdgIds = cms.vint32(),
    )

    process.genParticles2HepMCHiggsVtx = cms.EDProducer("GenParticles2HepMCConverter",
         genParticles = cms.InputTag("mergedGenParticles"),
         genEventInfo = cms.InputTag("generator"),
         signalParticlePdgIds = cms.vint32(25), ## for the Higgs analysis
    )


    process.particleLevel = cms.EDProducer("ParticleLevelProducer",
        src = cms.InputTag("genParticles2HepMC:unsmeared"),
        
        usePromptFinalStates = cms.bool(True), # for leptons, photons, neutrinos
        excludePromptLeptonsFromJetClustering = cms.bool(False),
        excludeNeutrinosFromJetClustering = cms.bool(True),
        
        particleMinPt  = cms.double(0.),
        particleMaxEta = cms.double(5.), # HF range. Maximum 6.0 on MiniAOD
        
        lepConeSize = cms.double(0.1), # for photon dressing
        lepMinPt    = cms.double(15.),
        lepMaxEta   = cms.double(2.5),
        
        jetConeSize = cms.double(0.4),
        jetMinPt    = cms.double(10.),
        jetMaxEta   = cms.double(999.),
        
        fatJetConeSize = cms.double(0.8),
        fatJetMinPt    = cms.double(170.),
        fatJetMaxEta   = cms.double(999.),
    )

    process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
       HepMCCollection = cms.InputTag('genParticles2HepMCHiggsVtx','unsmeared'),
       LHERunInfo = cms.InputTag('externalLHEProducer'),
       ProductionMode = cms.string('AUTO'),
    )


    # ##################### Tables for final output and docs ##########################
    process.rivetLeptonTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("particleLevel:leptons"),
        cut = cms.string(""),
        name= cms.string("GenDressedLepton"),
        doc = cms.string("Dressed leptons from Rivet-based ParticleLevelProducer"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        # externalVariables = cms.PSet(
        #     hasTauAnc = ExtVar(cms.InputTag("tautagger"),bool, doc="true if Dressed lepton has a tau as ancestor"),
        #     ),
        variables = cms.PSet(
            P4Vars,
            pdgId = Var("pdgId", int, doc="PDG id"), 
        )
    )

    process.rivetMetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("particleLevel:mets"),
        name = cms.string("MET"),
        doc = cms.string("MET from Rivet-based ParticleLevelProducer in fiducial volume abs(eta)<5"),
        singleton = cms.bool(True),  # there's always exactly one MET per event
        extension = cms.bool(True), # this is the main table
        variables = cms.PSet(
           fiducialGenPt  = Var("pt",  float, precision=10),
           fiducialGenPhi = Var("phi", float, precision=10),
        ),
    )

    process.genJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("slimmedGenJets"),
    cut = cms.string("pt > 10"),
    name = cms.string("GenJet"),
    doc  = cms.string("slimmedGenJets, i.e. ak4 Jets made with visible genparticles"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the genjets
    variables = cms.PSet(P4Vars,
	#anything else?
    )
    )

    process.metTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("slimmedMETs"),
    name = cms.string("MET"),
    doc = cms.string("slimmedMET, type-1 corrected PF MET"),
    singleton = cms.bool(True),  # there's always exactly one MET per event
    extension = cms.bool(False), # this is the main table for the MET
    variables = cms.PSet(PTVars,
       sumEt = Var("sumEt()", float, doc="scalar sum of Et",precision=10),
       covXX = Var("getSignificanceMatrix().At(0,0)",float,doc="xx element of met covariance matrix", precision=8),
       covXY = Var("getSignificanceMatrix().At(0,1)",float,doc="xy element of met covariance matrix", precision=8),
       covYY = Var("getSignificanceMatrix().At(1,1)",float,doc="yy element of met covariance matrix", precision=8),
       significance = Var("metSignificance()", float, doc="MET significance",precision=10),
       MetUnclustEnUpDeltaX = Var("shiftedPx('UnclusteredEnUp')-px()", float, doc="Delta (METx_mod-METx) Unclustered Energy Up",precision=10),
       MetUnclustEnUpDeltaY = Var("shiftedPy('UnclusteredEnUp')-py()", float, doc="Delta (METy_mod-METy) Unclustered Energy Up",precision=10),

    ),
)



    process.metMCTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = process.metTable.src,
    name = cms.string("GenMET"),
    doc = cms.string("Gen MET"),
    singleton = cms.bool(True),  
    extension = cms.bool(False),
    variables = cms.PSet(
       pt  = Var("genMET.pt",  float, doc="pt", precision=10),
       phi = Var("genMET.phi", float, doc="phi", precision=10),
    ),
)

#     process.patJetPartons = cms.EDProducer('HadronAndPartonSelector',
#     src = cms.InputTag("generator"),
#     particles = cms.InputTag("prunedGenParticles"),
#     partonMode = cms.string("Auto"),
#     fullChainPhysPartons = cms.bool(True)
# )
# 
#     process.genJetFlavourAssociation = cms.EDProducer("JetFlavourClustering",
#     jets = process.genJetTable.src,
#     bHadrons = cms.InputTag("patJetPartons","bHadrons"),
#     cHadrons = cms.InputTag("patJetPartons","cHadrons"),
#     partons = cms.InputTag("patJetPartons","physicsPartons"),
#     leptons = cms.InputTag("patJetPartons","leptons"),
#     jetAlgorithm = cms.string("AntiKt"),
#     rParam = cms.double(0.4),
#     ghostRescaling = cms.double(1e-18),
#     hadronFlavourHasPriority = cms.bool(False)
# )
#     process.genJetFlavourTable = cms.EDProducer("GenJetFlavourTableProducer",
#     name = process.genJetTable.name,
#     src = process.genJetTable.src,
#     cut = process.genJetTable.cut,
#     deltaR = cms.double(0.1),
#     jetFlavourInfos = cms.InputTag("slimmedGenJetsFlavourInfos"),
# )
