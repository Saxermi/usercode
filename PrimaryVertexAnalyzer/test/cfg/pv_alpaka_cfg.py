import FWCore.ParameterSet.Config as cms

# from Configuration.Eras.Era_Run3_cff import Run3

from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag

import sys, os
import subprocess, shlex


def root_path(f):
    # f is a file path as e.g. das would spit it out:  /store/...../XXX.root
    # prefers the local version if there is one
    if os.path.exists("/pnfs/psi.ch/cms/trivcat/store/user/werdmann" + f):
        return (
            "root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/werdmann"
            + f
        )
    elif f.startswith("/eos/cms/store"):
        return "root://cms-xrd-global.cern.ch//" + f[len("/eos/cms/") :]
    else:
        return "root://xrootd-cms.infn.it/" + f


timealgo = ""  # i.e. none

parameters = {  # can be overwritten by arguments of the same name
    "4D": cms.untracked.bool(False),
    "refit": cms.untracked.bool(False),
    "zMatchWindow_sigmas": cms.untracked.double(3.0),
    "compareCollections": cms.untracked.int32(0),
    "selNdof": cms.untracked.double(4.0),
    "selNdofWithBS": cms.untracked.double(2.0),
    "usefit": cms.untracked.bool(False),
    "use2file": cms.untracked.bool(False),
    "fill_track_histos": cms.untracked.bool(False),
    "fplo": cms.untracked.double(0),
    "chi2cutoff": cms.double(4.0),
    "verboseAnalyzer": cms.untracked.bool(True),
    "verboseProducer": cms.untracked.bool(False),
    "verboseClusterizer": cms.untracked.bool(False),
    "verboseClusterizer2D": cms.untracked.bool(False),
    "zdumpcenter": cms.untracked.double(0.0),
    "zdumpwidth": cms.untracked.double(20.0),
    "reco": cms.untracked.bool(False),
    "miniaod": cms.untracked.bool(False),
    "autodump": cms.untracked.int32(0),
    "nDump": cms.untracked.int32(1),
    "nDumpTracks": cms.untracked.int32(1),
    "uniquetrkweight": cms.double(0.8),
    "uniquetrkminp": cms.double(0.0),
    "zmerge": cms.double(1.0e-2),
    "coolingFactor": cms.double(0.6),
    "Tmin": cms.double(2.0),
    "Tpurge": cms.double(2.0),
    "Tstop": cms.double(0.5),
    "vertexSize": cms.double(0.006),
    "vertexSizeTime": cms.double(0.008),
    "d0CutOff": cms.double(3.0),
    "dzCutOff": cms.double(3.0),
    "zrange": cms.double(4),
    "convergence_mode": cms.int32(0),
    "delta_lowT": cms.double(1.0e-3),
    "delta_highT": cms.double(1.0e-2),
    "minValidStripHits": cms.int32(0),
    # track selection
    "maxNormalizedChi2": cms.double(10.0),
    "minPixelLayersWithHits": cms.int32(2),
    "minSiliconLayersWithHits": cms.int32(5),
    "maxD0Significance": cms.double(4.0),
    "minPt": cms.double(0.0),
    "maxEta": cms.double(2.4),
    "trackQuality": cms.string("any"),
    # track selection, experimental
    "maxDzError": cms.double(1.0),
    "maxD0Error": cms.double(1.0),
    # vertex selection
    "minNdof": cms.double(2.0),
    # clustering in blocks
    "runInBlocks": cms.bool(True),
    "block_size": cms.uint32(1024),
    "overlap_frac": cms.double(0.5),
}
info = "alpaka_run4_test"


process = cms.Process("RERECO", eras.Phase2)
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, "141X_mcRun4_realistic_v3", "")
parameters["maxEta"] = cms.double(4.0)


# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

# process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
# process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
##

process.load("Configuration.StandardSequences.RawToDigi_cff")  #
process.load("Configuration.StandardSequences.L1Reco_cff")  #
process.load("CommonTools.ParticleFlow.EITopPAG_cff")  #
process.load("Configuration.StandardSequences.AlCaRecoStreamsMC_cff")
process.load("Configuration.StandardSequences.Validation_cff")
process.load("DQMOffline.Configuration.DQMOfflineMC_cff")

process.load("Configuration.StandardSequences.Accelerators_cff")
process.load("HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi")


# Input files
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
        root_path(
            "/store/relval/CMSSW_15_1_0_pre1/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v1/2580000/c699f8a2-4146-4edb-8614-9bf8d188ffee.root"
        )
    ),
    secondaryFileNames=cms.untracked.vstring(),
)

# use tp
# process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi") # should i add this also to the src files? but it should be there already
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

print(process)
print("*************test********************")
print(process.tpClusterProducer)
process.theTruth = cms.Sequence(
    process.tpClusterProducer
    * process.quickTrackAssociatorByHits
    * process.trackingParticleRecoTrackAsssociation
)

process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True),
    #  TryToContinue=cms.untracked.vstring(
    #    "ProductNotFound"
    # ),  # temporary trying to fix the beamspot issue
)


process.Timing = cms.Service("Timing", summaryOnly=cms.untracked.bool(True))

process.load("RecoVertex.Configuration.RecoVertex_cff")

process.ak4CaloJetsForTrk.srcPVs = "unsortedOfflinePrimaryVertices"
process.vertexreco.remove(process.generalV0Candidates)
process.vertexreco.remove(process.inclusiveVertexing)


# Number of events to run
process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(1),
)

# Production metadata
process.configurationMetadata = cms.untracked.PSet(
    annotation=cms.untracked.string("PV nevts:10"),
    name=cms.untracked.string("Applications"),
    version=cms.untracked.string("$Revision: 1.19 $"),
)


# Output definition
process.FEVToutput = cms.OutputModule(
    "PoolOutputModule",
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string("GEN-SIM-DIGI-RECO"),
        filterName=cms.untracked.string(""),
    ),
    fileName=cms.untracked.string("testAlpaka_PU200.root"),  # output file name
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep *_tracksSoA_*_*",
        "keep *_beamSpotDevice_*_*",
        "keep *_vertexSoA_*_*",
        "keep *_vertexAoS_*_*",
    ),  # I.e., just drop everything and keep things in this module
    splitLevel=cms.untracked.int32(0),
)

# Endpath and output
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVToutput_step = cms.EndPath(process.FEVToutput)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete

process = customiseEarlyDelete(process)

import sys
import argparse

parser = argparse.ArgumentParser(
    prog=f"{sys.argv[0]} {sys.argv[1]} --",
    description="Test and validation of PrimaryVertexProducer_Alpaka",
)
parser.add_argument(
    "-b",
    "--backend",
    type=str,
    default="auto",
    help="Alpaka backend. Comma separated list. Possible options: cpu, gpu-nvidia, gpu-amd",
)
args = parser.parse_args()

# Set the backend for all jobs
process.options.accelerators = args.backend.split(",")
print(">>> Alpaka backend(s):", process.options.accelerators)

################################
## Now the plugins themselves ##
################################

# Convert reco::Track to portable Track
process.tracksSoA = cms.EDProducer(
    "PortableTrackSoAProducer@alpaka",
    TrackLabel=cms.InputTag("generalTracks"),
    BeamSpotLabel=cms.InputTag("offlineBeamSpot"),
    # For historical reasons, parameters are set up as in the legacy code of RecoVertex/PrimaryVertexProducer/python/OfflinePrimaryVertices_cfi.py
    TkFilterParameters=cms.PSet(
        maxNormalizedChi2=cms.double(10.0),
        minPixelLayersWithHits=cms.int32(2),
        minSiliconLayersWithHits=cms.int32(5),
        maxD0Significance=cms.double(4.0),
        maxD0Error=cms.double(1.0),
        maxDzError=cms.double(1.0),
        minPt=cms.double(0.0),
        maxEta=cms.double(4.0),  # cms.double(2.4),
        trackQuality=cms.string("any"),
        vertexSize=cms.double(0.006),
        d0CutOff=cms.double(3.0),
    ),
)

#
#
#
# Convert reco::BeamSpot to portable BeamSpot
# process.beamSpotSoA = cms.EDProducer(
#   "PortableBeamSpotSoAProducer@alpaka", BeamSpotLabel=cms.InputTag("offlineBeamSpot")
# )
process.beamSpotDevice = cms.EDProducer(
    "BeamSpotDeviceProducer@alpaka",
    src=cms.InputTag("offlineBeamSpot"),
    alpaka=cms.untracked.PSet(backend=cms.untracked.string("")),
)


process.vertexSoA = cms.EDProducer(
    "PrimaryVertexProducer_Alpaka@alpaka",
    TrackLabel=cms.InputTag("tracksSoA"),
    # BeamSpotLabel=cms.InputTag("beamSpotSoA"),
    BeamSpotLabel=cms.InputTag("beamSpotDevice"),
    blockOverlap=cms.double(0.50),
    blockSize=cms.int32(512),
    TkFitterParameters=cms.PSet(
        chi2cutoff=cms.double(2.5),
        minNdof=cms.double(0.0),
        useBeamSpotConstraint=cms.bool(True),
        maxDistanceToBeam=cms.double(1.0),
    ),
    TkClusParameters=cms.PSet(
        coolingFactor=cms.double(0.6),  # moderate annealing speed
        zrange=cms.double(
            4.0
        ),  # consider only clusters within 4 sigma*sqrt(T) of a track
        delta_highT=cms.double(1.0e-2),  # convergence requirement at high T
        delta_lowT=cms.double(1.0e-3),  # convergence requirement at low T
        convergence_mode=cms.int32(0),  # 0 = two steps, 1 = dynamic with sqrt(T)
        Tmin=cms.double(2.0),  # end of vertex splitting
        Tpurge=cms.double(2.0),  # cleaning
        Tstop=cms.double(0.5),  # end of annealing
        vertexSize=cms.double(0.006),  # added in quadrature to track-z resolutions
        d0CutOff=cms.double(3.0),  # downweight high IP tracks
        dzCutOff=cms.double(3.0),  # outlier rejection after freeze-out (T<Tmin)
        zmerge=cms.double(
            1e-2
        ),  # merge intermediat clusters separated by less than zmerge
        uniquetrkweight=cms.double(
            0.8
        ),  # require at least two tracks with this weight at T=Tpurge
        uniquetrkminp=cms.double(
            0.0
        ),  # minimal a priori track weight for counting unique tracks
    ),
)

process.vertexAoS = cms.EDProducer(
    "SoAToRecoVertexProducer",
    soaVertex=cms.InputTag("vertexSoA"),
    srcTrack=cms.InputTag("generalTracks"),
)

# modify the track filter if needed
tkFilterParameters = process.unsortedOfflinePrimaryVertices.TkFilterParameters.clone()

print("original trackFilterParameters (z-clustering)")
print(tkFilterParameters)

for par_name in (
    "maxNormalizedChi2",
    "minPixelLayersWithHits",
    "minSiliconLayersWithHits",
    "maxD0Significance",
    "maxD0Error",
    "maxDzError",
    "minPt",
    "maxEta",
    "trackQuality",
):
    try:
        default = getattr(tkFilterParameters, par_name)
        if default != parameters[par_name]:
            print(
                f"pv_cfg: changing tkFilter parameter '{par_name}' from '{default}' to '{parameters[par_name]}'"
            )
        setattr(tkFilterParameters, par_name, parameters[par_name])
    except ValueError:
        print(f"pv_cfg: attribute 'tkFilterParameters.{par_name}' not found ")


# analysis
# process.oldVertexAnalysis = cms.EDAnalyzer("PrimaryVertexAnalyzer4PU",
vcollections = ["offlinePrimaryVertices", "vertexAoS"]

process.oldVertexAnalysis = cms.EDAnalyzer(
    "TestAnalyzer",
    info=cms.untracked.string(info),
    f4D=cms.untracked.bool(False),
    frefit=cms.untracked.bool(False),
    selNdof=parameters["selNdof"],
    selNdofWithBS=parameters["selNdofWithBS"],
    beamSpot=cms.InputTag("offlineBeamSpot"),
    simG4=cms.InputTag("g4SimHits"),
    outputFile=cms.untracked.string("pv.root"),
    verbose=parameters["verboseAnalyzer"],
    veryverbose=cms.untracked.bool(False),
    recoTrackProducer=cms.untracked.string("generalTracks"),
    zmatch=cms.untracked.double(0.05),
    zMatchWindow_sigmas=parameters["zMatchWindow_sigmas"],
    autodump=parameters["autodump"],
    nDump=parameters["nDump"],
    nDumpTracks=parameters["nDumpTracks"],
    RECO=cms.untracked.bool(True),
    track_timing=cms.untracked.bool(False),
    TkFilterParameters=tkFilterParameters,
    trackingParticleCollection=cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackingVertexCollection=cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackAssociatorMap=cms.untracked.InputTag("trackingParticleRecoTrackAsssociation"),
    vertexAssociator=cms.untracked.InputTag("VertexAssociatorByPositionAndTracks"),
    lumiInfoTag=cms.untracked.InputTag("LumiInfo", "brilcalc"),
    useVertexFilter=cms.untracked.bool(False),
    compareCollections=parameters["compareCollections"],
    vertexRecoCollections=cms.VInputTag(vcollections),
)
process.analyze = cms.Path(process.theTruth * process.oldVertexAnalysis)


###################################
## Last, organize paths and exec ##
###################################

process.vertexing_task = cms.EndPath(
    process.tracksSoA + process.beamSpotDevice + process.vertexSoA + process.vertexAoS
)
# process.vertexing_task = cms.EndPath(process.beamSpotDevice)

process.content = cms.EDAnalyzer("EventContentAnalyzer")
# process.dump = cms.Path(process.content)

# .vertexing_task = cms.EndPath(
#   process.tracksSoA + process.beamSpotSoA + process.vertexSoA + process.vertexAoS
# )

process.schedule = cms.Schedule(
    process.vertexing_task,
    # process.dump,
    process.analyze,
    process.endjob_step,
    # process.FEVToutput_step
)
