# Alpaka-enabled version of pv_cfg.py with CLI parsing
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag
import sys, os, argparse

def root_path(f):
    if os.path.exists("/pnfs/psi.ch/cms/trivcat/store/user/werdmann" + f):
        return "root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/werdmann" + f
    elif f.startswith("/eos/cms/store"):
        return "root://cms-xrd-global.cern.ch//" + f[len("/eos/cms/"):] 
    else:
        return "root://xrootd-cms.infn.it/" + f

# CLI argument parser
parser = argparse.ArgumentParser(description="Alpaka PrimaryVertexProducer config")
parser.add_argument("--input", type=str, required=True, help="Comma-separated list of input ROOT files")
parser.add_argument("--nevent", type=int, default=10, help="Number of events to process")
parser.add_argument("--info", type=str, default="alpaka-test", help="Info string for analyzer")
args, unknown = parser.parse_known_args()

# Input files
files = args.input.split(",")
source_files = cms.untracked.vstring(*[root_path(f) for f in files])

process = cms.Process("RERECO", eras.Phase2)
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag = GlobalTag(process.GlobalTag, "141X_mcRun4_realistic_v3", "")
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')

process.source = cms.Source("PoolSource", fileNames=source_files)
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(args.nevent))
process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True)
)

process.Timing = cms.Service("Timing", summaryOnly=cms.untracked.bool(True))

# Alpaka Producers
process.tracksSoA = cms.EDProducer("PortableTrackSoAProducer@alpaka",
    TrackLabel=cms.InputTag("generalTracks"),
    BeamSpotLabel=cms.InputTag("offlineBeamSpot"),
    TkFilterParameters=cms.PSet(
        maxNormalizedChi2=cms.double(10.0),
        minPixelLayersWithHits=cms.int32(2),
        minSiliconLayersWithHits=cms.int32(5),
        maxD0Significance=cms.double(4.0),
        maxD0Error=cms.double(1.0),
        maxDzError=cms.double(1.0),
        minPt=cms.double(0.0),
        maxEta=cms.double(2.4),
        trackQuality=cms.string("any"),
        vertexSize=cms.double(0.006),
        d0CutOff=cms.double(3.0)
    )
)

process.beamSpotSoA = cms.EDProducer("PortableBeamSpotSoAProducer@alpaka",
    BeamSpotLabel=cms.InputTag("offlineBeamSpot")
)

process.vertexSoA = cms.EDProducer("PrimaryVertexProducer_Alpaka@alpaka",
    TrackLabel=cms.InputTag("tracksSoA"),
    BeamSpotLabel=cms.InputTag("beamSpotSoA"),
    blockOverlap=cms.double(0.5),
    blockSize=cms.int32(512),
    TkFitterParameters=cms.PSet(
        chi2cutoff=cms.double(2.5),
        minNdof=cms.double(0.0),
        useBeamSpotConstraint=cms.bool(True),
        maxDistanceToBeam=cms.double(1.0)
    ),
    TkClusParameters=cms.PSet(
        coolingFactor=cms.double(0.6),
        zrange=cms.double(4.0),
        delta_highT=cms.double(1.e-2),
        delta_lowT=cms.double(1.e-3),
        convergence_mode=cms.int32(0),
        Tmin=cms.double(2.0),
        Tpurge=cms.double(2.0),
        Tstop=cms.double(0.5),
        vertexSize=cms.double(0.006),
        d0CutOff=cms.double(3.0),
        dzCutOff=cms.double(3.0),
        zmerge=cms.double(1e-2),
        uniquetrkweight=cms.double(0.8),
        uniquetrkminp=cms.double(0.0)
    )
)

process.vertexAoS = cms.EDProducer("SoAToRecoVertexProducer",
    soaVertex=cms.InputTag("vertexSoA"),
    srcTrack=cms.InputTag("generalTracks")
)

# Analyzer
process.oldVertexAnalysis = cms.EDAnalyzer("TestAnalyzer",
    info=cms.untracked.string(args.info),
    f4D=cms.untracked.bool(False),
    frefit=cms.untracked.bool(False),
    selNdof=cms.double(4.0),
    selNdofWithBS=cms.double(2.0),
    beamSpot=cms.InputTag('offlineBeamSpot'),
    simG4=cms.InputTag("g4SimHits"),
    outputFile=cms.untracked.string("pv_alpaka.root"),
    verbose=cms.untracked.bool(True),
    veryverbose=cms.untracked.bool(False),
    recoTrackProducer=cms.untracked.string("generalTracks"),
    zmatch=cms.untracked.double(0.05),
    zMatchWindow_sigmas=cms.double(3.0),
    autodump=cms.untracked.int32(0),
    nDump=cms.untracked.int32(1),
    nDumpTracks=cms.untracked.int32(1),
    RECO=cms.untracked.bool(True),
    track_timing=cms.untracked.bool(False),
    TkFilterParameters=process.tracksSoA.TkFilterParameters,
    trackingParticleCollection=cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackingVertexCollection=cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackAssociatorMap=cms.untracked.InputTag("trackingParticleRecoTrackAsssociation"),
    vertexAssociator=cms.untracked.InputTag("VertexAssociatorByPositionAndTracks"),
    lumiInfoTag=cms.untracked.InputTag("LumiInfo", "brilcalc"),
    useVertexFilter=cms.untracked.bool(False),
    compareCollections=cms.untracked.int32(0),
    vertexRecoCollections=cms.VInputTag([
        cms.InputTag("vertexAoS")
    ])
)

process.analyze = cms.Path(process.oldVertexAnalysis)
process.vertexing = cms.Path(
    process.tracksSoA + process.beamSpotSoA + process.vertexSoA + process.vertexAoS
)
process.schedule = cms.Schedule(process.vertexing, process.analyze)
