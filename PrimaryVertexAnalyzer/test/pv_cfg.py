import FWCore.ParameterSet.Config as cms
#import FWCore.PythonUtilities.LumiList as LumiList
from Configuration.StandardSequences.Eras import eras
from Configuration.AlCa.GlobalTag import GlobalTag

import sys, os
import subprocess, shlex


def root_path(f):
    # f is a file path as e.g. das would spit it out:  /store/...../XXX.root
    # prefers the local version if there is one
    if os.path.exists( "/pnfs/psi.ch/cms/trivcat/store/user/werdmann" + f):
        return "root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/werdmann" + f
    elif f.startswith("/eos/cms/store"):
        return "root://cms-xrd-global.cern.ch//"+ f[len("/eos/cms/"):]
    else:
        return "root://xrootd-cms.infn.it/" + f



DO_VTX_RECO = True
args =""
nevent=-1
info = "test"
debug = False
eventsToProcess = "" #cms.untracked.VEventRange('1:4474')
zdumpcenter = 0
zdumpwidth = 100.
source_files = ()



timealgo = ""  # i.e. none

parameters={  # can be overwritten by arguments of the same name
  "4D": cms.untracked.bool(False),
  "refit": cms.untracked.bool(False),
  "zMatchWindow_sigmas": cms.untracked.double(3.0),  
  "compareCollections" : cms.untracked.int32(0),
  "selNdof": cms.untracked.double(4.0),
  "selNdofWithBS": cms.untracked.double(2.0),
  "usefit":cms.untracked.bool(False),
  "use2file":cms.untracked.bool(False),
  "fill_track_histos":cms.untracked.bool(False),
  "fplo":cms.untracked.double(0),
  "chi2cutoff":cms.double(4.), 
  "verboseAnalyzer":cms.untracked.bool(False),
  "verboseProducer":cms.untracked.bool(False),
  "verboseClusterizer":cms.untracked.bool(False),
  "verboseClusterizer2D":cms.untracked.bool(False),
  "zdumpcenter" : cms.untracked.double(0.0),
  "zdumpwidth" : cms.untracked.double(20.0),
  "reco":cms.untracked.bool(False),
  "miniaod":cms.untracked.bool(False),
  "autodump":cms.untracked.int32(0),
  "nDump":cms.untracked.int32(1),
  "nDumpTracks":cms.untracked.int32(1),
  "uniquetrkweight":cms.double(0.8),
  "uniquetrkminp":cms.double(0.0),
  "zmerge":cms.double(1.e-2),
  "coolingFactor":cms.double(0.6),
  "Tmin": cms.double(2.0),
  "Tpurge":cms.double(2.0),
  "Tstop": cms.double(0.5),
  "vertexSize" : cms.double(0.006),
  "vertexSizeTime" : cms.double(0.008),
  "d0CutOff" : cms.double(3.),
  "dzCutOff" : cms.double(3.),
  "zrange" : cms.double(4),
  "convergence_mode" : cms.int32(0),
  "delta_lowT" : cms.double(1.e-3),
  "delta_highT" : cms.double(1.e-2),
# track selection
  "maxNormalizedChi2":cms.double(10.0),
  "minPixelLayersWithHits":cms.int32(2),
  "minSiliconLayersWithHits":cms.int32(5),
  "maxD0Significance":cms.double(4.0), 
  "minPt":cms.double(0.0),
  "maxEta":cms.double(2.4),
  "trackQuality":cms.string("any"),
# track selection, experimental
  "maxDzError":cms.double(1.0),
  "maxD0Error":cms.double(1.0),
# vertex selection
  "minNdof": cms.double( 2.0 ),
# clustering in blocks
  "runInBlocks" : cms.bool(True),# default value run in blocks is true
  "block_size" : cms.uint32(256),
    "overlap_frac" : cms.double(0.0)
}

print(">>>>>>>>>>>>>>>>>>>>>>   this is pv_cfg.py   <<<<<<<<<<<<<<<<<<<<<<<<<")

# temporary fix, should not be needed
args=[]
for a in sys.argv:
    args += a.split()
        

for a in args:

    try:
        key, value = a.split("=")
    except ValueError:
        continue


    if key == "nevent" or key == "events":

        nevent = int(value)

    elif key == "input":

        print( "processing input directive, value=",value)
        inputfile = value
        if inputfile.endswith(".root"):
            # a single root file
            #files = [inputfile]
            #source_files =  cms.untracked.vstring( root_path(inputfile) )
            # a comma separated list of root files
            files = inputfile.split(",")
            source_files = cms.untracked.vstring(*tuple([root_path(f) for f in files]))
        else:
            # text file containing a list of filenames
            files = [l.strip() for l in open(inputfile).readlines()]
            source_files = cms.untracked.vstring(*tuple([root_path(f) for f in files]))


    elif key == "redopv" or key == "redoPV":
        DO_VTX_RECO = (value in ("True","yes","Yes") )

    elif key == "zdump":
        debug = True
        try:
            zdumpcenter,zdumpwidth = value.split(":")
            zdumpcenter = float(zdumpcenter)
            zdumpwidth = float(zdumpwidth)
            parameters["zdumpcenter"] = cms.untracked.double(zdumpcenter)
            parameters["zdumpwidth"] = cms.untracked.double(zdumpwidth)
            nevent = 1
        except ValueError:
            nevent = int(value)

    elif key == "debug":
        debug = True
        try:
            run,event,zdumpcenter,zdumpwidth = value.split(":")
            eventsToProcess = '%s:%s'%(run,event)
            zdumpcenter = float(zdumpcenter)
            zdumpwidth = float(zdumpwidth)
            parameters["zdumpcenter"] = cms.untracked.double(zdumpcenter)
            parameters["zdumpwidth"] = cms.untracked.double(zdumpwidth)
            nevent = 1
        except ValueError:
            nevent = int(value)
        
    elif key == "verbose":
        if (value in ("True","yes","Yes") ):
            parameters["verboseClusterizer"] = cms.untracked.bool(True)
            parameters["verboseAnalyzer"] = cms.untracked.bool(True)
            parameters["verboseProducer"] = cms.untracked.bool(True)
    elif key == "info":
        info  = value
    elif key in parameters.keys():
        typename = parameters[key].configTypeName()
        print( "setting parameter ", key, value, type(parameters[key]),"'%s'"%typename)
        if typename == "double":
            parameters[key] = cms.double( float(value) )
        elif type(parameters[key])==type(cms.string("")):
            parameters[key] = cms.string( value )
        elif typename == "untracked double":
            parameters[key] = cms.untracked.double( float(value) )
        elif typename == "int32":
            parameters[key] = cms.int32( int(value) )
        elif typename == "untracked int32":
            parameters[key] = cms.untracked.int32( int(value) )
        elif typename == "bool":
            parameters[key] = cms.bool( value == "True")
        elif typename == "untracked bool":
            parameters[key] = cms.untracked.bool( value == "True")
    else:
        print( "!! pv_cfg.py  :  unknown key ",key)



print( "pv_cfg.py")
print( "events         ",nevent)
print( "info           ",info)
print( "DO_VTX_RECO    ",DO_VTX_RECO)

if len(source_files) == 0:
    print( "empty source file list!")
else:
    print( "source files = ",source_files)
# 



# done processing the options, now construct the CMSSW process

#>>>  input file dependent !!
process = cms.Process("RERECO", eras.Phase2)
process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag = GlobalTag(process.GlobalTag, "141X_mcRun4_realistic_v1", '')
parameters["maxEta"] = cms.double(4.0)
#<<<  


# load standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000   # don't be too noisy
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreamsMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')




# use tp
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.theTruth = cms.Sequence(
    process.tpClusterProducer *
    process.quickTrackAssociatorByHits *
    process.trackingParticleRecoTrackAsssociation
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(nevent))


# Input source
process.source = cms.Source("PoolSource",
                            fileNames = source_files,
)

if not eventsToProcess == "":
    process.source.eventsToProcess = cms.untracked.VEventRange( eventsToProcess )



process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)





process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True)  )

process.load("RecoVertex.Configuration.RecoVertex_cff")

process.ak4CaloJetsForTrk.srcPVs = 'unsortedOfflinePrimaryVertices'
process.vertexreco.remove(process.generalV0Candidates)
process.vertexreco.remove(process.inclusiveVertexing)

# modify the track filter if needed
tkFilterParameters = process.unsortedOfflinePrimaryVertices.TkFilterParameters.clone()    

print( "original trackFilterParameters (z-clustering)")
print( tkFilterParameters)

for par_name in ("maxNormalizedChi2", "minPixelLayersWithHits", 
                 "minSiliconLayersWithHits", "maxD0Significance",
                 "maxD0Error", "maxDzError", "minPt", "maxEta", "trackQuality"):
    try:
        default =  getattr(tkFilterParameters, par_name)
        if default != parameters[par_name]:
            print( f"pv_cfg: changing tkFilter parameter '{par_name}' from '{default}' to '{parameters[par_name]}'")
        setattr(tkFilterParameters, par_name, parameters[par_name])
    except ValueError:
        print( f"pv_cfg: attribute 'tkFilterParameters.{par_name}' not found ")



# modify the clustering parameters if needed
tkClusParameters = process.unsortedOfflinePrimaryVertices.TkClusParameters.clone()    
print( "original z clustering parameters")
print( tkClusParameters)
tkDAClusParameters = tkClusParameters.TkDAClusParameters
for par_name in ("zrange", "convergence_mode", "delta_lowT", 
                 "Tmin", "Tpurge", "Tstop", "coolingFactor",
                 "uniquetrkminp", "uniquetrkweight",
                 "delta_highT", "vertexSize", "zmerge",
                 "runInBlocks", "block_size", "overlap_frac"):
    try:
        default =  getattr(tkDAClusParameters, par_name)
        if default != parameters[par_name]:
            print( "changing TkDAClusParameters parameter ", par_name, " from ", default, "  to ", parameters[par_name])
            setattr(tkDAClusParameters, par_name, parameters[par_name])
    except AttributeError:
        print( "pv_cfg: attribute TkDAClusParameters.", par_name , " not found ")

    setattr(tkDAClusParameters, "zdumpcenter", parameters[ "zdumpcenter"])  
    setattr(tkDAClusParameters, "zdumpwidth", parameters[ "zdumpwidth"])


print( "modified z clustering parameters")
print( tkClusParameters)


# configure the vertex producer to create new vertex collection called testVertices:test
test = cms.PSet(label=cms.string("test"),
                   algorithm=cms.string("AdaptiveVertexFitter"),
                   chi2cutoff = cms.double(2.5), 
                   minNdof=cms.double(0.0),
                   useBeamConstraint = cms.bool(False),
                   maxDistanceToBeam = cms.double(1.0)
              )


process.testVertices = cms.EDProducer(
    "PrimaryVertexProducer",

    verbose = parameters["verboseProducer"],
    TrackLabel = cms.InputTag("generalTracks"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    
    TkFilterParameters = tkFilterParameters,

    TkClusParameters = tkClusParameters,

    vertexCollections = cms.VPSet([test]),

    isRecoveryIteration = cms.bool(False),
    recoveryVtxCollection = cms.InputTag("")
)
process.MPV = cms.Path(process.testVertices)


# configure the analyzer:
# "offlinePrimaryVertices" are the vertices that are already present in the input file, ie. produced with the default configuration
# "testVertices:test"  are the vertices produced by re-running the producer with the modified configuration, e.g. runInBlocks=True

if DO_VTX_RECO:
    vcollections = ["offlinePrimaryVertices", "testVertices:test"]
else:
    vcollections = ["offlinePrimaryVertices"]  # alternative "offlinePrimaryVerticesWithBS" 


# analysis
#process.oldVertexAnalysis = cms.EDAnalyzer("PrimaryVertexAnalyzer4PU",
process.oldVertexAnalysis = cms.EDAnalyzer("TestAnalyzer",
    info=  cms.untracked.string(info),
    f4D = cms.untracked.bool(False),
    frefit = cms.untracked.bool(False),
    selNdof = parameters["selNdof"],                                    
    selNdofWithBS = parameters["selNdofWithBS"],                                    
    beamSpot = cms.InputTag('offlineBeamSpot'),
    simG4 = cms.InputTag("g4SimHits"),
    outputFile = cms.untracked.string("pv.root"),
    verbose = parameters["verboseAnalyzer"],
    veryverbose = cms.untracked.bool(False),
    recoTrackProducer = cms.untracked.string("generalTracks"),
    zmatch = cms.untracked.double(0.05),
    zMatchWindow_sigmas = parameters["zMatchWindow_sigmas"],
    autodump = parameters["autodump"],
    nDump = parameters["nDump"],
    nDumpTracks = parameters["nDumpTracks"],
    RECO = cms.untracked.bool(True),
    track_timing = cms.untracked.bool(False),
    TkFilterParameters = tkFilterParameters,
    trackingParticleCollection = cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackingVertexCollection = cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackAssociatorMap = cms.untracked.InputTag("trackingParticleRecoTrackAsssociation"),
    vertexAssociator = cms.untracked.InputTag("VertexAssociatorByPositionAndTracks"),
    lumiInfoTag = cms.untracked.InputTag("LumiInfo", "brilcalc"),
    useVertexFilter = cms.untracked.bool(False),
    compareCollections = parameters["compareCollections"],
    vertexRecoCollections = cms.VInputTag(vcollections)
)




process.analyze =  cms.Path( process.theTruth * process.oldVertexAnalysis )
#process.content = cms.EDAnalyzer("EventContentAnalyzer")
#process.dump = cms.Path(process.content)

   


print( "=============================================")
print( process.source)
print( "=============================================")

# Schedule definition
if DO_VTX_RECO:
    print( "running vertex reco and analysis, last message from pv_cfg")
    process.schedule = cms.Schedule(process.MPV, process.analyze)#,process.dump) 
else:
    print( "running analysis only, last message from pv_cfg")
    process.schedule = cms.Schedule(process.analyze)


print("pv_cfg.py done")
