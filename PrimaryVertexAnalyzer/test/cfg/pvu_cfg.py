import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList



from Configuration.StandardSequences.Eras import eras

# obosolete from RecoVertex.PrimaryVertexProducer.TkClusParameters_cff import DA_vectParameters

import sys, os
import subprocess, shlex

# run 4 full event, including timing
# cmsRun cfg/ttbar_cfg.py debug=1:7659:-1.4364:0.1 uniquetrkweight=-0.8

def root_path(f):
    # f is a file path as e.g. das would spit it out:  /store/...../XXX.root
    # prefers the local version if there is one
    if os.path.exists( "/pnfs/psi.ch/cms/trivcat/store/user/werdmann" + f):
        return "root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/werdmann" + f
    elif f.startswith("/eos/cms/store"):
        return "root://cms-xrd-global.cern.ch//"+ f[len("/eos/cms/"):]
    else:
        return "root://xrootd-cms.infn.it/" + f

DO_VTX_RECO = False
use_genHTFilter = False  #  for now, make this false later

args =""
nevent=-1
info = "test"
debug = False
eventsToProcess = "" #cms.untracked.VEventRange('1:4474')
zdumpcenter = 0
zdumpwidth = 100.
use_tp = True
use_hlt = False
use_2file_solution = False
use_lumiInfo = False

json = "goodList.json"
era = "Phase2"  # or Phase2 or Run2_2018,  see eras.__dict__.keys()
geometry = ""
globaltag = ""
source_files = ()
secondary_files =  None #cms.untracked.vstring("")

# twofile solutions
f1 = {}
f2 = {}


# the labels will be overridden when pvreco is done:
#trkTimesLabel = "tofPID:t0safe"          #  
#trkTimeResosLabel = "tofPID:sigmat0safe"
trkTimesLabel = "tofPID:t0"          #  
trkTimeResosLabel = "tofPID:sigmat0"
#trkTimesLabel = "trackExtenderWithMTD:tmtd" 
#trkTimeResosLabel = "trackExtenderWithMTD:sigmatmtd"
timealgo = ""  # i.e. none

parameters={  # can be overwritten by arguments of the same name
  "4D": cms.untracked.bool(False),
  "refit": cms.untracked.bool(False),
  "zMatchWindow_sigmas": cms.untracked.double(3.0),  
  "compareCollections" : cms.untracked.int32(0),
  "selNdof": cms.untracked.double(4.0),
  "selNdofWithBS": cms.untracked.double(2.0),
  "usefit":cms.untracked.bool(False),
#  "use_tp":cms.untracked.bool(False),
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
  "nDump":cms.untracked.int32(10),
  "nDumpTracks":cms.untracked.int32(1),
#  "mintrkweight":cms.untracked.double(0.5),
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
# temp
  "trackTimeQualityThreshold" : cms.untracked.double(0.8),
  "purge" : cms.untracked.int32(0),
  "rho0mode" : cms.untracked.int32(0),  # /nt, as before
  "mergenotc" : cms.untracked.bool(False),
  "mergeafterpurge" : cms.untracked.bool(False),
  "fillzmerge" : cms.untracked.double(0.0060), # default = vertexSize 
  "use_hitpattern" : cms.untracked.int32(0),
  "use_pt" : cms.untracked.int32(0),
  "runInBlocks" : cms.bool(False),
  "block_size" : cms.uint32(512),
  "overlap_frac" : cms.double(0.5)
}

print(">>>>>>>>>>>>>>>>>>>>>>   this is pvu_cfg.py   <<<<<<<<<<<<<<<<<<<<<<<<<")

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
            files = [inputfile]
            source_files =  cms.untracked.vstring( root_path(inputfile) )
        else:
            files = [l.strip() for l in open(inputfile).readlines()]
            source_files = cms.untracked.vstring(*tuple([root_path(f) for f in files]))

    elif key == "json":
        json = None if value in ("None","") else value

    elif key == "era":
        era = value

    elif key == "geometry":
        geometry = value

    elif key == "globaltag":
        globaltag = value

    elif key == "4D":
        parameters["4D"] = cms.untracked.bool(value == "True")
        print("4D = ", parameters["4D"])
        if value == "True":
            timealgo = "fromTracksPID"

    elif key == "redopv" or key == "redoPV":
        DO_VTX_RECO = (value in ("True","yes","Yes") )

    elif key == "use_genHTFilter":
        use_genHTFilter = (value in ("True","yes","Yes") )

    elif key in ("usetp", "use_tp", "useTP", "use_TP"):
        use_tp = (value in ("True", "yes", "Yes"))

    elif key in ("use2file",):
        use_2file_solution = (value in ("True", "yes", "Yes"))

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
        print("blablabla")
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
        print( "!! pvu_cfg.py  :  unknown key ",key)


if use_2file_solution and len(files) > 0:
    fpath = "/store" + files[0].split("/store")[-1]
    raw_files = subprocess.check_output(shlex.split('dasgoclient -query="parent file=%s"' % fpath)).strip().split(b'\n')
    print( raw_files)
    secondary_files = cms.untracked.vstring([root_path(f.decode('utf-8')) for f in raw_files])

print( "pvu_cfg.py")
print( "args           ",args)
print( "events         ",nevent)
print( "info           ",info)
print( "json           ",json)
print( "era            ",era)
print( "geometry       ",geometry)
print( "globaltag      ",globaltag)
print( "DO_VTX_RECO    ",DO_VTX_RECO)
print( "all parameters")
print( parameters)
if len(source_files) == 0:
    print( "empty source file list!")
else:
    print( "source files = ",source_files)
    print( "secondary files = ", secondary_files)

if era == "Phase2":
    process = cms.Process("RERECO", eras.Phase2)
    autotag = "auto:phase2_realistic"
    parameters["maxEta"] = cms.double(4.0)
elif era == "Run3":
    process = cms.Process("RERECO", eras.Run3)
    #autotag = "111X_mcRun3_2021_realistic_v4"
    #autotag = "113X_mcRun3_2021_realistic_v1"
    #autotag = "113X_mcRun3_2021_realistic_v4"
    #autotag = "auto:phase1_2021_realistic"
    autotag = "auto:phase1_2023_realistic"
    parameters["4D"] = cms.untracked.bool(False)
elif era == "Run2_2018":
    process = cms.Process("RERECO", eras.Run2_2018)
    autotag = "auto:run2_data"
    parameters["4D"] = cms.untracked.bool(False)
else:
    print( "unknown era", era)
    print( (era == "Run2_2018"), ">%s<"%era)
    sys.exit(1)
# for auto: tags see $CMSSW_RELEASE_BASE/python/Configuration/AlCa/autoCond.py
# to find the tag used for producing a file 
# edmProvDump root:....
#ESSource: GlobalTag RECO
#  @module_label: string tracked  = 'GlobalTag'
#  globaltag: string tracked  = '111X_mcRun3_2021_realistic_v4'
#
##module label e.g., "SimG4Objects" or "TrackProducer". 
##product instance label (=instance name ?) 


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000   # don't be too noisy
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreamsMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#

#detectorVersionDict = {
#    ("O8","T24","C17","M10","F6","I14") : "D86",
#    ("O9","T24","C17","M10","F6","I15") : "D88",
#    ("O9","T30","C17","M10","F6","I15") : "D91",
#    ("O9","T24","C18","M10","F6","I15") : "D92",
#    ("O9","T24","C19","M10","F6","I15") : "D93",
#    ("O9","T24","C20","M10","F8","I15") : "D94",
#    ("O9","T31","C17","M10","F6","I16") : "D95",
#    ("O9","T31","C18","M10","F6","I16") : "D96",
#    ("O9","T25","C17","M10","F6","I15") : "D97",
#    ("O9","T32","C17","M10","F6","I16") : "D98",
#    ("O9","T32","C18","M10","F6","I16") : "D99",
#}

if era.startswith("Run2"):
    process.load('Configuration.StandardSequences.GeometryDB_cff') #?
elif era == "Run3":
    process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
elif era == "Phase2" and geometry == "":
    process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')  
elif geometry in ("D76", "T21"):
    process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')  
elif geometry in ("D77", "T22"):
    process.load('Configuration.Geometry.GeometryExtended2026D77Reco_cff')  
elif geometry in ("D80", "T25"):
    process.load('Configuration.Geometry.GeometryExtended2026D80Reco_cff')  
elif geometry in ("D81", "T26"):
    process.load('Configuration.Geometry.GeometryExtended2026D81Reco_cff')  
elif geometry in ("D88",):
    process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')  
elif geometry in ("D89", "T28"):
    process.load('Configuration.Geometry.GeometryExtended2026D89Reco_cff')  
elif geometry in ("D90", "T29"):
    process.load('Configuration.Geometry.GeometryExtended2026D90Reco_cff')  
elif geometry in ("D95", ):
    process.load('Configuration.Geometry.GeometryExtended2026D95Reco_cff')
elif geometry in ("D96", "T31"):
    process.load('Configuration.Geometry.GeometryExtended2026D96Reco_cff')
elif geometry in ("D97", ):
    process.load('Configuration.Geometry.GeometryExtended2026D97Reco_cff')
elif geometry in ("D98", ):
    process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')
elif geometry in ("D99", "T32" ):
    process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')
elif geometry in ("D110",  ):
    process.load('Configuration.Geometry.GeometryExtended2026D110Reco_cff')
    
# $CMSSW_RELEASE_BASE/src/Configuration/Geometry/python/GeometryExtended2026D95Reco_cff.py
print("*"*80,"\n use_tp",use_tp)

if use_tp:
    process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
    process.theTruth = cms.Sequence(
        process.tpClusterProducer *
        process.quickTrackAssociatorByHits *
        process.trackingParticleRecoTrackAsssociation
    )
if False:
    # see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTrackMCTruth
    process.load("Validation.RecoTrack.TrackValidation_cff")
    process.theTruth = cms.Sequence(
        process.tpClusterProducer +
        process.quickTrackAssociatorByHits +
        process.trackingParticleRecoTrackAsssociation
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(nevent))


# Input source
if secondary_files is None:
    process.source = cms.Source("PoolSource",
                                fileNames = source_files,
                                #eventsToProcess = cms.untracked.VEventRange()
    )
else:
    process.source = cms.Source("PoolSource",
                                fileNames = source_files,
                                secondaryFileNames = secondary_files
                                #inputCommands = cms.untracked.vstring('drop *_muonReducedTrackExtras_*_*')
                                #eventsToProcess = cms.untracked.VEventRange()
                            )


if not eventsToProcess == "":
    process.source.eventsToProcess = cms.untracked.VEventRange( eventsToProcess )


if json is not None  and os.path.exists(json) and (len(source_files)==0 or source_files[0].find("SIM")<0):
    process.source.lumisToProcess = LumiList.LumiList(filename = json).getVLuminosityBlockRange()


# no HLT for MC
if len(source_files)>0 and source_files[0].find("SIM") >= 0:
    use_hlt = False



#process.Tracer = cms.Service("Tracer")

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)



from Configuration.AlCa.GlobalTag import GlobalTag
if globaltag == "":
    print("using autotag ", autotag)
    process.GlobalTag = GlobalTag(process.GlobalTag, autotag, '')
else:
    print("using globaltag ", globaltag)
    process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')


if DO_VTX_RECO:
    #trkTimesLabel = "tofPID4DnoPID:t0safe" # as opposed to "tofPID:t0safe" , only exist after pvreco
    #trkTimeResosLabel = "tofPID4DnoPID:sigmat0safe"
    process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True)  )

#process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
#process.load("RecoLocalCalo.HcalRecProducers.HBHEIsolatedNoiseReflagger_cfi")
process.load("RecoVertex.Configuration.RecoVertex_cff")

#if DO_VTX_RECO:
#    process.caloTowerForTrk.hbheInput = 'reducedHcalRecHits:hbhereco'
#    process.caloTowerForTrk.hoInput = 'reducedHcalRecHits:horeco'
#    process.caloTowerForTrk.hfInput = 'reducedHcalRecHits:hfreco'
#    process.caloTowerForTrk.ecalInputs = ['reducedEcalRecHitsEB', 'reducedEcalRecHitsEE']

print( process.vertexreco)
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
            print( f"pvu_cfg: changing tkFilter parameter '{par_name}' from '{default}' to '{parameters[par_name]}'")
        setattr(tkFilterParameters, par_name, parameters[par_name])
    except ValueError:
        print( f"pvu_cfg: attribute 'tkFilterParameters.{par_name}' not found ")



# modify the clustering parameters if needed
tkClusParameters = process.unsortedOfflinePrimaryVertices.TkClusParameters.clone()    
print( "original z clustering parameters")
print("---------------------")
print( tkClusParameters)
print("---------------------")
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
        print( "pvu_cfg: attribute TkDAClusParameters.", par_name , " not found ")

    setattr(tkDAClusParameters, "zdumpcenter", parameters[ "zdumpcenter"])  
    setattr(tkDAClusParameters, "zdumpwidth", parameters[ "zdumpwidth"])


print( "modified z clustering parameters")
print("---------------------")
print( tkClusParameters)
print("---------------------")

print("XXX1")
coll = process.unsortedOfflinePrimaryVertices.vertexCollections[0]
print(coll)
print("XXX2")
vertex_time_parameters = process.unsortedOfflinePrimaryVertices.vertexCollections[0].vertexTimeParameters.clone()
print("XXX3")
print(vertex_time_parameters)
print("XXX4")

adaptive = cms.PSet(label=cms.string("adaptive"),
                   algorithm=cms.string("AdaptiveVertexFitter"),
                   chi2cutoff = cms.double(2.5), 
                   minNdof=cms.double(0.0),
                   useBeamConstraint = cms.bool(False),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( 
                       algorithm = cms.string(timealgo),
                       fromTracksPID = cms.PSet(
                           Tstart = cms.double(256.)
                       )
                   )
                   )
test = cms.PSet(label=cms.string("test"),
                   algorithm=cms.string("AdaptiveVertexFitter"),
                   chi2cutoff = cms.double(2.5), 
                   minNdof=cms.double(0.0),
                   useBeamConstraint = cms.bool(False),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( 
                       algorithm = cms.string(timealgo),
                       fromTracksPID = cms.PSet(
                           Tstart = cms.double(-256.)
                       )
                   )
                   )
adaptiveWithBS = cms.PSet(label=cms.string("adaptiveWithBS"),
                   algorithm=cms.string("AdaptiveVertexFitter"),
                   chi2cutoff = cms.double(2.5),
                   minNdof=cms.double(0.0),
                   useBeamConstraint = cms.bool(True),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( algorithm = cms.string(timealgo))
                   )
chisquared = cms.PSet(label=cms.string("chisquared"),
                   algorithm=cms.string("AdaptiveChisquareVertexFitter"),
                   chi2cutoff = cms.double(2.5),
                   mintrkweight = cms.double(0.0),
                   minNdof=cms.double(0.0),
                   useBeamConstraint = cms.bool(False),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( algorithm = cms.string(timealgo))
                   )  
chisquaredWithBS = cms.PSet(label=cms.string("chisquaredWithBS"),
                   algorithm=cms.string("AdaptiveChisquareVertexFitter"),
                   chi2cutoff = cms.double(2.5),
                   minNdof=cms.double(0.0),
                   useBeamConstraint = cms.bool(True),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( algorithm = cms.string(timealgo))
                   )  
cluster = cms.PSet(label=cms.string("cluster"),
                   algorithm=cms.string(""),
                   chi2cutoff = cms.double(0),
                   minNdof=cms.double(0.0),
                   useBeamConstraint = cms.bool(False),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( algorithm = cms.string(timealgo))
                   )  
MPVBS = cms.PSet(label=cms.string("MPVBS"),
                   algorithm=cms.string("MultiPrimaryVertexFitter"),
                   chi2cutoff = cms.double(2.5),
                   mintrkweight = cms.double(0.2),
                   minNdof=cms.double(2.0),
                   useBeamConstraint = cms.bool(True),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( algorithm = cms.string(timealgo))
                   )
AMPVBS = cms.PSet(label=cms.string("AMPVBS"),
                   algorithm=cms.string("AMultiPrimaryVertexFitter"),
                   chi2cutoff = cms.double(5.0),
                   mintrkweight = cms.double(0.2),
                   minNdof=cms.double(2.0),
                   useBeamConstraint = cms.bool(True),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( algorithm = cms.string(timealgo))
                   )
weightedmean = cms.PSet(label=cms.string("weightedmean"),
                   algorithm=cms.string("WeightedMeanFitter"),
                   minNdof= cms.double(0.0),
                   chi2cutoff = cms.double(2.5),
                   useBeamConstraint = cms.bool(False),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( algorithm = cms.string(timealgo))
                   )
weightedmeanWithBS = cms.PSet(label=cms.string("weightedmeanWithBS"),
                   algorithm=cms.string("WeightedMeanFitter"),
                   minNdof= cms.double(0.0),
                   chi2cutoff = cms.double(2.5),
                   useBeamConstraint = cms.bool(True),
                   maxDistanceToBeam = cms.double(1.0),
                   vertexTimeParameters = cms.PSet( algorithm = cms.string(timealgo))
                   )



process.testVertices = cms.EDProducer(
    "PrimaryVertexProducer",

    verbose = parameters["verboseProducer"],#cms.untracked.bool(True),
    TrackLabel = cms.InputTag("generalTracks"),
    beamSpotLabel = cms.InputTag("offlineBeamSpot"),
    
    TkFilterParameters = tkFilterParameters,

    TkClusParameters = tkClusParameters,

    #vertexCollections = cms.VPSet([adaptive, chisquared, chisquaredWithBS, adaptiveWithBS, MPVBS, AMPVBS,weightedmean, weightedmeanWithBS]),
    vertexCollections = cms.VPSet([adaptive, test, chisquared]),

    isRecoveryIteration = cms.bool(False),
    recoveryVtxCollection = cms.InputTag("")
)
process.MPV = cms.Path(process.testVertices)


# the  list "vertexCollections" defined above is for the producer
# the list "vcollections" defined below is for the analyzer and may be different
#
#
#vcollections = ["offlinePrimaryVertices", "offlinePrimaryVerticesWithBS"]
#vcollections = ["offlinePrimaryVertices","offlinePrimaryVerticesWithBS","offlineMultiPrimaryVertices","offlineMultiPrimaryVertices:WithBS"]
#vcollections = ["offlinePrimaryVertices","offlineMultiPrimaryVertices"]
#vcollections = ["testVertices:adaptive","testVertices:chisquared","offlinePrimaryVertices","testVertices:cluster","testVertices:MPVBS","testVertices:weightedmean","offlinePrimaryVerticesWithBS"]
#vcollections = ["testVertices:adaptive","testVertices:weightedmean"]
#vcollections = ["testVertices:MPVBS","testVertices:chisquared","testVertices:chisquaredWithBS","offlinePrimaryVerticesWithBS","offlinePrimaryVertices4D","offlinePrimaryVertices"]
#vcollections = ["testVertices:MPVBS","testVertices:AMPVBS","testVertices:chisquared","testVertices:chisquaredWithBS","offlinePrimaryVerticesWithBS","testVertices:adaptiveWithBS","testVertices:adaptive","offlinePrimaryVertices","testVertices:weightedmean","testVertices:weightedmeanWithBS"]
vcollections = ["offlinePrimaryVertices","testVertices:adaptive","testVertices:test"]
   
print(">"*40)
print("zMatchWindow_sigmas=",parameters["zMatchWindow_sigmas"])
print("f4D=",parameters["4D"])
print("<"*40)

# analysis
process.oldVertexAnalysis = cms.EDAnalyzer("PrimaryVertexAnalyzer4PU",
    info=  cms.untracked.string(info),
    f4D = parameters["4D"],
    frefit = parameters["refit"],
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
    RECO = parameters["reco"],
    MINIAOD = parameters["miniaod"],
    use_tp = cms.untracked.bool(use_tp),
    fill_track_histos = parameters["fill_track_histos"],
    track_timing = cms.untracked.bool(True),
    TkFilterParameters = tkFilterParameters,
    trackingParticleCollection = cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackingVertexCollection = cms.untracked.InputTag("mix", "MergedTrackTruth"),
    trackAssociatorMap = cms.untracked.InputTag("trackingParticleRecoTrackAsssociation"),
    TrackTimesLabel = cms.untracked.InputTag(trkTimesLabel),  
    TrackTimeResosLabel = cms.untracked.InputTag(trkTimeResosLabel),
    TrackTimeQualityMapLabel = cms.untracked.InputTag("mtdTrackQualityMVA:mtdQualMVA"),
    TrackTimeQualityThreshold = cms.untracked.double( parameters["trackTimeQualityThreshold"].value()),
    vertexAssociator = cms.untracked.InputTag("VertexAssociatorByPositionAndTracks"),
    lumiInfoTag = cms.untracked.InputTag("LumiInfo", "brilcalc"),
    useVertexFilter = cms.untracked.bool(False),
    compareCollections = parameters["compareCollections"],
    vertexRecoCollections = cms.VInputTag(vcollections)
)

# lumi, 
# from /RecoLuminosity/LumiProducer/test/TestLumiProducerFromBrilcalc_cfg.py
# https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#LumiCMSSW
# so how do I get csv file, back to running brilcalc
if use_lumiInfo:
    process.LumiInfo = cms.EDProducer('LumiProducerFromBrilcalc',
                                      lumiFile = cms.string("./lumiFile.csv"),
                                      throwIfNotFound = cms.bool(False),
                                      doBunchByBunch = cms.bool(False))


##process.p = cms.Path(process.LumiInfo*process.test)

#--------------------------------------------------------------------------
# Trigger Selection

process.hltSelection = cms.EDFilter("HLTHighLevel",
                               TriggerResultsTag =
                               cms.InputTag("TriggerResults","","HLT"),
                               HLTPaths = cms.vstring(["HLT_ZeroBias_v*"]),
                               eventSetupPathsKey = cms.string(''),
                               andOr = cms.bool(True),
                               throw = cms.bool(False)  )

process.genHTFilter = cms.EDFilter("GenHTFilter",
   src = cms.InputTag("ak4GenJets"), #GenJet collection as input
   jetPtCut = cms.double(30.0), #GenJet pT cut for HT
   jetEtaCut = cms.double(4.5), #GenJet eta cut for HT
   genHTcut = cms.double(1000.0) #genHT cut
)


print( "============================================")

if use_tp and not use_genHTFilter:
    print ("using tracking particles, no genHTFilter")
    process.analyze =  cms.Path( process.theTruth * process.oldVertexAnalysis )
elif not use_genHTFilter:
    print ("using no tracking particles, no genHTFilter")
    process.analyze =  cms.Path( process.oldVertexAnalysis )
elif use_tp and use_genHTFilter:
    print( "using genHTFilter")
    process.analyze =  cms.Path( process.genHTFilter * process.theTruth * process.oldVertexAnalysis )
elif use_hlt:
    print( "using hlt seletion")
    process.analyze =  cms.Path( process.LumiInfo * process.hltSelection * process.oldVertexAnalysis )
else:
    print( "no fiters, no tp")
    process.analyze =  cms.Path( process.LumiInfo * process.oldVertexAnalysis )


print( "=============================================")
print( process.source)
print( "=============================================")

if False:

    # only dump event content, don't run the analyzer
    process.content = cms.EDAnalyzer("EventContentAnalyzer")
    process.dump = cms.Path(process.content)
    process.schedule = cms.Schedule(process.dump)

else:

    # Schedule definition
    if DO_VTX_RECO:
        print( "running vertex reco and analysis, last message from pvu_cfg")
        #process.schedule = cms.Schedule(process.MPV, process.dump,process.analyze) 
        process.schedule = cms.Schedule(process.MPV, process.analyze) 
        #process.schedule = cms.Schedule(process.MPV) #vertex reco only 
    else:
        print( "running analysis only, last message from pvu_cfg")
        #process.schedule = cms.Schedule(process.dump, process.analyze)
        process.schedule = cms.Schedule(process.analyze)
print("pvu_fg.py done")
