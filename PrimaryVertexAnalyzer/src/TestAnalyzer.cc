#include "usercode/PrimaryVertexAnalyzer/interface/TestAnalyzer.h"
#include "usercode/PrimaryVertexAnalyzer/interface/FFA.h"
//  std::cout << "XDBG " << __func__ << " : " << __LINE__ << std::endl;

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "MagneticField/Engine/interface/MagneticField.h"

// reco track and vertex
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// simulated vertices,..., add <use name=SimDataFormats/Vertex> and <../Track>
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include <SimDataFormats/Track/interface/SimTrack.h>
#include <SimDataFormats/Track/interface/SimTrackContainer.h>

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// Lumi
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiInfo.h"
// obsolete?
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "DataFormats/Common/interface/ConditionsInEdm.h"

// AOD et al
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//generator level + CLHEP
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"

// TrackingParticle
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

// fit
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/TrajectoryParametrization/interface/TrajectoryStateExceptions.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

// Root
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TProfile.h>
#include <TRandom.h>
//#include <TEfficiency.h>

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

// cluster stufff
//#include "DataFormats/TrackRecoTrack.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"  // starting with CMSSW_11
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

// some of this is probably redundant or unused
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
//#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

// TCDS / reset orbit info
#include "DataFormats/TCDS/interface/TCDSRecord.h"

#include <sstream>
#include <fstream>

#include <assert.h>
#include <algorithm>


#include<chrono>

using namespace edm;
using namespace reco;
using namespace std;
//
// constants, enums and typedefs
//
//
// static data member definitions
//
namespace {
  template <typename T, size_t N>
  std::array<T, N + 1> makeLogBins(const double min, const double max) {
    const double minLog10 = std::log10(min);
    const double maxLog10 = std::log10(max);
    const double width = (maxLog10 - minLog10) / N;
    std::array<T, N + 1> ret;
    ret[0] = std::pow(10, minLog10);
    const double mult = std::pow(10, width);
    for (size_t i = 1; i <= N; ++i) {
      ret[i] = ret[i - 1] * mult;
    }
    return ret;
  }
}  // namespace

//
// constructors and destructor
//
TestAnalyzer::TestAnalyzer(const ParameterSet& iConfig)
    : verbose_(iConfig.getUntrackedParameter<bool>("verbose", false)),
      //do_pixel_cluster_analysis_( iConfig.getUntrackedParameter<bool>( "doPixelClusterAnalysis", false ) ),
      //do_pixel_cluster_analysis_all_layers_( iConfig.getUntrackedParameter<bool>( "doPixelClusterAnalysisAllLayers", false ) ),
      do_vertex_analysis_(iConfig.getUntrackedParameter<bool>("doVertexAnalysis", true)),
      doMatching_(iConfig.getUntrackedParameter<bool>("matching", false)),
      analyzeLS_(iConfig.getUntrackedParameter<int>("LS", -1)),
      DEBUG_(false),
      eventcounter_(0),
      dumpcounter_(0),
      ndump_(0),
      ndump_tracks_(0),
      run_(0),
      luminosityBlock_(0),
      event_(0),
      bunchCrossing_(0),
      orbitNumber_(0),
      fBfield_(0.),
      simUnit_(1.0),    // cm everywhere starting with CMSSW_1_2_x ??
      simtUnit_(1.e9),  // nanoseconds  in rec, do the same for sim
      zmatch_(iConfig.getUntrackedParameter<double>("zmatch", 0.0500)),
      sigmaZoverride_(iConfig.getUntrackedParameter<double>("sigmaZ", 0.0)),
      sigmaZ_(0),
      wxy2_(0.),
      outputFile_(iConfig.getUntrackedParameter<std::string>("outputFile")),
      theTrackFilter(iConfig.getParameter<edm::ParameterSet>("TkFilterParameters")),
      recoTrackProducer_(iConfig.getUntrackedParameter<std::string>("recoTrackProducer")),
      vecPileupSummaryInfoToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag(std::string("slimmedAddPileupInfo")))),
      genParticleCollection_Token_(consumes<GenParticleCollection>(edm::InputTag(std::string("genParticles")))),
      recoTrackCollectionToken_(consumes<reco::TrackCollection>(
          edm::InputTag(iConfig.getUntrackedParameter<std::string>("recoTrackProducer")))),
      recoBeamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      edmView_recoTrack_Token_(consumes<edm::View<reco::Track>>(
          edm::InputTag(iConfig.getUntrackedParameter<std::string>("recoTrackProducer")))),
      edmSimVertexContainerToken_(consumes<edm::SimVertexContainer>(iConfig.getParameter<edm::InputTag>("simG4"))),
      edmSimTrackContainerToken_(consumes<edm::SimTrackContainer>(iConfig.getParameter<edm::InputTag>("simG4"))),
      edmHepMCProductToken_(consumes<edm::HepMCProduct>(
          edm::InputTag(std::string("generatorSmeared")))),
      trackingParticleCollectionToken_(consumes<TrackingParticleCollection>(
          iConfig.getUntrackedParameter<edm::InputTag>("trackingParticleCollection"))),
      trackingVertexCollectionToken_(
          consumes<TrackingVertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackingVertexCollection"))),
      simToRecoAssociationToken_(
          consumes<reco::SimToRecoCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociatorMap"))),
      recoToSimAssociationToken_(
          consumes<reco::RecoToSimCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociatorMap"))),
      vertexAssociatorToken_(consumes<reco::VertexToTrackingVertexAssociator>(
          iConfig.getUntrackedParameter<edm::InputTag>("vertexAssociator"))),
      lumiDetailsToken_(consumes<LumiDetails, edm::InLumi>(edm::InputTag("lumiProducer"))),
      lumiSummaryToken_(consumes<LumiSummary, edm::InLumi>(edm::InputTag("lumiProducer"))),
      lumiInfoToken_(consumes<LumiInfo>(iConfig.getUntrackedParameter<edm::InputTag>("lumiInfoTag"))),
      pdtToken_(esConsumes()),
      transientTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      trackerTopologyToken_(esConsumes())
{

  // v34 and higher, avoid biases
  minNumberOfRecTrks_ = 0.;  // FIXME make this configurable or maybe even better obsolete (are these empty BXs?)
  minNumberOfSelTrks_ = 0.;  // FIXME make this configurable

  // definition of visible vertices
  etaMaxVisible_ = 4.0;
  ptMinVisible_ = 0.2;
  numTrkHitsVisible_ = 4;

  fill_track_histos_ = iConfig.getUntrackedParameter<bool>("fill_track_histos", false);
  selNdofNoBS_ = iConfig.getUntrackedParameter<double>("selNdof", 4.);
  std::cout << "TestAnalyzer: selNDof_(noBS) = " << selNdofNoBS_ << std::endl;
  selNdofWithBS_ = iConfig.getUntrackedParameter<double>("selNdofWithBS", 7.);
  std::cout << "TestAnalyzer: selNDofWithBS_ = " << selNdofWithBS_ << std::endl;
  selNdof_ = selNdofNoBS_; // to be changed later according to the collection name
  ndof0trk_ = 0.;

  min_trk_in_vtx_weight_ = 0.2;



  
  reco_vertex_collections_ = iConfig.getParameter<std::vector<edm::InputTag>>("vertexRecoCollections");

  for (auto const& l : reco_vertex_collections_) {
    auto l_encode = l.encode();
    std::replace(l_encode.begin(), l_encode.end(), ':', '_');
    std::cout << "TestAnalyzer: vertex collection [" << l_encode << "]" << std::endl;
    vertexCollectionLabels_.push_back(l_encode);

    auto token = edm::EDGetTokenT<reco::VertexCollection>(consumes<reco::VertexCollection>(edm::InputTag(l)));
    vertexCollectionTokens_[l_encode] = token;
  }

  // open output file to store histogram
  rootFile_ = TFile::Open(outputFile_.c_str(), "RECREATE");
  info_ = new TObjString("info");
  info_->SetString(iConfig.getUntrackedParameter<std::string>("info", "").c_str());
  build_ = new TObjString("build");
  cout << "build = " << Form("%s %s", __DATE__, __TIME__) << endl;
  build_->SetString(Form("%s %s", __DATE__, __TIME__));

  veryverbose_ = iConfig.getUntrackedParameter<bool>("veryverbose", false);
  if (verbose_) {
    cout << "TestAnalyzer: veryverbose = " << veryverbose_ << endl;
    cout << "    extra level of verbosity " << endl;
    cout << endl;
  }

  if (!do_vertex_analysis_)
    return;

  if (verbose_) {
    cout << "TestAnalyzer: zmatch=" << zmatch_ << endl;
    cout << "sigmaZ = " << sigmaZoverride_ << endl;
    cout << "     if 0.0 : use the value from the beamspot" << endl;
    cout << endl;
  }

  nPUmin_ = std::abs(iConfig.getUntrackedParameter<int>("PUmin", 0));
  nPUmax_ = std::abs(iConfig.getUntrackedParameter<int>("PUmax", 1000000));
  if (verbose_) {
    cout << "nPUMin = " << nPUmin_ << endl;
    cout << "nPUMax = " << nPUmax_ << endl;
    cout << "     in MC, only analyze events with nPUmin <  N(simvertex) < nPUmax " << endl;
    cout << endl;
  }

  useVertexFilter_ = iConfig.getUntrackedParameter<bool>("useVertexFilter", false);

  nEventSummary_ = iConfig.getUntrackedParameter<int>("eventSummaries", 0);
  ndump_ = iConfig.getUntrackedParameter<int>("nDump", 0);
  ndump_tracks_ = iConfig.getUntrackedParameter<int>("nDumpTracks", 0);

  nCompareCollections_ =
      iConfig.getUntrackedParameter<int>("compareCollections", 0);  //-1= compare all, >0 dump n events

  zmatch_ = iConfig.getUntrackedParameter<double>("zmatch", 0.0500);
  if (verbose_) {
    cout << "TestAnalyzer: zmatch =" << zmatch_ << endl;
    cout << "     cut-off for matching sim to reco by z" << endl;
    cout << endl;
  }

  zWosMatchMax_ = iConfig.getUntrackedParameter<double>("zMatchMax", 1.);
  if (verbose_) {
    cout << "TestAnalyzer: zMatchMax = " << zWosMatchMax_ << endl;
    cout << "     cut-off for insane recvertex <-> simvertex  matches" << endl;
    cout << "     (TrackingParticles, matching by weight/sigma^2) " << endl;
    cout << "     default is 1 cm, configurable for exotic events where" << endl;
    cout << "     all tracks appear far away from the vertex   " << endl;
    // such as LongLivedChi0ToNuLL_MSquark-1000_MChi-148_TuneZ2Star_8TeV-pythia6
    cout << endl;
  }
  
  zwindow_sigmas_ = iConfig.getUntrackedParameter<double>("zMatchWindow_sigmas", 3.0);
  if (verbose_) {
    cout << "TestAnalyzer: zwindow_sigmas_ = " << zwindow_sigmas_ << endl;
    cout << "     for z-coordinate based recvertex <-> simvertex  matchin" << endl;
    cout << "     window size in multiples of  the reconstructed sigma(z)" << endl;
    cout << endl;
  }

  RECO_ = iConfig.getUntrackedParameter<bool>("RECO", false);
  if (verbose_) {
    cout << "TestAnalyzer: RECO = " << RECO_ << endl;
    cout << "      use RECO information (pixel hits and some trackextra)" << endl;
    cout << endl;
  }

  autoDumpCounter_ = iConfig.getUntrackedParameter<int>("autodump", 0);
  if (verbose_) {
    cout << "TestAnalyzer: autodump = " << autoDumpCounter_ << endl;
    cout << "      dump detailed information about the first <autodump> events " << endl;
    cout << endl;
  }

  // extras
  extraInfoToken_ = consumes<std::vector<float>>(edm::InputTag("testVertices","extraInfo"));
  clusteringCPUtimeToken_ = consumes<float>(edm::InputTag("testVertices","clusteringCPUTime"));
  

  trkhiptmin_ = 3.0;
  trkloptmax_ = 1.0;
  trkcentraletamax_ = 1.5;
  trackAssociatorMin_ = 0.5;

  /* initialize counters */
  eventcounter_ = 0;
  emptyeventcounter_ = 0;
  dumpcounter_ = 0;
  eventSummaryCounter_ = 0;
  nEventNsel_ = 0;

  currentLS_ = -1;

  // profiling
  timers_.clear();
  
}

TestAnalyzer::~TestAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete rootFile_;
}

//
// member functions
//

std::map<std::string, TH1*> TestAnalyzer::bookVertexHistograms(TDirectory * dir) {
    std::map<std::string, TH1*> h;  // will be returned
    
    /* Efficiency */
    dir->mkdir("efficiency")->cd();  // Create the "efficiency" directory and navigate into it

    // Create the histogram for reconstruction efficiency
    TH1I* hist_reco_matches = new TH1I("SEEfficiency", " SE Reconstruction Efficiency", 2, 0, 2);
    
    // Add the histogram to the map using addn
    addn(h, hist_reco_matches);
    std::cout << "Creating histogram: SEEfficiency" << std::endl;

    // Set bin labels for better visualization
    hist_reco_matches->GetXaxis()->SetBinLabel(1, "Matched vertex");
    hist_reco_matches->GetXaxis()->SetBinLabel(2, "Unmatched vertex");

    // Set histogram style: fill color and bar width
    hist_reco_matches->SetFillColor(kOrange + 2);  // Set fill color for the bars
    hist_reco_matches->SetBarWidth(0.8);           // Adjust the bar width

    //set y axis title
    hist_reco_matches->GetYaxis()->SetTitle("Count ");
    hist_reco_matches->SetMinimum(0);
    hist_reco_matches->SetMaximum(150);

    // new histogramm 
    // Definition of the 2D histogram, which shows the simulated vs recon position (on z-axis) in one plot for SE
    TH2F* SERecoVsSimZPosition = new TH2F("SERecoVsSimZPosition", "SE Reconstructed vs. True Z-Position", 100, -10, 10, 100, -10, 10);
    addn(h, SERecoVsSimZPosition );

    // another new histogramm
    // Definition of the 2D histogram, which shows the simulated vs recon position (on z-axis) in one plot for PU
    TH2F* PURecoVsSimZPosition = new TH2F("PURecoVsSimZPosition", "PU Reconstructed vs. True Z-Position", 100, -10, 10, 100, -10, 10);
    addn(h, PURecoVsSimZPosition);
    PURecoVsSimZPosition ->Draw("COLonl"); // should add a legend, doesent work tough might remove or find alternative. t pallet axis should work but dont have time for this now
    // for now just manually add in TBrowser

    // New histogramm
    // SE tracks purity
    TH1F *SETracksPurity = new TH1F("SETracksPurity", "SE Tracks Purity", 100, 0, 101);
    addn(h, SETracksPurity);

    // New histogramm
    // SE tracks Efficiency
    TH1F *SETracksEfficiency = new TH1F("SETracksEfficiency", "SE Tracks Efficiency", 100, 0, 101);
    addn(h, SETracksEfficiency);

    // New histogramm
    // PU tracks purity
    TH1F *PUTracksPurity = new TH1F("PUTracksPurity", "PU Tracks Purity", 100, 0, 101);
    addn(h, PUTracksPurity);

    // New histogramm
    // PU tracks Efficiency
    TH1F *PUTracksEfficiency = new TH1F("PUTracksEfficiency", "PU Tracks Efficiency", 100, 0, 101);
    addn(h, PUTracksEfficiency);

    //new histogram
    //definition of 1 H hist
    // SE Resolution
    TH1F *SEResolution = new TH1F("SEResolution","SE Reconstructed vs. True Z-Position position difference", 200, -0.01, 0.01);
    addn(h, SEResolution);
    //set bar width to 0.1
     SEResolution->SetBarWidth(0.1);   
      //new histogram
    //definition of 1 H hist
    // PU Resolution
    TH1F *PUResolution = new TH1F("PUResolution","PU Reconstructed vs. True Z-Position position difference", 200, -0.2, 0.2);
    addn(h, PUResolution);
    //set bar width to 0.1
     PUResolution->SetBarWidth(0.1);
    //new histogram
    //definition of an 2H histogram
    // PU Confusion Matrix
      TH1F *PUConfusionMatrixCategorialC1 = new TH1F("reco_vs_true_z_position_hist_categorial_c1", "Freq. one reconstructed per simulated vertex;Vertex Z Position;Frequency", 100, -30, 30);
      addn(h, PUConfusionMatrixCategorialC1);
      // new histogram
      // definition of an 2H histogram
      TH1F *PUConfusionMatrixCategorialC2 = new TH1F("reco_vs_true_z_position_hist_categorial_c2", "Freq. mult reconstructed per simulated vertex;Vertex Z Position;Frequency", 100, -30, 30);
      addn(h, PUConfusionMatrixCategorialC2);
      // new histogram
      // definition of an 2H histogram
      TH1F *PUConfusionMatrixCategorialC3 = new TH1F("reco_vs_true_z_position_hist_categorial_c3", "Freq. fake vertices;Vertex Z Position;Frequency", 100, -30, 30);
      addn(h, PUConfusionMatrixCategorialC3);

     // new histogram
      // this histogramm shows distance betweeen point in 3 d space and plane
      // we can later adapt this histogram to for example only show the distance between the border of the subspace and the fake vertices etc
      TH1F *TrueE3DDistanceToPlane = new TH1F("True_3D_point_to_plane_distance", "Distance between 3d point to plane ", 100, 0, 30);
      addn(h, TrueE3DDistanceToPlane);



      // Return to the base directory to maintain proper organization
      dir->cd();

      // If you plan to add more histograms or profiles later, follow a similar structure:
      // - Create or navigate into the appropriate directory using mkdir and cd.
      // - Use addn to add histograms to the map h.
      // - Ensure proper cleanup by returning to the base directory using dir->cd().
      cout << "The bookVertexHistograms method was run sucesssfully ";
      return h; // Return the map containing all histograms
}


void TestAnalyzer::bookTrackHistograms(const char * directory_name)
{
}

void TestAnalyzer::bookSimPVHistograms(const char * directory_name)
{
}

void TestAnalyzer::bookEventHistograms(const char * directory_name)
{
}


void TestAnalyzer::get_luminosity_infos(const edm::Event& iEvent) {
}

void TestAnalyzer::get_particle_data_table(const edm::EventSetup& iSetup) {
  try {
    //iSetup.getData(pdt_);
    pdt_ = iSetup.getHandle(pdtToken_);
  } catch (const Exception&) {
    std::cout << "Some problem occurred with the particle data table. This may not work !" << std::endl;
  }
}

bool TestAnalyzer::get_beamspot_data(const edm::Event& iEvent) {
  // load beam spot, fill some histograms
  // return false if no usable beamspot data was found

  if (iEvent.getByToken(recoBeamSpotToken_, recoBeamSpotHandle_)) {
    vertexBeamSpot_ = *recoBeamSpotHandle_;
    wxy2_ = pow(vertexBeamSpot_.BeamWidthX(), 2) + pow(vertexBeamSpot_.BeamWidthY(), 2);
    wx_ = vertexBeamSpot_.BeamWidthX();
    wy_ = vertexBeamSpot_.BeamWidthY();
    sigmaZ_ = vertexBeamSpot_.sigmaZ(); // unsicherheit von Beamspot in z-Achse
    sigmaT_ = sigmaZ_ / 2.998e1;  // c in cm/ns
    if (sigmaZ_ < 0.1) {
      reportEvent("crazy small sigma Z");
      return false;
    }

    if (sigmaZ_ > 9) {
      reportEvent("huge sigma Z");
      return false;
    }

    if (std::abs(vertexBeamSpot_.z0()) > 2.) {
      reportEvent("suspicious beamspot position  %6.3f", vertexBeamSpot_.z0()); // ist das weil bs im zentrum des Detektors sein soll?
      return false;
    }

  } else {
    sigmaZ_ = 0;
    reportEvent("no beamspot found, aborting");
    return false;
  }

  // abort without useful beamspot data
  if (sigmaZ_ == 0) { // in case somthing above went wrong
    return false;
  }

  return true;
}



// get tracks
/********************************************************************************************************/
bool TestAnalyzer::get_reco_and_transient_tracks(const edm::EventSetup& iSetup,
                                                             const edm::Event& iEvent,
                                                             Tracks& tracks) {
/********************************************************************************************************/
  // requires beamspot
  if (!iEvent.getByToken(edmView_recoTrack_Token_, tracks.trackCollectionH)) {
    if (verbose_) {
      report_counted( " no reco tracks found, bailing out", 10);
    }
    return false;
  }

  const View<reco::Track>* recTrks = tracks.trackCollectionH.product();

  edm::Handle<edm::ValueMap<float>> trackTimesH;
  edm::Handle<edm::ValueMap<float>> trackTimeResosH;


  std::vector<reco::TransientTrack> t_tks;
  t_tks = (*theB_).build(tracks.trackCollectionH, vertexBeamSpot_);
  report_counted("TestAnalyzer::get_reco_and_transient_tracks timing requested but not found",1);


  /* fill private track container and make it globally available */
  for (View<reco::Track>::size_type i = 0; i < recTrks->size(); ++i) {

    const reco::TransientTrack* tt = &(t_tks[i]);
    unsigned int key1 = t_tks[i].trackBaseRef().key();
    RefToBase<reco::Track> trb(tracks.trackCollectionH, i);
    unsigned int key2 = trb.key();
    // paranoia is in bloom
    if (key1 != key2) {
      cout << "get_reco_and_transient_tracks : key confusion" << endl;
    }

    auto tk = MTrack(i, &(recTrks->at(i)), t_tks[i], key1, false);
    tk._MTD_timeerror = 1.e10;

    tk._selected = theTrackFilter(*tt);

    tracks.push_back(tk);
  }
  if (recTrks->size() < minNumberOfRecTrks_) {
    emptyeventcounter_++;
    if (emptyeventcounter_ < 100) {
      cout << "Event without tracks skipped   run= " << run_ << " ls = " << luminosityBlock_ << " BX=" << bunchCrossing_
           << "  instBXLumi_=" << 1.e3 * instBXLumi_ << " expected PU = " << lumiPU_ << "    selected tracks"
           << recTrks->size() << std::endl;
    }
    if (emptyeventcounter_ == 100) {
      cout << "This is the last message of this kind" << endl;
    }
    return false;
  }

  return true;
}
/********************************************************************************************************/


   
/********************************************************************************************************/
void  TestAnalyzer::fill_track_to_vertex_pointers(Tracks& tracks) {
/********************************************************************************************************/
  for (auto label : vertexCollectionLabels_) {

    if (recVtxs_[label] == NULL) continue;
    
    unsigned int iv = 0;
    for (auto recv : *recVtxs_[label]) {
      if (recv.tracksSize() > 0){
	for (trackit_t t = recv.tracks_begin(); t != recv.tracks_end(); t++) {
          auto & tk = tracks.from_key(t->key());
          const auto weight = recv.trackWeight(*t);
          assert(tk.key() == t->key());
          if( tk.get_recv(label) == NO_RECVTX ||  tk.get_weight(label) < weight){
            tk._recv[label] = iv;
            tk._weight[label] =  weight;
            assert(tracks.from_key(t->key()).get_weight(label) == weight);
            assert(tracks.from_key(t->key()).get_recv(label) == iv);
          }
	}
      }
      iv++;
    }
  }
}
/********************************************************************************************************/






/************************************************************************************************************/
void TestAnalyzer::getSimEvents_pu(PileupSummaryInfo& puInfo, std::vector<SimEvent>& simEvents)
/********************************************************************************************************
 * append the pile-up vertices to the simEvent list 
 * hard scatter assumed to be already filled in the first entry
 * Remove or ignore this for the moment, as we only will use Signal-events
 ************************************************************************************************************/
{
    report_counted("adding simEvts from puInfo", 1);
    if(simEvents.size() != 1 ){
      cout << "getSimEvents_pu : Warning !!!!  the size of simEvents is " <<  simEvents.size()  << " instead of 1" << endl;
    }

    for (unsigned int i = 0; i < puInfo.getPU_zpositions().size(); i++) {
      double t = 0;
      if (puInfo.has_times()) {
	t = puInfo.getPU_times()[i];  // appears to be in ns already
      }
      if (veryverbose_) {
	cout << setw(4) << i << "  z=" << setw(8) << fixed << setprecision(4) << puInfo.getPU_zpositions()[i]
	     << "  t=" << setw(7) << fixed << setprecision(4 ) << t 
	     << "  pthat= " << scientific << puInfo.getPU_pT_hats()[i] 
	     << " sumpt hi" << puInfo.getPU_sumpT_highpT()[i]                     // 0
	     << " sumpt lo" << puInfo.getPU_sumpT_lowpT()[i]                      // 0
	     << " ntrk hi" << fixed << puInfo.getPU_ntrks_highpT()[i]             // nonsense
	     << " ntrk lo" << fixed << puInfo.getPU_ntrks_lowpT()[i]              // nonsense
	     << endl;
      }

      SimEvent e(i+1);
      e.z = puInfo.getPU_zpositions()[i];
      e.x = vertexBeamSpot_.x(e.z);
      e.y = vertexBeamSpot_.y(e.z);
      e.t = t;

      e.type = FROM_PU_SUMMARY;  // partial
      e.pt_hat = puInfo.getPU_pT_hats()[i];

      simEvents.push_back(e);
    }
}
/********************************************************************************************************/






/********************************************************************************************************/
bool TestAnalyzer::get_MC_truth(const edm::Event& iEvent,
                                            Tracks& tracks,
                                            bool bPuInfo,
                                            PileupSummaryInfo& puInfo,
                                            std::vector<SimEvent>& simEvt)
/********************************************************************************************************/
{
  //  std::string mcproduct = "generator";       // starting with 3_1_0 pre something
  tracking_truth_available_ = false;         // meaning trackingparticles

  // genParticles for AOD et al:  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate
  Handle<reco::GenParticleCollection> genParticlesH;
  bool bgenParticles = iEvent.getByToken(genParticleCollection_Token_, genParticlesH);
  if (bgenParticles) report_counted("found genParticles", 1);

  Handle<HepMCProduct> evtMC;              // generator level
  Handle<SimVertexContainer> simVtxs;
  Handle<SimTrackContainer> simTrks;
  bool bSimVtxs = false;
  bool bSimTrks = false;
  try {
    bSimVtxs = iEvent.getByToken(edmSimVertexContainerToken_, simVtxs);
    bSimTrks = iEvent.getByToken(edmSimTrackContainerToken_, simTrks);
  } catch (...) {
    cout << "no sim tracks found" << endl;
    MC_ = false;
  }
  
  if (bSimTrks) report_counted("found simTracks", 1);
  if (bSimVtxs) report_counted("found simVtxs", 1);

  edm::Handle<TrackingParticleCollection> TPCollectionH;
  edm::Handle<TrackingVertexCollection> TVCollectionH;
  bool gotTP = iEvent.getByToken(trackingParticleCollectionToken_, TPCollectionH); //retrieves simulated particles
  bool gotTV = iEvent.getByToken(trackingVertexCollectionToken_, TVCollectionH);  //retrieves simulated vertices

  if(gotTP)
    {
      report_counted("found tracking particles", 1);
    }
  else
    {
      report_counted("no tracking particles found", 1);
    }      
  if(gotTV) report_counted("found tracking vertices", 1);

  MC_ |= gotTP;
  trkidx2tp_.clear();

  

  if (gotTP)
    {
      edm::Handle<reco::RecoToSimCollection> recoToSimH;
      iEvent.getByToken(recoToSimAssociationToken_, recoToSimH);
      tp_r2s_ = recoToSimH.product();
      tracking_truth_available_ = true;
      simEvt = getSimEvents_tp(TPCollectionH, tracks);
    }
  else
    {
      // look for mc info for the signal vertex

      if( bSimTrks && bSimVtxs ){
	simEvt = getSimEvents_simtrks(simTrks, simVtxs, tracks);
	//tsim = TestAnalyzer::getSimTrkParameters(simTrks, simVtxs, simUnit_);
	MC_ = true;
      }

      else if (iEvent.getByToken(edmHepMCProductToken_, evtMC)) {
	if(verbose_){
	  cout << "all we have is  hepmc" << endl;
	}
	//simEvt = getSimEvents_hepmc(evtMC); FIXME  convert obsolete getsimpv
	MC_ = true;
      }

      else if (bgenParticles){
	if(verbose_){
	  cout << "making sim events from genparticles" << endl;
	  //simEvt = getSimEvents_hepmc(evtMC); FIXME  convert obsolete getsimpv
	  MC_= false;
	}

      }

      //  throw in some pile-up
      if (bPuInfo) {
	report_counted("filling simEvents from signal MC  + puInfo", 1);
	getSimEvents_pu(puInfo, simEvt);
	if (verbose_) {
	  cout << "PileupSummaryInfo  nPU=" << puInfo.getPU_NumInteractions() << endl;
	}
      }
    }
  return MC_;
}//get_MC_truth
/***********************************************************************************/


void TestAnalyzer::beginJob() {
  matchsummaries_ = 4;   // number of match summaries to be dumped

  MC_ = false;

  reports_.clear();

  /* histogram booking */
  rootFile_->cd();
  
  // vertex collections
  for (std::vector<std::string>::const_iterator vCollection = vertexCollectionLabels_.begin();
       vCollection != vertexCollectionLabels_.end();
       vCollection++) {
    cout << "TestAnalyzer: booking histograms for collection " << *vCollection << endl;
    std::string sdir = *vCollection;
    if (sdir == "offlineSlimmedPrimaryVertices") {
      sdir = "offlinePrimaryVertices";
      cout << "histogram directory changed from offlineSlimmedPrimaryVertices to offlinePrimaryVertices " << endl;
    }
    TDirectory* dir = rootFile_->mkdir(sdir.c_str());
    dir->cd();
    std::string s = *vCollection;
    
    histograms_[s] = bookVertexHistograms(dir);
  }

  // others
  rootFile_->cd();
  bookTrackHistograms("tracks");
  bookSimPVHistograms("simpvs");
  bookEventHistograms("event");

  rootFile_->cd(); // warum 3 mal in dieser Funktion
}




void TestAnalyzer::endJob() {
  std::cout << flush;
  std::cout << "this is void TestAnalyzer::endJob() " << std::endl;


  /* write the rootfile */
    
  rootFile_->cd();

  std::cout << "Info=" << info_->String() << std::endl;
  std::cout << flush;

  if (info_->GetString().Length() > 0) {
    info_->Write("Info");
  }
  if (build_->GetString().Length() > 0) {
    build_->Write("Build");
  }
  
  
  cout << "empty histograms:" << endl;
  int nempty = 0;
  for (std::map<std::string, TH1*>::const_iterator hist = hTrk.begin(); hist != hTrk.end(); hist++) {
    if (hist->second->GetEntries() == 0) {
      cout << hist->first << " ";
      nempty++;
    }
  }
  cout << " found " << nempty << endl;

  rootFile_->Write();

  /* print out collected reports at the end of the logfile */

  std::cout << "*******************************************************" << std::endl;
  std::cout << "reports " << reports_.size() << std::endl;
  for (auto r : reports_) {
    std::cout << r << endl;
  }
  std::cout << "*******************************************************" << std::endl;

  /* print out counted messages */
  report_counted("summary", 0);

   /* poor mans profiling */
  for (auto it : pvtimers_ ){
    std::cout << setw(50) << it.first   << std::setw(15) << std::fixed << std::setprecision(3) << (it.second).duration * 1e-6<< " s" << "  counter = " << (it.second).counter <<std::endl;
  }

    
  std::cout << flush;
  std::cerr << "TestAnalyzer::endJob: done" << std::endl;
  std::cout << flush;
}

// helper functions

bool TestAnalyzer::getPuInfo(const Event& iEvent, PileupSummaryInfo& puInfo) {
  /*
    loads the puInfo for MC events and  fills pseudo "lumi" info
   */

  simPU_ = 0.;

  try {
    //get the pileup information
    Handle<vector<PileupSummaryInfo>> puInfoH;
    if (iEvent.getByToken(vecPileupSummaryInfoToken_, puInfoH)) {
      // ignore out-of-time pile-up
      for (unsigned int puInfo_ite = 0; puInfo_ite < (*puInfoH).size(); ++puInfo_ite) {
        if ((*puInfoH)[puInfo_ite].getBunchCrossing() == 0) {
          puInfo = (*puInfoH)[puInfo_ite];
          break;
        }
      }

      lumiPU_ = 1. + puInfo.getTrueNumInteractions();  // this is the expected number of interactions
      avglumiPU_ = lumiPU_;
      instBXLumi_ = lumiPU_ / sigma_pp_;
      avginstBXLumi_ = lumiPU_ / sigma_pp_;
      simPU_ = puInfo.getPU_NumInteractions();  // this is the actual number of interactions

      if (verbose_) {
        cout << "  puInfo size = " << puInfo.getPU_NumInteractions() << endl;
        cout << "  puInfo true = " << puInfo.getTrueNumInteractions() << endl;
        if (veryverbose_) {
          for (unsigned int i = 0; i < puInfo.getPU_zpositions().size(); i++) {
            cout << "pile-up " << setw(3) << i << ")"
                 << " z= " << setw(8) << setprecision(4) << fixed << puInfo.getPU_zpositions()[i];

            cout << "  pt_hat= " << setw(6) << setprecision(1) << fixed
                 << puInfo.getPU_pT_hats()[i]
                 //<<" ntrks_lowPt=" << setw(6) << setprecision(1) << fixed << puInfo.getPU_ntrks_lowpT()[i]
                 //<< " sumpT_lowPT " << setw(6) << setprecision(1) << fixed << puInfo.getPU_sumpT_lowpT()[i]
                 << endl;
          }
        }
      }

      return true;

    } else {
      return false;
    }

  } catch (edm::Exception const& ex) {
    if (!ex.alreadyPrinted()) {
      std::cout << ex.what() << " Maybe data?" << std::endl;
    }
  }
  return false;
}


std::vector<TestAnalyzer::SimPart> TestAnalyzer::getSimTrkParameters(
    edm::Handle<edm::SimTrackContainer>& simTrks, edm::Handle<edm::SimVertexContainer>& simVtcs, double simUnit) {
  std::vector<SimPart> tsim;
  if (simVtcs->begin() == simVtcs->end()) {
    if (verbose_) {
      cout << "  TestAnalyzer::getSimTrkParameters  no simvtcs" << endl;
    }
    return tsim;
  }
  if (verbose_) {
    cout << "  TestAnalyzer::getSimTrkParameters simVtcs n=" << simVtcs->size() << endl;
    cout << "  TestAnalyzer::getSimTrkParameters 1st position" << setw(8) << setprecision(4)
         << simVtcs->begin()->position() << endl;
  }
  double t0 = simVtcs->begin()->position().e();

  for (edm::SimTrackContainer::const_iterator t = simTrks->begin(); t != simTrks->end(); ++t) {
    if (t->noVertex()) {
      std::cout << "simtrk has no vertex" << std::endl;
    } else {
      // get the vertex position
      //HepLorentzVector v=(*simVtcs)[t->vertIndex()].position();
      math::XYZTLorentzVectorD v((*simVtcs)[t->vertIndex()].position().x(),
                                 (*simVtcs)[t->vertIndex()].position().y(),
                                 (*simVtcs)[t->vertIndex()].position().z(),
                                 (*simVtcs)[t->vertIndex()].position().e());
      int pdgCode = t->type();

      if (pdgCode == -99) {
        // such entries cause crashes, no idea what they are
        std::cout << "funny particle skipped  , code=" << pdgCode << std::endl;
      } else {
        double Q = 0;  //double Q=HepPDT::theTable().getParticleData(pdgCode)->charge();
        if ((pdgCode == 11) || (pdgCode == 13) || (pdgCode == 15) || (pdgCode == -211) || (pdgCode == -2212) ||
            (pdgCode == -321) || (pdgCode == -3222) || (pdgCode == 3312)) {
          Q = -1;
        } else if ((pdgCode == -11) || (pdgCode == -13) || (pdgCode == -15) || (pdgCode == 211) || (pdgCode == 2212) ||
                   (pdgCode == 321) || (pdgCode == 3222) ||(pdgCode == -3312)) {
          Q = 1;
        } else {
          //std::cout << pdgCode << " " <<std::endl;
        }
        math::XYZTLorentzVectorD p(t->momentum().x(), t->momentum().y(), t->momentum().z(), t->momentum().e());
        if ((Q != 0) && (p.pt() > 0.1) && (fabs(t->momentum().eta()) < 5.0) && fabs(v.z() * simUnit < 20) &&
            (sqrt(v.x() * v.x() + v.y() * v.y()) < 10.)) {

	  // FIXME, switch to the constructors of SimPart here
	  
          double x0 = v.x() * simUnit;
          double y0 = v.y() * simUnit;
          double z0 = v.z() * simUnit;
          double kappa = -Q * 0.002998 * fBfield_ / p.pt();
          double D0 = x0 * sin(p.phi()) - y0 * cos(p.phi()) - 0.5 * kappa * (x0 * x0 + y0 * y0);
          double q = sqrt(1. - 2. * kappa * D0);
          double s0 = (x0 * cos(p.phi()) + y0 * sin(p.phi())) / q;
          double s1;
          if (fabs(kappa * s0) > 0.001) {
            s1 = asin(kappa * s0) / kappa;
          } else {
            double ks02 = (kappa * s0) * (kappa * s0);
            s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
          }
          SimPart sp;  //ParameterVector par;
	  sp.pt = 0;
          sp.par[reco::TrackBase::i_qoverp] = Q / p.P();
          sp.par[reco::TrackBase::i_lambda] = M_PI / 2. - p.theta();
          sp.par[reco::TrackBase::i_phi] = p.phi() - asin(kappa * s0);
          sp.par[reco::TrackBase::i_dxy] = -2. * D0 / (1. + q);
          sp.par[reco::TrackBase::i_dsz] = z0 * sin(p.theta()) - s1 * cos(p.theta());

          sp.pdg = pdgCode;
          if (v.t() - t0 < 1e-15) {
            sp.type = 0;  // primary
          } else {
            sp.type = 1;  //secondary
          }

          // now get zpca  (get perigee wrt beam)
          //double x1 = x0; double y1 = y0;
          double x1 = x0 - vertexBeamSpot_.x(z0);
          double y1 = y0 - vertexBeamSpot_.y(z0);

          D0 = x1 * sin(p.phi()) - y1 * cos(p.phi()) - 0.5 * kappa * (x1 * x1 + y1 * y1);
          q = sqrt(1. - 2. * kappa * D0);
          s0 = (x1 * cos(p.phi()) + y1 * sin(p.phi())) / q;
          if (fabs(kappa * s0) > 0.001) {
            s1 = asin(kappa * s0) / kappa;
          } else {
            double ks02 = (kappa * s0) * (kappa * s0);
            s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
          }
          sp.ddcap = -2. * D0 / (1. + q);
          sp.zdcap = z0 - s1 / tan(p.theta());
          sp.zvtx = z0;
          sp.xvtx = x0;
          sp.yvtx = y0;
	  sp.tvtx = 0;
	  sp.charge = Q;
	  sp.ldec =0 ;

          sp.phi = p.phi();
          sp.eta = p.eta();
	  

          tsim.push_back(sp);
        }
      }
    }  // has vertex
  }    //for loop
  return tsim;
}

std::vector<TestAnalyzer::SimPart> TestAnalyzer::getSimTrkParameters(
    const Handle<reco::GenParticleCollection> genParticles) {
  std::vector<SimPart> tsim;
  double xp = 0, yp = 0, zp = -99;
  for (size_t i = 0; i < genParticles->size(); ++i) {
    const GenParticle& gp = (*genParticles)[i];
    int pdgCode = gp.pdgId();
    int st = gp.status();

    if ((st == 1) && (xp == 0) && (yp == 0) && (zp == -99)) {
      xp = gp.vx();
      yp = gp.vy();
      zp = gp.vz();
    }
    if (pdgCode == -99) {
      // such entries cause crashes, no idea what they are
      std::cout << "funny particle skipped  , code=" << pdgCode << std::endl;
    } else {
      double Q = gp.charge();
      if ((st == 1) && (Q != 0) && (gp.pt() > 0.1) && (fabs(gp.eta()) < 5.0) && fabs(gp.vz() < 20) &&
          (sqrt(gp.vx() * gp.vx() + gp.vy() * gp.vy()) < 10.)) {
        double x0 = gp.vx();
        double y0 = gp.vy();
        double z0 = gp.vz();
        double kappa = -Q * 0.002998 * fBfield_ / gp.pt();
        double D0 = x0 * sin(gp.phi()) - y0 * cos(gp.phi()) - 0.5 * kappa * (x0 * x0 + y0 * y0);
        double q = sqrt(1. - 2. * kappa * D0);
        double s0 = (x0 * cos(gp.phi()) + y0 * sin(gp.phi())) / q;
        double s1;
        if (fabs(kappa * s0) > 0.001) {
          s1 = asin(kappa * s0) / kappa;
        } else {
          double ks02 = (kappa * s0) * (kappa * s0);
          s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
        }
        SimPart sp;  //ParameterVector par;
        sp.par[reco::TrackBase::i_qoverp] = Q / gp.p();
        sp.par[reco::TrackBase::i_lambda] = M_PI / 2. - gp.theta();
        sp.par[reco::TrackBase::i_phi] = gp.phi() - asin(kappa * s0);
        sp.par[reco::TrackBase::i_dxy] = -2. * D0 / (1. + q);
        sp.par[reco::TrackBase::i_dsz] = z0 * sin(gp.theta()) - s1 * cos(gp.theta());

        sp.pdg = pdgCode;
        double t = sqrt(pow(gp.vx() - xp, 2) + pow(gp.vy() - yp, 2) + pow(gp.vz() - zp, 2));
        if (t < 1e-6) {
          sp.type = 0;  // primary
        } else {
          sp.type = 1;  //secondary
        }

        // now get zpca  (get perigee wrt beam)
        double x1 = x0;
        double y1 = y0;
        D0 = x1 * sin(gp.phi()) - y1 * cos(gp.phi()) - 0.5 * kappa * (x1 * x1 + y1 * y1);
        q = sqrt(1. - 2. * kappa * D0);
        s0 = (x1 * cos(gp.phi()) + y1 * sin(gp.phi())) / q;
        if (fabs(kappa * s0) > 0.001) {
          s1 = asin(kappa * s0) / kappa;
        } else {
          double ks02 = (kappa * s0) * (kappa * s0);
          s1 = s0 * (1. + ks02 / 6. + 3. / 40. * ks02 * ks02 + 5. / 112. * pow(ks02, 3));
        }
        sp.ddcap = -2. * D0 / (1. + q);
        sp.zdcap = z0 - s1 / tan(gp.theta());
        sp.zvtx = z0;
        sp.xvtx = x0;
        sp.yvtx = y0;

        sp.phi = gp.phi();
        sp.eta = gp.eta();

        tsim.push_back(sp);
      }
    }
  }  //for loop
  return tsim;
}







bool TestAnalyzer::match(const ParameterVector& a, const ParameterVector& b) {
  double dqoverp = a(0) - b(0);
  double dlambda = a(1) - b(1);
  double dphi = a(2) - b(2);
  double dsz = a(4) - b(4);
  if (dphi > M_PI) {
    dphi -= M_2_PI;
  } else if (dphi < -M_PI) {
    dphi += M_2_PI;
  }
  //  return ( (fabs(dqoverp)<0.2) && (fabs(dlambda)<0.02) && (fabs(dphi)<0.04) && (fabs(dsz)<0.1) );
  return ((fabs(dqoverp) < 0.2) && (fabs(dlambda) < 0.02) && (fabs(dphi) < 0.04) && (fabs(dsz) < 1.0));
}

bool TestAnalyzer::isResonance(const HepMC::GenParticle* p) {
  double ctau = (pdt_->particle(std::abs(p->pdg_id())))->lifetime();
  //std::cout << "isResonance   " << p->pdg_id() << " " << ctau << std::endl;
  return (ctau > 0) && (ctau < 1e-6);
}

bool TestAnalyzer::isFinalstateParticle(const HepMC::GenParticle* p) {
  return (!p->end_vertex() && p->status() == 1);
}

bool TestAnalyzer::isCharged(const HepMC::GenParticle* p) {
  const ParticleData* part = pdt_->particle(p->pdg_id());
  if (part) {
    return part->charge() != 0;
  } else {
    // the new/improved particle table doesn't know anti-particles
    return pdt_->particle(-p->pdg_id()) != 0;
  }
}

double TestAnalyzer::vertex_pxy(const reco::Vertex& v) {
  double z = v.z();
  double vxx = v.covariance(iX, iX) + pow(vertexBeamSpot_.BeamWidthX(), 2);
  double vyy = v.covariance(iY, iY) + pow(vertexBeamSpot_.BeamWidthY(), 2);
  ;
  double vxy = v.covariance(iX, iY);
  double dx = v.x() - vertexBeamSpot_.x(z) ;
  double dy = v.y() - vertexBeamSpot_.y(z) ;
  double D = vxx * vyy - vxy * vxy;
  double c2xy = pow(dx, 2) * vyy / D + pow(dy, 2) * vxx / D - 2 * dx * dy * vxy / D;
  return TMath::Prob(c2xy, 2);
}

double TestAnalyzer::vertex_sumw(const reco::Vertex& v) {
  double wsum = 0.;
  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    wsum += v.trackWeight(*t);
  }
  return wsum;
}

double TestAnalyzer::vertex_aptsum(const reco::Vertex& v) {
  double aptsum = 0.;
  double waptsum = 0.;
  double wsum = 0.;
  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    aptsum += t->get()->pt();
    double w = v.trackWeight(*t);
    wsum += w;
    waptsum += t->get()->pt() * w;
  }
  return aptsum;
}


bool  TestAnalyzer::uv_analysis(const MVertex & vtx, const SimEvent & simevt, double & du, double & dv, double & uError, double & vError, double &pol){
  du = 0;
  dv = 0;
  uError = 1.;
  vError = 1.;
  double pt_max = 0, px_max = 0, py_max = 0;
  double sumpxy = 0, sumpxx=0, sumpyy = 0;
  for (auto tk : vtx.tracks){
    double px = tk->px();
    double py = tk->py();
    sumpxy += px * py;
    sumpxx += px * px;
    sumpyy += py * py;
    if(tk->pt() > pt_max){
      pt_max = tk->pt();
      px_max = tk->px();
      py_max = tk->py();
    }
  }
  double nx = 0, ny = 0; 
  if(sumpxx != sumpyy){
    double tan2phi = 2 * sumpxy / (sumpxx - sumpyy);
    double phi0 = 0.5 * atan(tan2phi);
    double T0 = sumpyy*cos(phi0)*cos(phi0) - 2*sumpxy * cos(phi0) * sin(phi0) + sumpxx * sin(phi0) * sin(phi0);
    double phi1 = phi0 + 0.5 * 3.14159;
    double T1 = sumpyy*cos(phi1)*cos(phi1) - 2*sumpxy * cos(phi1) * sin(phi1) + sumpxx * sin(phi1) * sin(phi1);
    // make (nx, ny) point into the thrust direction (minimizes perpendicular momentum**2 sum)
    // a small T0/T1 should correspond to a pencil-like track list, T0 ~ T1 more spherical
    if(T0 < T1){
      nx = cos(phi0);
      ny = sin(phi0);
      pol = T0/T1;
    }else{
      nx = cos(phi1);
      ny = sin(phi1);
      pol = T1/T0;
    }
  }else{
    if (pt_max == 0) {
      return false;
    }
    nx = px_max / pt_max;
    ny = py_max / pt_max;
  }

  // (nx,ny) = thrust direction, minimizes the sum of perpendicular momentum**2 sum
  // u = (x,y) component along the thrust axis
  // v = (x,y) component perpendicular to the thrust axis
  double vxx = vtx.covariance(iX, iX);
  double vyy = vtx.covariance(iY, iY);
  double vxy = vtx.covariance(iX, iY);
  double dx = vtx.x() - simevt.x;
  double dy = vtx.y() - simevt.y;
  du = dx * nx + dy * ny; 
  dv = dx * ny - dy * nx;
  uError = sqrt(nx * nx * vxx + ny * ny * vyy + 2 * nx * ny * vxy);
  vError = sqrt(ny * ny * vxx + nx * nx * vyy - 2 * nx * ny * vxy);
  return true;
}




double TestAnalyzer::vertex_ptmax2(const reco::Vertex& v) {
  double ptmax1 = 0;
  double ptmax2 = 0;

  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    if (v.trackWeight(*t) > min_trk_in_vtx_weight_) {
      double pt = t->get()->pt();
      if (pt > ptmax1) {
        ptmax2 = ptmax1;
        ptmax1 = pt;
      } else if (pt > ptmax2) {
        ptmax2 = pt;
      }
    }
  }
  return ptmax2;
}



double TestAnalyzer::vertex_sumpt2(const reco::Vertex& v) {
  double sumpt2 = 0.;

  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    // no cut on the track weight here!
    double pt = t->get()->pt();
    sumpt2 += pt * pt;
  }
  return sqrt(sumpt2);
}

double TestAnalyzer::vertex_sumpt(const reco::Vertex& v) {
  double sumpt = 0.;

  for (trackit_t t = v.tracks_begin(); t != v.tracks_end(); t++) {
    // no cut on the track weight here!
    sumpt += t->get()->pt();
  }
  return sumpt;
}

double TestAnalyzer::vertex_r(const reco::Vertex& v) {
  double z = v.z();
  double dx = v.x() - vertexBeamSpot_.x(z) ;
  double dy = v.y() - vertexBeamSpot_.y(z) ;
  return sqrt(dx * dx + dy * dy);
}

bool TestAnalyzer::select(const MVertex& v, int level) {
  /* level
     0  !isFake  && ndof>selNdof_  (default)
     1  !isFake  && ndof>selNdof_ && prob > 0.01 
     2  !isFake  && ndof>selNdof_ && prob > 0.01 && ptmax2 > 0.4
  */
  if (v.isRecoFake())
    return false;
  if ((level == 0) && (v.ndof() >= selNdof_))
    return true;
  if ((level == 1) && (v.ndof() >= selNdof_) && (v.pxy(vertexBeamSpot_) > 0.01))
    return true;
  if ((level == 2) && (v.ndof() >= selNdof_) && (v.pxy(vertexBeamSpot_) > 0.01) && (v.ptmax2() > 0.4))
    return true;
  if ((level == 3) && (v.ndof() >= selNdof_) && (v.ptmax2() < 0.4))
    return true;
  return false;
}

bool TestAnalyzer::select(const reco::Vertex& v, int level) {
  /* level
     0  !isFake  && ndof>4  (default)
     1  !isFake  && ndof>4 && prob > 0.01 
     2  !isFake  && ndof>4 && prob > 0.01 && ptmax2 > 0.4
  */
  if (v.isFake())
    return false;
  if ((level == 0) && (v.ndof() > selNdof_))
    return true;
  if ((level == 1) && (v.ndof() > selNdof_) && (vertex_pxy(v) > 0.01))
    return true;
  if ((level == 2) && (v.ndof() > selNdof_) && (vertex_pxy(v) > 0.01) && (vertex_ptmax2(v) > 0.4))
    return true;
  if ((level == 3) && (v.ndof() > selNdof_) && (vertex_ptmax2(v) < 0.4))
    return true;
  return false;
}








/****************************************************************************/
void TestAnalyzer::fillVertexHistos(std::map<std::string, TH1*>& h,
                                                const std::string& vtype,
                                                const MVertex & v,
                                                Tracks& tracks,
                                                const double deltaz,
                                                const bool verbose)
/****************************************************************************/
{
  if(v.isRecoFake()) return;
  // delta z = z_this - z_other
  timer_start("fillVertexHistos");
  Fill(h, vtype + "/ptmax2", v.ptmax2());
  double sumpt2 = v.sumpt2();
  double sumpt = v.sumpt();  // not weighted, unlike sumapt
  Fill(h, vtype + "/sumpt", sumpt);   // code duplication, FIXME
  Fill(h, vtype + "/sumpt2", sumpt2);
  Fill(h, vtype + "/logsumpt", log(sumpt)/log(10.));
  Fill(h, vtype + "/logsumpt2", log(sumpt2)/log(10.));
  Fill(h, vtype + "/sumpt2vssumpt", sumpt, sumpt2);
  if(sumpt >0){
    Fill(h, vtype + "/sumpt2oversumpt", sumpt2/sumpt);
    Fill(h, vtype + "/sumpt2oversumptvssumpt2", sumpt2, sumpt2/sumpt);
  }
  Fill(h, vtype + "/nseltrkvtx", float(v.tracksSize()));



  timer_stop("fillVertexHistos");
}
/****************************************************************************/




void TestAnalyzer::printRecVtxs(const MVertexCollection & vcollection, std::string title) {
  int ivtx = 0;

  std::cout << std::endl << title << "  nv=" << vcollection.vtxs.size() << "" << std::endl;

  for (auto v : vcollection.vtxs){
    string vtype = " recvtx";
    if (v.isRecoFake()) {
      vtype = " \033[31mfake\033[0m  ";
    } else if (v.ndof() == -3) {
      vtype = " event   ";
    }else if (!select(v)){
      vtype = " \033[31mreject\033[0m";
    }
    std::cout << "vtx " << std::setw(3) << std::setfill(' ') << ivtx++ << vtype
	      << " #trk " << std::fixed  << std::setprecision(4) << std::setw(3) << v.tracksSize()
	      << " chi2 " << std::fixed << std::setw(5) << std::setprecision(1) << v.chi2()
	      << " sumw " << std::fixed << std::setw(6) << std::setprecision(2) << v.sumw()
	      << " ndof " << std::fixed << std::setw(6) << std::setprecision(2) << v.ndof()
	      << "  x=" << std::setw(7) << std::fixed << std::setprecision(4)
	      << v.x() << "+/-"  << std::setw(6) << v.xError()                                                // <<  std::endl
              << "  y=" << std::setw(7) << v.y() << "+/-" << std::setw(6) << v.yError()  //<< std::endl
              << "  z=" << std::setw(8) << v.z() << "+/-" << std::setw(6) << v.zError();
    std::cout << "  dxy= ("
              << std::setw(7) << std::fixed << std::setprecision(4) << v.x() - vertexBeamSpot_.x(v.z()) << ","
              << std::setw(7) << std::fixed << std::setprecision(4) << v.y() - vertexBeamSpot_.y(v.z()) << ")"
              << " pxy= " << std::setw(6) << std::fixed << std::setprecision(4) << v.pxy(vertexBeamSpot_)
              << " ptsum= " << std::setw(6)  << std::fixed << std::setprecision(1) << v.sumabspt()
      //<< " sumpt2= " << std::setw(6)  << std::fixed << std::setprecision(1) << v.sumpt2()
              << " ptmax2= " << std::setw(5) << std::fixed << std::setprecision(2) << v.ptmax2() << std::endl;
  }
  std::cout << std::endl;
}



void TestAnalyzer::printRecVtxs(const reco::VertexCollection* recVtxs, std::string title) {
  int ivtx = 0;
  std::cout << std::endl;
  std::cout << std::endl << title << "  nv=" << recVtxs->size() << "" << std::endl;

  for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
    string vtype = " recvtx  ";
    if (v->isFake()) {
      vtype = " fake   ";
    } else if (v->ndof() == -3) {
      vtype = " event   ";
    }
    std::cout << "vtx " << std::setw(3) << std::setfill(' ') << ivtx++ << vtype
	      << " #trk " << std::fixed  << std::setprecision(4) << std::setw(3) << v->tracksSize()
	      << " chi2 " << std::fixed << std::setw(5) << std::setprecision(1) << v->chi2()
	      << " sumw " << std::fixed << std::setw(6) << std::setprecision(2) << vertex_sumw(*v)
	      << " ndof " << std::fixed << std::setw(6) << std::setprecision(2) << v->ndof()
	      << "  x=" << std::setw(7) << std::fixed << std::setprecision(4)
	      << v->x() << " +/-"  << std::setw(6) << v->xError()                                                // <<  std::endl
              << "  y=" << std::setw(7) << v->y() << " +/-" << std::setw(6) << v->yError()  //<< std::endl
              << "  z=" << std::setw(8) << v->z() << " +/-" << std::setw(6) << v->zError();
    std::cout << "  dxy = ("
              << std::setw(7) << std::fixed << std::setprecision(4) << v->x() - vertexBeamSpot_.x(v->z()) << ","
              << std::setw(7) << std::fixed << std::setprecision(4) << v->y() - vertexBeamSpot_.y(v->z()) << ")"
              << " pxy= " << std::setw(7) << std::fixed << std::setprecision(4) << vertex_pxy(*v)
              << "  r= " << std::setw(7) << std::fixed << std::setprecision(4) << vertex_r(*v)
              << " ptsum= " << std::setw(6)  << std::fixed << std::setprecision(1) << vertex_aptsum(*v)
              << " sumpt2= " << std::setw(6)  << std::fixed << std::setprecision(1) << vertex_sumpt2(*v)
              << " ptmax2= " << std::setw(5) << std::fixed << std::setprecision(2) << vertex_ptmax2(*v) << std::endl;
  }
  std::cout << std::endl;
}




void TestAnalyzer::printSimVtxs(const Handle<SimVertexContainer> simVtxs) {
  int i = 0;
  for (SimVertexContainer::const_iterator vsim = simVtxs->begin(); vsim != simVtxs->end(); ++vsim) {
    if (vsim->position().x() * vsim->position().x() + vsim->position().y() * vsim->position().y() < 1.) {
      std::cout << i++ << ")" << std::scientific << " evtid=" << vsim->eventId().event() << ","
                << vsim->eventId().bunchCrossing() << " sim x=" << vsim->position().x() * simUnit_
                << " sim y=" << vsim->position().y() * simUnit_ << " sim z=" << vsim->position().z() * simUnit_
                << " sim t=" << vsim->position().t() << " parent=" << vsim->parentIndex() << std::endl;
    }
  }
}

void TestAnalyzer::printSimTrks(const Handle<SimTrackContainer> simTrks) {
  std::cout << " simTrks   type, (momentum), vertIndex, genpartIndex" << std::endl;
  int i = 1;
  for (SimTrackContainer::const_iterator t = simTrks->begin(); t != simTrks->end(); ++t) {
    std::cout << i++ << ")" << t->eventId().event() << "," << t->eventId().bunchCrossing() << (*t)
              << " index=" << (*t).genpartIndex();
    std::cout << std::endl;
  }
}

void TestAnalyzer::printRecTrks(const View<reco::Track>& recTrks) {
  cout << "printRecTrks" << endl;
  int i = 1;
  for (auto t = recTrks.begin(); t != recTrks.end(); ++t) {
    //    reco::TrackBase::ParameterVector  par = t->parameters();

    cout << endl;
    cout << "Track " << i << " ";
    i++;
    //enum TrackQuality { undefQuality=-1, loose=0, tight=1, highPurity=2, confirmed=3, goodIterative=4, qualitySize=5};
    if (t->quality(reco::TrackBase::loose)) {
      cout << "loose ";
    };
    if (t->quality(reco::TrackBase::tight)) {
      cout << "tight ";
    };
    if (t->quality(reco::TrackBase::highPurity)) {
      cout << "highPurity ";
    };
    if (t->quality(reco::TrackBase::confirmed)) {
      cout << "confirmed  ";
    };
    if (t->quality(reco::TrackBase::goodIterative)) {
      cout << "goodIterative  ";
    };
    cout << endl;

    TransientTrack tk = theB_->build(&(*t));
    tk.setBeamSpot(vertexBeamSpot_);
    double ipsig = 0;
    if (tk.stateAtBeamLine().isValid()) {
      ipsig = tk.stateAtBeamLine().transverseImpactParameter().significance();
    } else {
      ipsig = -1;
    }

    cout << Form("pt=%8.3f phi=%6.3f eta=%6.3f z=%8.4f  dz=%6.4f, ipsig=%6.1f",
                 t->pt(),
                 t->phi(),
                 t->eta(),
                 t->vz(),
                 t->dzError(),
                 ipsig)
         << endl;

    cout << Form(" found=%6d  lost=%6d   chi2/ndof=%10.3f ", t->found(), t->lost(), t->normalizedChi2()) << endl;
    cout << endl;
  }
}

/********************************************************************************************************/
// helpers for z-sorting
namespace {
  bool lt(const std::pair<double, unsigned int>& a, const std::pair<double, unsigned int>& b) {
    return a.first < b.first;
  }
}  // namespace
/********************************************************************************************************/










/********************************************************************************************************/
void TestAnalyzer::getTc(
    const std::vector<MTrack>& tracks, double& Tc, double& chsq, double& dzmax, double& dztrim) 
/********************************************************************************************************/
{
  if (tracks.size() < 2) {
    Tc = -1;
    chsq = -1;
    dzmax = -1;
    dztrim = -1;
    return;
  }

  double sumw = 0, sumwz = 0, sumww = 0, sumwwz = 0, sumwwzz = 0;
  double zmin = 1e10, zmin1 = 1e10, zmax1 = -1e10, zmax = -1e10;
  double m4 = 0, m3 = 0, m2 = 0, m1 = 0, m0 = 0;

  for(auto tk : tracks){
    double tantheta = tan(tk.theta());
    double z = tk.z();
    double dz2 = pow(tk.dzError(), 2) + pow(vertexBeamSpot_.BeamWidthX() / tantheta, 2);
    double w = 1. / dz2;  // take p=1
    sumw += w;
    sumwz += w * z;
    sumwwz += w * w * z;
    ;
    sumwwzz += w * w * z * z;
    sumww += w * w;
    m0 += w;
    m1 += w * z;
    m2 += w * z * z;
    m3 += w * z * z * z;
    m4 += w * z * z * z * z;
    if (dz2 < pow(0.1, 2)) {
      if (z < zmin1) {
        zmin1 = z;
      }
      if (z < zmin) {
        zmin1 = zmin;
        zmin = z;
      }
      if (z > zmax1) {
        zmax1 = z;
      }
      if (z > zmax) {
        zmax1 = zmax;
        zmax = z;
      }
    }
  }
  double z = sumwz / sumw;
  double a = sumwwzz - 2 * z * sumwwz + z * z * sumww;
  double b = sumw;
  if (tracks.size() > 1) {
    chsq = (m2 - m0 * z * z) / (tracks.size() - 1);
    Tc = 2. * a / b;
  } else {
    chsq = 0;
    Tc = 0;
  }
  dzmax = zmax - zmin;
  dztrim = zmax1 - zmin1;  // truncated
}//getTc
/********************************************************************************************************/




/********************************************************************************************************/
bool TestAnalyzer::select_truthMatchedTrack(const edm::RefToBase<reco::Track> track, TrackingParticleRef& tpr)const

/********************************************************************************************************/
// for a reco track select the matching tracking particle, always use this function to make sure we
// are consistent
// after calling select_truthMatchedTrack, tpr may have changed its value
// to get the TrackingParticle from the TrackingParticleRef, use ->get();
{
  if (tp_r2s_->find(track) == tp_r2s_->end()) {
    return false;
  } else {
    double f = -1e10;
    TrackingParticleRef tf;
    std::vector<std::pair<TrackingParticleRef, double>> tp = (*tp_r2s_)[track];
    for (auto it = tp.begin(); it != tp.end(); ++it) {
      if (it->second > f) {
        tf = it->first;
        f = it->second;
      }
    }
    if (f > trackAssociatorMin_) {
      tpr = tf;
      return true;
    }
  }
  return false;
}//select_truthMatchedTrack
/********************************************************************************************************/




/********************************************************************************************************/
void TestAnalyzer::printTruthMatchValues(edm::RefToBase<reco::Track> track)

/********************************************************************************************************/
// print the tracking truth assocation values for a reco track
{
  if (tp_r2s_->find(track) == tp_r2s_->end()) {
    return;
  }

  cout << "[" << setw(4) << setprecision(2) << fixed;
  std::vector<std::pair<TrackingParticleRef, double>> tp = (*tp_r2s_)[track];
  unsigned int n = 0;
  for (auto it = tp.begin(); it != tp.end(); ++it) {
    if (n > 0)
      cout << ",";
    cout << it->second;
    if (++n > 5)
      break;
  }
  cout << "]";

  return;
}//printTruthMatchValues
/********************************************************************************************************/




/********************************************************************************************************/
std::vector<edm::RefToBase<reco::Track>> TestAnalyzer::getTruthMatchedVertexTracks_obsolete(const reco::Vertex& v,
                                                                                               double min_weight) const
// for rec vertex v get a list of tracks for which truth matching is available
/********************************************************************************************************/
{
  std::vector<edm::RefToBase<reco::Track>> b;
  TrackingParticleRef tpr;

  for (trackit_t tv = v.tracks_begin(); tv != v.tracks_end(); tv++) {
    if (v.trackWeight(*tv) >= min_weight) {
      if (select_truthMatchedTrack(*tv, tpr)) {
        b.push_back(*tv);
      }
    }
  }
  return b;
}
/********************************************************************************************************/



/********************************************************************************************************/
std::vector<TestAnalyzer::SimEvent> TestAnalyzer::getSimEvents_miniaod(
											       const edm::Event & iEvent, 
											       Tracks & tracks
											       )
/********************************************************************************************************/
{
  //typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> GlobalPoint; in DataFormats/Math/interface/Point3D.h
  Handle<GenEventVertex> genParticlesXyz0Handle;
  iEvent.getByToken(theGenParticlesXyz0Token_, genParticlesXyz0Handle);
  GenEventVertex  genXyz0 = (*genParticlesXyz0Handle.product());
  Handle<float> genParticlesT0Handle;
  iEvent.getByToken(theGenParticlesT0Token_, genParticlesT0Handle);
  float t0 = (*genParticlesT0Handle.product());
  if(verbose_){
    cout << "getSimEvents_miniaod  genXyz0.z()="<< genXyz0 <<   "  t0=" << t0 << endl;
  }

  Handle<edm::View<pat::PackedGenParticle> > genParticlesHandle;
  iEvent.getByToken(theGenParticlesToken_, genParticlesHandle);
  edm::View<pat::PackedGenParticle> gen = (*genParticlesHandle.product());
  if(verbose_){
    cout << "getSimEvents_miniaod  gen.size()="<< gen.size() << endl;
  }


  // the "important" particles
  Handle<edm::View<reco::GenParticle> > prunedGenParticlesHandle;
  iEvent.getByToken(thePrunedGenParticlesToken_, prunedGenParticlesHandle);
  edm::View<reco::GenParticle> pruned = (*prunedGenParticlesHandle.product());

  std::map<const reco::GenParticle *, unsigned int> pm;
  if(verbose_){
    cout << "getSimEvents_miniaod  pruned.size()="<< pruned.size() << endl;
    /*
    if(pruned.size()>0){
      for(unsigned int p = 0; p < pruned.size(); p++){
	const reco::GenParticle * pp = &pruned[p];
	pm[pp] = p;
	cout << "pruned " << p << " pdgid " << pruned[p].pdgId() << "  vertex = " << pruned[p].vertex() << "  " << &(pruned[p]) << endl;
      }
    }
    */
  }


  /*
  edm::Handle<GenEventInfoProduct> genEvtInfoH;
  iEvent.getByLabel( "generator", genEvtInfoH );
  cout << "qScale = " << genEvtInfoH->qScale();  
  */

  vector<SimEvent> simEvents;
  SimEvent e(0);
  e.type = FROM_WHATEVER; // partial
  e.x = genXyz0.x();
  e.y = genXyz0.y();
  e.z = genXyz0.z();
  e.t = t0;

  e.pt_hat = 0;
  for(size_t it = 0; it < gen.size(); it++){
    auto & cand = gen[it];
    e.pt_hat += cand.energy();  // not really

    // get the parent vertex (is it the parent's decay or production vertex?)
    double vx = e.x;
    double vy = e.y;
    double vz = e.z;
    bool prompt = cand.isPromptFinalState() || cand.isDirectPromptTauDecayProductFinalState();
    if(prompt){
      // use the event vertex as the production point
    }
    if((cand.numberOfMothers() > 0) && (!prompt)){
      auto parent = cand.mother(0); // politically correct, parent in "pruned"
      //std::cout << "cand  " <<  cand.pdgId() << "   parent= " <<  cand.mother(0) << "  pdgid=" << parent->pdgId() << " " << parent->vertex() << " this=" << &cand << std::endl;
      if(!((parent->vx() == 0) && (parent->vy()==0) && (parent->vz()==0))){
	vx = parent->vx();
	vy = parent->vy();
	vz = parent->vz();
	if(pow(vx-e.x,2) + pow(vy-e.y,2) + pow(vz-e.z, 2) > pow(10e-4,2)){
	  prompt = false;
	}
      }
    }

      /*
      cout << "gen " << it << "  pdg " <<  cand.pdgId() <<  "(" << parentpdgid << ") " << " charge= " <<cand.charge() << " pt=" << cand.pt()
	 <<  " phi=" << cand.phi()
	 <<  " theta=" << cand.theta()
	 << " vtx " << cand.vertex() 
	 << " z " << cand.vertex().z()             // all 0
	 << " vx="  << cand.vx()
	 << " vy="  << cand.vy()
	 << " vz="  << cand.vz()
	 << " dxy="  << cand.dxy()
	 << " dz="  << cand.dzError()                  
	 << " chi2="  << cand.vertexChi2()         // always 0
      	 << " ndau"  << cand.numberOfDaughters()  // always = 0
      	 << " nmo" << cand.numberOfMothers()      //  always =1
	 << " isPromptFinalState() " << cand.isPromptFinalState()
	 << " prompt = " << prompt
	 << endl;
      */

    if ( (fabs(cand.eta()) < etaMaxVisible_) && (cand.charge() != 0) ){
      e.sumpt += cand.pt();
      e.sumpt2 += pow(cand.pt(), 2);
      if(cand.pt() > ptMinVisible_){
	e.nGenTrk ++;
	e.pxvis += cand.px();
	e.pyvis += cand.py();
	e.ptvis += cand.pt();
	// we probably can't rely on the production vertex info
	// that makes track parameter matching kind of a gamble
	auto p = SimPart(0, prompt, vx, vy, vz, cand, fBfield_);
	e.parts.push_back(p);
      }
    }
  }
  
  simEvents.push_back(e);
  return simEvents;
}//getSimEvents_miniaod
/********************************************************************************************************/



/********************************************************************************************************/
std::vector<TestAnalyzer::SimEvent> TestAnalyzer::getSimEvents_simtrks
(
 const Handle<SimTrackContainer> simTrks,
 const Handle<SimVertexContainer> simVtxs,
 Tracks& tracks
 )
// SimTrack + SimVertex info for the signal vertex only, old
// assumes that all sim tracks belong to the signal vertex
/********************************************************************************************************/{
  report_counted("TestAnalyzer::getSimEvents from SimTracks + SimVertexs",1);
  SimVertexContainer::const_iterator vsim = simVtxs->begin(); 
  vector<SimEvent> simEvt;
  SimEvent e(0); 
  e.type = FROM_HEPMC;
  e.eventId = vsim->eventId();
  e.nChTP = 0;
  e.ptvis = 0;
  e.sumpt = 0;
  e.sumpt2 = 0;
  e.pxvis = 0;
  e.pyvis = 0;
  e.sumpt2 = 0;
  e.dzmin = 999.;
  e.x = vsim->position().x() * simUnit_;
  e.y = vsim->position().y() * simUnit_;
  e.z = vsim->position().z() * simUnit_;
  e.t = vsim->position().t() * simtUnit_;

  simEvt.push_back(e);
  for (unsigned int i = 0; i < tracks.size(); i++) {
    MTrack& t = tracks[i];
    simEvt[0].trkidx.push_back(i);
    if (t.selected()) {
      simEvt[0].rtk.push_back(t);
    }

    // fill the matching info for the MTrack
    t._matched = 0;
    t._simEvt = &simEvt[0];
    t._zsim = simEvt[0].z;
    t._is_primary = true;
  }

  return simEvt;
}//getSimEvents_simtrks
/********************************************************************************************************/






/********************************************************************************************************/
std::vector<TestAnalyzer::SimEvent> TestAnalyzer::getSimEvents_tp(
    edm::Handle<TrackingParticleCollection> TPCollectionH,
    Tracks& tracks)
// collect information about simulated collisions from tracking particles / vertices
/********************************************************************************************************/
{
  const TrackingParticleCollection* simTracks = TPCollectionH.product();

  vector<SimEvent> simEvt;
  map<EncodedEventId, unsigned int> eventIdToEventMap;
  map<EncodedEventId, unsigned int>::iterator id;

  bool dumpTP = false;
  for (TrackingParticleCollection::const_iterator it = simTracks->begin(); it != simTracks->end(); it++) {
    if (fabs(it->parentVertex().get()->position().z()) > 100.)
      continue;  // skip funny entries @ -13900

    unsigned int event = 0;  //note, this is no longer the same as eventId().event()
    id = eventIdToEventMap.find(it->eventId());
    // skip out of time pile-up, irrelevant for tracking
    if (it->eventId().bunchCrossing() != 0)
      continue;
    //
    if (id == eventIdToEventMap.end()) {
      // new event here
      event = simEvt.size();
      SimEvent e(event);
      e.type = FROM_TRACKING_TRUTH;  //full
      e.eventId = it->eventId();
      e.nChTP = 0;
      e.sumpt = 0;
      e.sumpt2 = 0;
      e.dzmin = 999.;
      const TrackingVertex* parentVertex = it->parentVertex().get();
      if (DEBUG_) {
        cout << "getSimEvents_tp: ";
        cout << it->eventId().bunchCrossing() << "," << it->eventId().event() << " z=" << it->vz() << " "
             << parentVertex->eventId().bunchCrossing() << "," << parentVertex->eventId().event()
             << " z=" << parentVertex->position().z() << endl;
      }
      if (it->eventId() == parentVertex->eventId()) {
        e.x = parentVertex->position().x();
        e.y = parentVertex->position().y();
        e.z = parentVertex->position().z();
        e.t = parentVertex->position().t() * simtUnit_;
      } else {
        e.x = it->vx();
        e.y = it->vy();
        e.z = it->vz();
        e.t = 0;  // FIXME timing
      }
      simEvt.push_back(e);
      eventIdToEventMap[e.eventId] = event;
    } else {
      event = id->second;
    }

    simEvt[event].tp.push_back(&(*it));
    if ((fabs(it->eta()) < etaMaxVisible_) && (it->charge() != 0) && (it->numberOfTrackerHits() > numTrkHitsVisible_)) {
      if (it->pt() > ptMinVisible_) {
        simEvt[event].nChTP++;
	simEvt[event].ptvis += it->pt();
	simEvt[event].pxvis += it->pt() * cos(it->phi());
	simEvt[event].pyvis += it->pt() * sin(it->phi());
      }

      simEvt[event].sumpt2 += pow(it->pt(), 2);  // should keep track of decays ?
      simEvt[event].sumpt += it->pt();
    }
  }

  if (dumpTP) {
    for (auto it = simTracks->begin(); it != simTracks->end(); it++) {
      std::cout << *it << std::endl;
    }
  }

  // collect reco tracks  and truth matched tracks for each simevent

  for (unsigned int i = 0; i < tracks.size(); i++) {
    MTrack& tk = tracks[i];

    if (select_truthMatchedTrack(tracks.ref(i), tk._tpr)) {
      if (eventIdToEventMap.find(tk._tpr->eventId()) == eventIdToEventMap.end()) {
        if (tk._tpr->eventId().bunchCrossing() == 0) {
          cout << "Bug in getSimEvents, missing event ? " << tk._tpr->eventId().bunchCrossing() << ","
               << tk._tpr->eventId().event() << endl;
        } else if (DEBUG_) {
          cout << "track without SimEvent " << tk._tpr->eventId().bunchCrossing() << "," << tk._tpr->eventId().event()
               << endl;
        }
        //break;
        continue;
      }

      unsigned int event = eventIdToEventMap[tk._tpr->eventId()];
      simEvt[event].trkidx.push_back(i);
      const TrackingVertex* parentVertex = tk._tpr->parentVertex().get();
      double vx = parentVertex->position().x();  // problems with tpr->vz()
      double vy = parentVertex->position().y();
      double vz = parentVertex->position().z();
      double ipdist = sqrt(pow(simEvt[event].x - vx, 2) + pow(simEvt[event].y - vy, 2) + pow(simEvt[event].z - vz, 2)) *
                      1.e4;  // in um

      // fill the matching info for the MTrack
      tk._matched = event;
      tk._simEvt = &simEvt[event];
      tk._zsim = simEvt[event].z;
      if (ipdist < 10)  // originated less than 10 um from the interaction point
      {
        tk._is_primary = true;
      }

      if (tk.selected()) {
        simEvt[event].trkref.push_back(tracks.ref(tk.index()));  // deprecated ?
	
        simEvt[event].rtk.push_back(tk);
        if (tk.is_primary())
          simEvt[event].rtkprim.push_back(tk);
        if ((tk.ip() / tk.dip()) < 4) // should this cut be configurable? or at least a global constant?
          simEvt[event].rtkprimsel.push_back(tk);
      }
      
    } else {
      // track not truth matched
      tk._matched = NOT_MATCHED_TK_SIM;
      tk._simEvt = NULL;
      tk._zsim = 0;
      tk._tsim = 0;
      tk._is_primary = false;
    }

  }  // end of track loop


  for (unsigned int i = 0; i < simEvt.size(); i++) {
    if (simEvt[i].rtkprim.size() > 0) {
      getTc(simEvt[i].rtkprimsel, simEvt[i].Tc, simEvt[i].chisq, simEvt[i].dzmax, simEvt[i].dztrim);
      simEvt[i].zfit = -99;
    } else {
      simEvt[i].Tc = 0;
      simEvt[i].chisq = 0;
      simEvt[i].dzmax = 0;
      simEvt[i].dztrim = 0;
      simEvt[i].zfit = -99;
      simEvt[i].tfit = -99;
    }

    auto dzmin = -1.;
    for (uint vj = 0; vj < simEvt.size(); ++vj) {
      if(i == vj) continue;
      auto const dz = std::abs(simEvt[i].z - simEvt[vj].z);
      if(dzmin < 0. or dz < dzmin) dzmin = dz;
    }
    simEvt[i].dzmin = dzmin < 0. ? 999. : dzmin;

    if (DEBUG_) {
      cout << setw(3) << i << " )   nTP=" << setw(4) << simEvt[i].tp.size() << "   z=" << setw(8) << setprecision(4)
           << fixed << simEvt[i].z << "    recTrks=" << setw(3) << simEvt[i].rtk.size() << "    recTrksPrim=" << setw(3)
           << simEvt[i].rtkprim.size() << "    allTrks=" << setw(3) << simEvt[i].trkidx.size() << endl;
    }
  }

  return simEvt;
}//getSimEvents_tp
/********************************************************************************************************/








/********************************************************************************************************/
reco::VertexCollection* TestAnalyzer::vertexFilter(Handle<reco::VertexCollection> pvs, bool filter)
/********************************************************************************************************/
{
  reco::VertexCollection* pv = new reco::VertexCollection;
  if (filter) {
    // ptmin filter
    for (reco::VertexCollection::const_iterator ipv = pvs->begin(); ipv != pvs->end(); ipv++) {
      double ptmin = 0;
      for (trackit_t tv = ipv->tracks_begin(); tv != ipv->tracks_end(); tv++) {
        double pt = tv->get()->pt();
        if (pt > ptmin) {
          ptmin = pt;
        }
      }
      if (ptmin > 0.5) {
        pv->push_back(*ipv);
      }
    }
  } else {
    for (reco::VertexCollection::const_iterator ipv = pvs->begin(); ipv != pvs->end(); ipv++) {
      pv->push_back(*ipv);
    }
  }
  return pv;
}// vertexFilter
/********************************************************************************************************/







/********************************************************************************************************/
void TestAnalyzer::analyze(const Event& iEvent, const EventSetup& iSetup)
/********************************************************************************************************/
{
  eventcounter_++;
  run_ = iEvent.id().run();
  luminosityBlock_ = iEvent.luminosityBlock();
  event_ = iEvent.id().event();
  bunchCrossing_ = iEvent.bunchCrossing();
  orbitNumber_ = iEvent.orbitNumber();
  if (sigmaZoverride_ > 0)
    sigmaZ_ = sigmaZoverride_;
  MC_ = false;
  dumpThisEvent_ = false;
  forceDump_ = false;   // use with caution
  lsglb_ = 0;


  edm::Handle<std::vector<float>> extraInfoHandle;
  iEvent.getByToken(extraInfoToken_, extraInfoHandle);
  
  if(extraInfoHandle.isValid()){
    std::cout << "************************ extra ***************************" << extraInfoHandle->size()<< std::endl;
    for(auto f : *extraInfoHandle){
      std::cout << f << std::endl;
    }
  }else{
    std::cout << "************************ no extras ***************************" << std::endl;
  }


  edm::Handle<float> clusteringCPUtimeHandle;
  iEvent.getByToken(clusteringCPUtimeToken_, clusteringCPUtimeHandle);
  if (clusteringCPUtimeHandle.isValid()){
    std::cout << " clustering time is " << *clusteringCPUtimeHandle << std::endl;
  }
  
  // in case we wanted to analyze a specific lumi block only
  if ((analyzeLS_ >= 0) && !(luminosityBlock_ == analyzeLS_))
    return;
  
  // get lumi info
  if(run_ > 1) get_luminosity_infos(iEvent);

  if (verbose_) {
    std::cout << endl
              << "TestAnalyzer::analyze   event counter=" << eventcounter_ << " Run=" << run_
              << "  LumiBlock " << luminosityBlock_ << "  event  " << event_ << " bx=" << bunchCrossing_
              << "  orbit=" << orbitNumber_ << std::endl;
  }

  //Retrieve tracker topology from geometry, need for track cluster infos
  //auto const & tTopoH = iSetup.getHandle(trackerTopologyToken_);
  //tTopo_ = tTopoH.product();

  if (!do_vertex_analysis_)
    return;

  
  // load pile-up info and select events in [nPUMin_, nPUmax_] if requested
  PileupSummaryInfo puInfo;
  bool bPuInfo = getPuInfo(iEvent, puInfo);
  if (bPuInfo) {
    report_counted("found PU info", 1);
  }else{
    report_counted("no PU info found", 1);
  }   


  // select pu window if requested
  assert(puInfo.getPU_zpositions().size() == simPU_);

  if (bPuInfo && ((puInfo.getPU_zpositions().size() < nPUmin_) || (puInfo.getPU_zpositions().size() > nPUmax_))) {
    if (verbose_) {
      cout << "skipping event, out of pu window  npu=" << puInfo.getPU_zpositions().size() << endl;
    }
    return;
  }


   // get the beamspot into global "vertexBeamSpot_"
  if (!get_beamspot_data(iEvent)) {
    return;
  }


  // load vertex collections
  for (auto label : vertexCollectionLabels_){

    if (label.find("WithBS") != std::string::npos){
      selNdof_ = selNdofWithBS_;
      ndof0trk_ = 0;
    }else{
      selNdof_ = selNdofNoBS_;
      ndof0trk_ = -3.;
    }

    Handle<reco::VertexCollection> recVtxsH;
    if (iEvent.getByToken(vertexCollectionTokens_[label], recVtxsH)) {
      recVtxs_[label] = vertexFilter(recVtxsH, useVertexFilter_);
      analyzeVertexRecoCPUTime(histograms_[label], recVtxs_[label], label);
    } else {
      recVtxs_[label] = NULL;
      cout << "collection " << label << " not found " << endl;
    }
    
  }


  // load the tracks (or in case of miniaod tracks and vertices)
  Tracks tracks;
  
  //  get transient track builder and b-field
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB_);
  theB_ = iSetup.getHandle(transientTrackBuilderToken_);
  fBfield_ = ((*theB_).field()->inTesla(GlobalPoint(0., 0., 0.))).z();


  if (!get_reco_and_transient_tracks(iSetup, iEvent, tracks)) {
    
    // clean up and exit if nothing is found
    for (auto label : vertexCollectionLabels_){
      delete recVtxs_[label];
    }
    return;  //FIXME some things can be done without tracks
  }


  // collect MC information
  std::vector<SimEvent> simEvt;
  if (run_ < 1000) {
    get_particle_data_table(iSetup);
      
    MC_ = get_MC_truth(iEvent, tracks, bPuInfo, puInfo, simEvt);

  }

  fill_track_to_vertex_pointers(tracks);  // must do this here in case refitting rejects vertices



  // analyze the vertex collections
  
  for( auto label : vertexCollectionLabels_){

    if (recVtxs_[label] == NULL){
      report_counted("skipping empty vertex collection " + label,5);
      continue;
    }
    
    if (verbose_) {
      cout << "**** analyzing " << label << "  with size " << recVtxs_[label]->size() << endl;
    }

    // create MVertexCollections
    vertexes_[label] = MVertexCollection(label, recVtxs_[label], tracks);
    auto & vertexes = vertexes_[label];
    auto & histos = histograms_[label];

    // replace by set_ndof_globals(label);
    if (label.find("WithBS") !=std::string::npos){
      selNdof_ = selNdofWithBS_;
      ndof0trk_ = 0;
    }else{
      selNdof_ = selNdofNoBS_;
      ndof0trk_ = -3.;
    }

    if (MC_) {

      if(tracking_truth_available_){
	
	timer_start("tp-matching");
	tpmatch(vertexes, simEvt, tracks);    // fill the rec-to-sim matching matrix (vertex.wos[simevent])
	wos_match(vertexes, simEvt, tracks);  // make a one-to-one assignment
	timer_stop("tp-matching");

        timer_start("analyzeVertexCollectionTP");
	analyzeVertexCollectionTP(histos, vertexes, tracks, simEvt, label);
	timer_stop("analyzeVertexCollectionTP");

	analyzeVertexCollectionZmatching(histos, vertexes, simEvt, label, zwindow_sigmas_);

      }
    }
  }

  if (do_vertex_analysis_ && ((nCompareCollections_ < 0) || (eventcounter_ <= nCompareCollections_))) {
    cout << " comparing collections " << endl;
    compareCollections(simEvt);
  }

  // print summary info
  if ((dumpThisEvent_ && (dumpcounter_ < ndump_)) || (veryverbose_ && (eventcounter_ < ndump_)) ||
      (autoDumpCounter_-- > 0) || (forceDump_))
    {
      std::cout  << "  dumpThisEvent_ " << dumpThisEvent_
		 << "  dumpcounter "  << dumpcounter_ << " ndump " << ndump_
		 << "  veryverbose " << veryverbose_
		 << "  autoDumpCounter " << autoDumpCounter_
		 << "  forceDump " << forceDump_
		 << std::endl;
      dumpEventSummary(simEvt, tracks);
    }

  // clean up
  for ( auto label : vertexCollectionLabels_ ){
    delete recVtxs_[label];
  }
  
  if (verbose_) {
    std::cout << std::endl << " End of TestAnalyzer " << std::endl;
  }
}
// end of "analyze"

/***************************************************************************************/





/***************************************************************************************/
void TestAnalyzer::reportEvent(const char* message, bool dump) {
  stringstream s;
  s << "run:event:ls:bx = " << run_ << ":" << setw(10) << setfill('0') << event_ << ":" << setw(4) << setfill('0')
    << luminosityBlock_ << ":" << setw(4) << setfill('0') << bunchCrossing_ << " PU=" << setfill(' ') << setw(5)
    << setprecision(1) << fixed << lumiPU_ << " avg PU=" << setw(5) << setprecision(1) << avglumiPU_
    << " sigma_z=" << setw(6) << setprecision(2) << fixed << sigmaZ_ << "    ###   " << message;
  cout << " report  " << s.str() << endl;
  reports_.push_back("[report] " + s.str());
  dumpThisEvent_ |= dump;
}
/***************************************************************************************/


/***************************************************************************************/
void TestAnalyzer::report_counted(const std::string msg, const int max_count)
{
  if (msg == "summary")
    {
      for(auto m : counted_messages_)
	{
	  std::cout << std::setw(8) << m.second.first << " " << m.first << std::endl;
	}
      return;
    }

  report_counted(msg, "", max_count);
  /*
  if ( counted_messages_.find(msg) == counted_messages_.end() )
    {
      counted_messages_[msg] = std::make_pair<unsigned int, unsigned int>(1, max_count);
    }
  else
    {
      counted_messages_[msg].first++;
    }    

  if (counted_messages_[msg].first <= counted_messages_[msg].second)
    {
    std::cout << msg << std::endl;
    }
  */
}
/***************************************************************************************/



/***************************************************************************************/
void TestAnalyzer::report_counted(const std::string msg, const std::string msg2, const int max_count)
{
 
  if ( counted_messages_.find(msg) == counted_messages_.end() )
    {
      counted_messages_[msg] = std::make_pair<unsigned int, unsigned int>(1, max_count);
    }
  else
    {
      counted_messages_[msg].first++;
    }    

  if (counted_messages_[msg].first <= counted_messages_[msg].second)
    {
      if (msg2 == ""){
	std::cout << msg << std::endl;
      }else{
	std::cout << msg << " " << msg2 << std::endl;
      }
    }
}
/***************************************************************************************/


/***************************************************************************************/
void TestAnalyzer::reportVertex(const Vertex& v, const char* message, bool dump) {
  reportEvent(message, dump);
  double dx = v.x() - vertexBeamSpot_.x(v.z());
  double dy = v.y() - vertexBeamSpot_.y(v.z());
  std::cout << " #trk " << std::fixed << std::setprecision(4) << std::setw(3) << v.tracksSize() << " chi2 "
            << std::fixed << std::setw(5) << std::setprecision(1) << v.chi2() << " ndof " << std::fixed << std::setw(6)
            << std::setprecision(2) << v.ndof() 
	    << " x "  << std::setw(8) << std::fixed << std::setprecision(4) << v.x() << " dx " << std::setw(6) << v.xError()
            << " y "  << std::setw(8) << v.y() << " dy " << std::setw(6) << v.yError()
            << " z "  << std::setw(8) << v.z() << " dz " << std::setw(6) << v.zError() ;

  std::cout  << " r  " << std::setw(8)<< std::setprecision(4) << sqrt(dx * dx + dy * dy);
  cout << endl;
}
/***************************************************************************************/




/***************************************************************************************/
void TestAnalyzer::dumpEventSummary(std::vector<SimEvent> & simEvt, Tracks & tracks)
/***************************************************************************************/
{
// print summary info
    cout << endl
         << "Event dump " << dumpcounter_ << endl
         << "event counter=" << eventcounter_ << " Run=" << run_ << "  LumiBlock " << luminosityBlock_ << "  event  "
         << event_ << " bx=" << bunchCrossing_ << "(" << firstBXinTrain_ << ")"
         << " orbit=" << orbitNumber_ << " PU= " << lumiPU_ << " sigma_z= " << sigmaZ_
         << " zbeam = " << vertexBeamSpot_.z0() << std::endl;

    cout << "  dumpThisEvent_ = " << dumpThisEvent_ << "  dumpcounter_ = " << dumpcounter_ << "  ndump_ = " << ndump_
         << "  verbose_ = " << verbose_ << " eventcounter_ = " << eventcounter_ << "  autoDumpCounter_"
         << autoDumpCounter_ << "  forceDump_ = " << forceDump_ << endl;

    dumpcounter_++;
    forceDump_ = false;
    bool bsdumped = false;

    for ( auto label : vertexCollectionLabels_ ){

      if (recVtxs_[label] == NULL) continue;
      
      cout << " dumping collection " << label  << "   run:event =" << run_ <<":"<< event_ << endl;
      auto & vertexes = vertexes_[label];
      printRecVtxs( vertexes );
      
      // redo the matching for the dump (invalid if another collection has been analyzed)
      if(tracking_truth_available_){
	tpmatch(vertexes, simEvt, tracks);
	wos_match(vertexes, simEvt, tracks);
      }

      if (matchsummaries_ > 0)
	  printMatchingSummary(vertexes, simEvt, label);
      
      
      
      if (!bsdumped &&(dumpcounter_ < 2)) {
	cout << "beamspot " << vertexBeamSpot_ << endl;
	bsdumped = true;
      }
      
      if (verbose_)
      cout << endl << endl;
    }

    matchsummaries_--;

}//dumpEventSummary
/***************************************************************************************/




/***************************************************************************************/
void TestAnalyzer::printEventSummary_notp(MVertexCollection& vtxs,
						      Tracks& tracks,
						      vector<SimEvent>& simEvt,
						      const string message)
// make a readable summary of rec and sim vertices without tracking particle information
// only using z-ordering
/***************************************************************************************/
{
  bool show_matching = true;
  cout << endl << "printEventSummary  (notp)" << endl;
  cout << "---------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   " << message << endl;
  cout << "run " << run_ << "    event " << event_ << endl;
  cout << "---------------------------" << endl;

  vector<pair<double, unsigned int>> zrecv;
  for (unsigned int idx = 0; idx < vtxs.size(); idx++) {
    zrecv.push_back(make_pair(vtxs(idx).z(), idx));
  }
  stable_sort(zrecv.begin(), zrecv.end(), lt);

  // same for simulated vertices
  vector<pair<double, unsigned int>> zsimv;
  for (unsigned int idx = 0; idx < simEvt.size(); idx++) {
    zsimv.push_back(make_pair(simEvt[idx].z, idx));
  }
  stable_sort(zsimv.begin(), zsimv.end(), lt);

  double zmatch = 0.05; // just for printing purposes, not used for anything else
  cout.precision(4);

  cout << "   sim z";
  cout << "     pT^       ";
  cout << "     rec z";
  if(show_matching) cout << "[  w ]";
  cout   << "                                                        #trk" << endl;

  unsigned int idxrec = 0;         // move these two pointers throught the lists
  unsigned int idxsim = 0;         // to show similar (in z) vertices side-by-side
  while ((idxrec < vtxs.size()) || (idxsim < simEvt.size())) {
    string signal = " ";
    string tag = " ";
    if ((idxsim < simEvt.size()) && (zsimv[idxsim].second == 0)) {
      signal = "*";
    }
    if ((idxrec < vtxs.size()) && (zrecv[idxrec].second == 0)) {
      tag = "*";
    }

    double ndof = 0, ptmax2 = 0, wmatch=0, sumpt=0; // pxy=0
    if (idxrec < vtxs.size()) {
      ndof = vtxs(zrecv[idxrec].second).ndof();
      //pxy = vtxs(zrecv[idxrec].second).pxy(vertexBeamSpot_);
      ptmax2 = vtxs(zrecv[idxrec].second).ptmax2();
      sumpt = vtxs(zrecv[idxrec].second).sumpt();
      wmatch = vtxs(zrecv[idxrec].second).sigwntfrac;
    }

    if ((idxrec < vtxs.size()) && (idxsim < simEvt.size()) &&
        (std::abs(zrecv[idxrec].first - zsimv[idxsim].first) <
         std::max(zmatch, 3 * vtxs(zrecv[idxrec].second).zError())) &&
        (((idxsim + 1) == simEvt.size()) ||
         (fabs(zrecv[idxrec].first - zsimv[idxsim].first) < fabs(zrecv[idxrec].first - zsimv[idxsim + 1].first))) &&
        (((idxrec + 1) == vtxs.size()) ||
         (fabs(zrecv[idxrec].first - zsimv[idxsim].first) < fabs(zrecv[idxrec + 1].first - zsimv[idxsim].first)))) {
      // rec and sim  on this line
      cout << setw(8) << fixed << setprecision(4) << zsimv[idxsim].first << signal;                      // sim z
      cout << setw(6) << setprecision(1) << fixed << simEvt[zsimv[idxsim].second].pt_hat;                // sim pt_hat
      cout << "   <->    ";
      cout << setw(8) << setprecision(4) << fixed << zrecv[idxrec].first << tag;                         //rec z
      if(show_matching){
	if(wmatch > 0.01){ cout << "[" << setw(4) << setprecision(2) << wmatch << "]";}
	else{cout << "      ";}
      }
      cout << " (ndof=" << fixed << setw(6) << setprecision(1) << ndof                                   // more
	//<< ", p="  << fixed<< setw(6) << setprecision(4) << pxy 
	   << ",PT="  << fixed<< setw(6) << setprecision(1) << sumpt 
           << ", ptmax2="  << fixed<< setw(4) << setprecision(1) << ptmax2 << ")";

      if (simEvt[zsimv[idxsim].second].is_visible()) {
        cout << "            (" << fixed << setw(3) << simEvt[zsimv[idxsim].second].nGenTrk << ")";
      }

      if (zsimv[idxsim].second == 0) {
        if (tag == " ") {
          cout << "  signal vertex not tagged" << endl;
        } else {
          cout << "  signal vertex found and tagged" << endl;
        }
      } else {
        cout << endl;
      }

      idxrec++;
      idxsim++;

    } else if (((idxrec < vtxs.size()) && (idxsim < simEvt.size()) && (zrecv[idxrec].first < zsimv[idxsim].first)) ||
               ((idxrec < vtxs.size()) && (idxsim == simEvt.size()))) {
      // rec only on this line
      cout << "                         ";                                                // no sim
      cout << setw(8) << setprecision(4) << fixed << zrecv[idxrec].first << tag;           // rec z
      if(show_matching){
	if(wmatch > 0.01){ cout << "[" << setw(4) << setprecision(2) << wmatch << "]";}
	else{cout << "      ";}
      }
      cout << " (ndof=" << fixed << setw(6) << setprecision(1) << ndof;                    // more
      //cout << ", p=" << setw(6) << setprecision(4) << pxy;
      cout << ",PT="  << fixed << setw(6) << setprecision(1) << sumpt;
      cout << ", ptmax2=" << setw(4) << setprecision(1) << ptmax2 << ")";
      cout << "   fake " << endl;
      idxrec++;

    } else if (((idxrec < vtxs.size()) && (idxsim < simEvt.size()) && (zrecv[idxrec].first > zsimv[idxsim].first)) ||
               ((idxrec == vtxs.size()) && (idxsim < simEvt.size()))) {
      // sim only on this line
      cout << setw(8) << setprecision(4) <<fixed << zsimv[idxsim].first << signal;                       // sim z
      cout << setw(6) << setprecision(1) << fixed << simEvt[zsimv[idxsim].second].pt_hat;                // sim pt_hat
      cout << "          ";
      cout << "          ";
      if (simEvt[zsimv[idxsim].second].is_visible()) {
	if(show_matching) cout << "      ";
        if (zsimv[idxsim].second == 0) {
          cout << "                                       lost signal vertex" << endl;
        } else {
          cout << "                                       lost PU  ";
          if (simEvt[zsimv[idxsim].second].is_visible()) {
            cout << "(" << setw(3) << simEvt[zsimv[idxsim].second].nGenTrk << ")";
            //   << "(" << setw(5) << setprecision(3) << simEvt[zsimv[idxsim].second].pt_hat << ")"
          }
          cout << endl;
        }
      } else {
        cout << endl;
      }
      idxsim++;
    } else {
      cout << "what else?" << endl;
      break;
    }
  }
  std::cout << std::endl;
}//printEventSummary_notp
/***************************************************************************************/









/***************************************************************************************/
void TestAnalyzer::compareCollections(vector<SimEvent>& simEvt) {
  /***************************************************************************************
prints a side-by-side comparison of the selected vertex collections
simulated vertices are shown if they are available
***************************************************************************************/

  bool debug = false;

  unsigned int ncoll = vertexCollectionLabels_.size();

  cout << "----compareCollections------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   ";
  cout << "run " << run_ << "    event " << event_ << endl;
  cout << "----------------------------------------------" << endl;

  cout << setw(15) << " simulated [ ch, tk]";

  for (unsigned int icol = 0; icol < ncoll; icol++) {
    string h = vertexCollectionLabels_[icol];
    if (h.size() > 17) {
      h = h.substr(h.size() - 17, h.size());
    } else {
      h.resize(17);
    }
    cout << setw(18) << h;
  }
  cout << endl;

  int differences = 0;

  // build the table

  vector<pair<double, unsigned int*>> row;

  // unmatched rec vertices from all collections
  for (unsigned int icol = 0; icol < ncoll; icol++) {
    auto vtxs = vertexes_[vertexCollectionLabels_[icol]];
    for (unsigned int idx = 0; idx < vtxs.size(); idx++) {
      if (vtxs(idx).isRecoFake()) continue;
      if (vtxs(idx).sim == NOT_MATCHED_VTX_SIM) {
          unsigned int* v = new unsigned int[ncoll + 1];
          for (unsigned int j = 0; j < ncoll + 1; j++) {
            v[j] = NOT_MATCHED_VTX_REC;
          }
          v[icol] = idx;
          row.push_back(make_pair(vtxs(idx).z(), v));
      }
    }
  }

  if (row.size() > 1) {
    if (debug) {
      cout << "dump    size=" << row.size() << endl;
      for (unsigned int irow = 0; irow < row.size(); irow++) {
        cout << setw(2) << irow << ")";
        cout << setw(8) << setprecision(4) << fixed << row[irow].first;
        if (row[irow].second == NULL)
          continue;
        for (unsigned int i = 0; i < ncoll + 1; i++) {
          cout << setw(6) << row[irow].second[i];
        }
        cout << endl;
      }
      cout << endl;
    }

    // join rows
    int join = 0;
    while (join >= 0) {
      if (row.size() > 1) {
        stable_sort(row.begin(), row.end());
      }

      if (debug) {
        cout << "dump before joining  size=" << row.size() << endl;
        for (unsigned int irow = 0; irow < row.size(); irow++) {
          cout << setw(2) << irow << ")";
          cout << setw(8) << setprecision(4) << fixed << row[irow].first;
          if (!(row[irow].second == NULL)) {
            for (unsigned int i = 0; i < ncoll + 1; i++) {
              cout << setw(6) << row[irow].second[i];
            }
          }
          cout << endl;
        }
        cout << endl;
      }

      double dzmin = 0.1;
      join = -1;
      for (unsigned int irow = 0; irow < row.size() - 1; irow++) {
        if ((row[irow].second == NULL) || (row[irow + 1].second == NULL))
          continue;
        if ((row[irow + 1].first - row[irow].first) < dzmin) {
          bool joinable = true;
          for (unsigned int i = 0; i < ncoll + 1; i++) {
            if ((row[irow].second[i] != NOT_MATCHED_VTX_REC) && (row[irow + 1].second[i] != NOT_MATCHED_VTX_REC))
              joinable = false;
          }
          if (joinable) {
            join = irow;
            dzmin = fabs(row[irow + 1].first - row[irow].first);
          }
        }
      }

      if (join >= 0) {
        if (debug)
          cout << "joining " << join << endl;
        if ((row[join].second == NULL) || (row[join + 1].second == NULL)) {
          cout << " bad join=" << join << endl;
        }
        if (join >= int(row.size())) {
          cout << " row pointers screwed up " << join << "   " << row.size() << endl;
        }
        //join
        for (unsigned int i = 0; i < ncoll + 1; i++) {
          if ((row[join].second[i] == NOT_MATCHED_VTX_REC) && (row[join + 1].second[i] != NOT_MATCHED_VTX_REC)) {
            row[join].second[i] = row[join + 1].second[i];
            if (i == ncoll)
              row[join].first = row[join + 1].first;
          }
        }

        //row z
        if (row[join].second[ncoll] == NOT_MATCHED_VTX_REC) {
          double zrow = 0;
          int nv = 0;
          for (unsigned int i = 0; i < ncoll; i++) {
            int iv = row[join].second[i];
            if (iv > int(vertexes_[vertexCollectionLabels_[i]].size())) {
              cout << "illegal vertex index " << iv << "    join=" << join << endl;
            }else{
	      if (iv >= 0) {
		zrow += vertexes_[vertexCollectionLabels_[i]].at(iv).z();
		nv++;
	      }
	    }
          }
          if (nv > 0) {
            row[join].first = zrow / nv;
          } else {
            // hmmm
          }
        }
        //delete swallowed row
        if (debug)
          cout << "deleting row " << join + 1 << "  row size= " << row.size() << "  ncoll= " << ncoll << endl;
        delete[] row[join + 1].second;
        row[join + 1].second = NULL;
        row.erase(row.begin() + (join + 1));

        if (debug) {
          cout << "dump after joining  " << join << " with " << join + 1 << endl;
          for (unsigned int irow = 0; irow < row.size(); irow++) {
            cout << setw(2) << irow << ")";
            cout << setw(8) << setprecision(4) << fixed << row[irow].first;
            if (!(row[irow].second == NULL)) {
              for (unsigned int i = 0; i < ncoll + 1; i++) {
                cout << setw(6) << row[irow].second[i];
              }
            }
            cout << endl;
          }
          cout << endl;
        }
      }
    }

    if (debug) {
      cout << "dump after joining  size=" << row.size() << endl;
      for (unsigned int irow = 0; irow < row.size(); irow++) {
        cout << setw(2) << irow << ")";
        cout << setw(8) << setprecision(4) << fixed << row[irow].first;
        if (!(row[irow].second == NULL)) {
          for (unsigned int i = 0; i < ncoll + 1; i++) {
            cout << setw(6) << row[irow].second[i];
          }
        }
        cout << endl;
      }
      cout << endl;
    }

  }  // handle un-matched vertices and simpv's

  // fill in sim vertices and matched rec vertices
  unsigned int suppressed = 0;
  for (unsigned int idx = 0; idx < simEvt.size(); idx++) {
    if (simEvt[idx].nChTP > 0) {
      unsigned int* v = new unsigned int[ncoll + 1];
      for (unsigned int j = 0; j < ncoll + 1; j++) {
        v[j] = NOT_MATCHED_VTX_REC;
      }
      for (unsigned int jcol = 0; jcol < ncoll; jcol++) {
        auto vtxs = vertexes_[vertexCollectionLabels_[jcol]];
	int i = NOT_MATCHED_VTX_REC;
	for (unsigned int j = 0; j < vtxs.size(); j++) {
	  if (vertexes_[vertexCollectionLabels_[jcol]][j].sim == idx) {
	    i = j;
	  }
	}
	v[jcol] = i;
      }
      v[ncoll] = idx;  // the sim vertex
      row.push_back(make_pair(simEvt[idx].z, v));
    } else {
      suppressed++;
    }
  }

  if (row.size() > 1) {
    stable_sort(row.begin(), row.end());
  }

  if (debug) {
    cout << "final dump  size=" << row.size() << endl;
    for (unsigned int irow = 0; irow < row.size(); irow++) {
      cout << setw(2) << irow << ")";
      cout << setw(8) << setprecision(4) << fixed << row[irow].first;
      if (!(row[irow].second == NULL)) {
        for (unsigned int i = 0; i < ncoll + 1; i++) {
          cout << setw(6) << row[irow].second[i];
        }
      }
      cout << endl;
    }
    cout << endl;
    cout << "done" << endl;
  }

  // readable dump
  for (unsigned int irow = 0; irow < row.size(); irow++) {
    if (row[irow].second == NULL)
      continue;

    double z = row[irow].first;

    unsigned int* v = row[irow].second;
    unsigned int idx0 = v[ncoll];

    if (idx0 == NOT_MATCHED_VTX_REC) {
      cout << "%                    ";
    } else {
      // sim vertex
      cout << fixed << setw(10) << setprecision(4) << z << " [" << setw(3) << simEvt[idx0].nChTP << "," << setw(3)
	   << simEvt[idx0].rtk.size() << "]";
      if (idx0 == 0) {
	cout << "*";
      } else {
	cout << " ";
      }
    }

    // count collections that  have a rec vertex for this sim vertex (for reporting)
    unsigned int nrec = 0;
    for (unsigned int jcol = 0; jcol < ncoll; jcol++) {
      if (v[jcol] != NOT_MATCHED_VTX_REC) {
        nrec++;
      }
    }
    if ((nrec > 0) && (nrec < ncoll)) {
      differences++;
    }

    for (unsigned int jcol = 0; jcol < ncoll; jcol++) {
      unsigned int idx = v[jcol];
      if (idx != NOT_MATCHED_VTX_REC) {
        auto vtxs = vertexes_[vertexCollectionLabels_[jcol]];
        if (idx < vtxs.size()){
          cout << fixed << setw(10) << setprecision(4) << vtxs(idx).z() << " (" << setw(5) << setprecision(1)
               << vtxs(idx).ndof() << ")";
        } else {
          cout << "                  ";
        }
      } else {
        if (idx0 == NOT_MATCHED_VTX_REC) {
          // no sim vertex, "fake not not found" (found by others)
          cout << "       +          ";
        } else {
          if (nrec == 0) {
            cout << "       -          ";
          } else {
            cout << "      ---         ";  // missed one (found by others)
          }
        }
      }
    }
    cout << endl;
  }

  for (unsigned int irow = 0; irow < row.size(); irow++) {
    if (!(row[irow].second == NULL)) {
      delete[] row[irow].second;
      row[irow].second = NULL;
    }
  }

  if (differences > 0) {
    cout << "the collections differ,  " << differences << "  differences " << endl;
  }
  if (suppressed > 0) {
    cout << suppressed << "  sim vertices without tracker hits suppressed " << endl;
  }
}
/***************************************************************************************/








/*****************************************************************************************************/
std::string TestAnalyzer::formatMatchList(const std::map<unsigned int, double>& valuemap,
                                                      unsigned int nfield,
                                                      bool sim) 
/*****************************************************************************************************/
{
  // sort
  std::vector<std::pair<double, unsigned int>> values;
  for (auto it : valuemap) {
    values.push_back(make_pair(it.second, it.first));
  }
  stable_sort(values.rbegin(), values.rend());  // reverse

  std::ostringstream s;
  unsigned int n = 0;
  while ((n < nfield) && (n < values.size())) {
    if (n > 0) {
      s << ",";
    }

    s << setw(5) << setprecision(1) << fixed << values[n].first;
    if (sim) {
      s << "(" << setw(3) << setprecision(0) << fixed << values[n].second << ")";
    } else {
      s << "[" << setw(3) << setprecision(0) << fixed << values[n].second << "]";
    }
    n++;
  }

  if (n < values.size()) {
    s << "..";
  } else {
    while (n < nfield) {
      s << "           ";
      n++;
    }
    s << "  ";
  }

  return s.str();
}
/*************formatMatchList**************************************************************************/






/***************************************************************************************/
void TestAnalyzer::printMatchingSummary(MVertexCollection & vertexes,
                                                    vector<SimEvent>& simEvt,
                                                    const string message)
// dump details about the matching, tracking particles only!
/***************************************************************************************/
{
  cout << endl << "printMatchingSummary  (simEvt)" << endl;
  cout << "---------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   " << message << endl;
  cout << "run " << run_ << "    event " << event_ << endl;

  cout << "     z        t   [sim]ntrk   wnt(rec)                             z        t   (rec) ndof   wnt [sim]" << endl;
  //  dump matched and unmatched vertices here, then sort in z
  vector<pair<double, std::string>> row;
  vector<unsigned int> dumped_rec;
  
  // unmatched rec
  for (unsigned int idxrec = 0; idxrec < vertexes.size(); idxrec++) {
    if (vertexes(idxrec).isRecoFake()) continue;
    if (vertexes(idxrec).sim == NOT_MATCHED_VTX_SIM) {
      std::ostringstream s;
      s << "%" << setw(61) << " " << setw(9) << setprecision(4) << fixed << vertexes(idxrec).z() << "," << setw(7)
        << setprecision(3) << fixed << vertexes(idxrec).t() << " " << setw(1) << "(" << setw(3) << idxrec << ") "
	<< setprecision(1) << fixed << setw(5) << vertexes(idxrec).ndof() << " "
        << formatMatchList(vertexes(idxrec).wnt, 3, false);
      row.push_back(make_pair(vertexes(idxrec).z(), s.str()));
      dumped_rec.push_back(idxrec);
    }
  }

  // unmatched sim
  for (unsigned int idxsim = 0; idxsim < simEvt.size(); idxsim++) {
    if ((simEvt[idxsim].is_visible()) && (simEvt[idxsim].rec == NOT_MATCHED_VTX_REC)) {
      std::ostringstream s;
      if (idxsim == 0){
	s << "*";
      }else{
	s << " ";
      }
      s << setw(8) << setprecision(4) << fixed << simEvt[idxsim].z << "," << setw(7) << setprecision(3) << fixed
        << simEvt[idxsim].t << " " << setw(1) << "[" << setw(3) << idxsim << "]"
	<< setw(4) << simEvt[idxsim].rtk.size() << " "
        << formatMatchList(simEvt[idxsim].wnt, 3, true)
	<< "%";
      if (idxsim == 0 ) {
	s << "      signal vertex not matched   ";
      }
      row.push_back(make_pair(simEvt[idxsim].z, s.str()));
    }
  }

  // matched rec + sim
  for (unsigned int idxsim = 0; idxsim < simEvt.size(); idxsim++) {
    if (simEvt[idxsim].rec != NOT_MATCHED_VTX_REC) {
      unsigned int idxrec = simEvt[idxsim].rec;
      if( vertexes(idxrec).sim != idxsim ){ cout << "crash and burn :  " << idxrec << "," << vertexes(idxrec).sim << "," << idxsim << endl;}
      assert(vertexes(idxrec).sim == idxsim);

      std::ostringstream s;
      if (idxsim == 0){
	s << "*";
      }else{
	s << " ";
      }
      s << setw(8) << setprecision(4) << fixed << simEvt[idxsim].z << "," << setw(7) << setprecision(3) << fixed
        << simEvt[idxsim].t << " " << setw(1) << "[" << setw(3) << idxsim << "]"
	<< setw(4) << simEvt[idxsim].rtk.size() << " "
        << formatMatchList(simEvt[idxsim].wnt, 3, true) << setw(9) << setprecision(4) << fixed
        << vertexes(idxrec).z() << "," << setw(7) << setprecision(3) << fixed << vertexes(idxrec).t() << setw(1)
        << "(" << setw(4) << idxrec << ") "
	<< setprecision(1) << fixed << setw(5) << vertexes(idxrec).ndof() << " "
	<< formatMatchList(vertexes(idxrec).wnt, 3, false);

      dumped_rec.push_back(idxrec);
      row.push_back(make_pair(vertexes(idxrec).z(), s.str()));
    }
  }

  // what have we missed?
  // unmatched sim
  for (unsigned int idxrec = 0; idxrec < vertexes.size(); idxrec++) {
    if (vertexes(idxrec).isRecoFake()) continue;
    if (std::find(dumped_rec.begin(), dumped_rec.end(), idxrec) == dumped_rec.end()){
      unsigned int idxsim = vertexes(idxrec).sim;
      
      cout << "missed in printMatchingSummary  "  << idxrec <<  "  " << vertexes(idxrec).z() << "  match=" << vertexes(idxrec).sim;
      std::ostringstream s;
      s << "@" << setw(61) << " " << setw(9) << setprecision(4) << fixed << vertexes(idxrec).z() << "," << setw(7)
        << setprecision(3) << fixed << vertexes(idxrec).t() << " " << setw(1) << "(" << setw(3) << idxrec << ") "
        << formatMatchList(vertexes(idxrec).wnt, 3, false);
      row.push_back(make_pair(vertexes(idxrec).z(), s.str()));
      if(idxsim != NOT_MATCHED_VTX_SIM){
	std::ostringstream s;
	s << "@" <<setw(8) << setprecision(4) << fixed << simEvt[idxsim].z << "," << setw(7) << setprecision(3) << fixed
        << simEvt[idxsim].t << " " << setw(1) << "[" << setw(3) << idxsim << "]"
	  << setw(4) << simEvt[idxsim].rtk.size() << " "
        << formatMatchList(simEvt[idxsim].wnt, 3, true);
	row.push_back(make_pair(simEvt[idxsim].z, s.str()));
	cout << s.str();
	cout << "  matched to  " << simEvt[idxsim].rec;
      }
      cout <<endl;
    }
  } 
  //cout << "printMatchingSummary rows=" << row.size() << endl;
  stable_sort(row.begin(), row.end());
  for (auto r : row) {
    cout << r.second << endl;
  }
}
/*************printMatchingSummary******************************************************/







/***************************************************************************************/
void TestAnalyzer::printEventSummary_tp(std::map<std::string, TH1*>& h,
  					    MVertexCollection & vtxs,
                                                 Tracks& tracks,
                                                 vector<SimEvent>& simEvt,
                                                 const string message)
  // make a readable summary of the vertex finding if the TrackingParticles are availabe
/***************************************************************************************/  
{
  if (simEvt.size() == 0) {
    return;
  }

  if (eventSummaryCounter_++ > nEventSummary_){
    report_counted(message + " count limit exceeded", 1);
    return;
  }
  // sort vertices in z ... for nicer printout

  vector<pair<double, unsigned int>> zrecv;
  for (unsigned int idx = 0; idx < vtxs.size(); idx++) {
    if (vtxs(idx).isRecoFake()) continue;
    zrecv.push_back(make_pair(vtxs(idx).z(), idx));
  }
  stable_sort(zrecv.begin(), zrecv.end(), lt);

  // same for simulated vertices
  vector<pair<double, unsigned int>> zsimv;
  for (unsigned int idx = 0; idx < simEvt.size(); idx++) {
    zsimv.push_back(make_pair(simEvt[idx].z, idx));
  }
  stable_sort(zsimv.begin(), zsimv.end(), lt);

  double zmin = -50;
  double zmax = 50;
  bool fillHistograms = true;

  if (message == "signal split") {
    zmin = simEvt[0].z - 1.0;
    zmax = simEvt[0].z + 1.0;
    fillHistograms = false;
  }

  cout << endl << "printEventSummary (tracking truth)" << endl;
  cout << "---------------------------" << endl;
  cout << "event counter = " << eventcounter_ << "   " << message << endl;
  cout << "run " << run_ << "    event " << event_ << endl;
  cout << "---------------------------" << endl;

  map<unsigned int, unsigned int> rvmatch;  // reco vertex matched to sim vertex  (sim to rec)
  map<unsigned int, unsigned int> svmatch;  // sim vertex matched to rec vertex  (rec to sim)
  map<unsigned int, double> nmatch;         // highest number of truth-matched tracks of ev found in a recvtx
  map<unsigned int, double> wnmatch;        // highest sum of weights of truth-matched tracks of ev found in a recvtx
  map<unsigned int, double> purity;  // highest purity of a rec vtx (i.e. highest number of tracks from the same simvtx)
  map<unsigned int, double> wpurity;  // same for the sum of weights

  for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
    svmatch[itrec->second] = NOT_MATCHED_VTX_SIM;
    purity[itrec->second] = 0.;
    wpurity[itrec->second] = 0.;
  }

  if (zrecv.size() > 20) {

    cout << " full dump dropped because we have more than 20 vertices   nrec = " << zrecv.size() << endl;
    
    
    for (vector<pair<double, unsigned int>>::iterator itsim = zsimv.begin(); itsim != zsimv.end(); itsim++) {
      // itsim->first = z of the simvx
      // itsim->second= index of the simvtx
      if ((itsim->first < zmin) || (itsim->first > zmax))
        continue;
      SimEvent* ev = &(simEvt[itsim->second]);
      rvmatch[itsim->second] = NOT_MATCHED_VTX_REC;

      nmatch[itsim->second] = 0;   // highest number of truth-matched tracks of ev found in a recvtx
      wnmatch[itsim->second] = 0;  // highest sum of weights of truth-matched tracks of ev found in a recvtx

      // compare this sim vertex to all recvertices:
      for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
        // itrec->first  = z coordinate of the recvtx
        // itrec->second = index of the recvtx
        unsigned int irecvtx = itrec->second;
	const auto v = & vtxs(irecvtx);

        // count tracks found in both, sim and rec
        double n = 0, wt = 0;
	for( auto tv : v->tracks ){
          if (ev->hasTrack(tv)) {
            n++;
            wt += v->trackWeight(tv);
          }
        }

        // consider for matching if reasonably close in z
        double deltaz = fabs(v->z() - itsim->first);
        if ((deltaz < (5 * v->zError())) && (deltaz < 0.5)) {
          // match by number of tracks
          if (n > nmatch[itsim->second]) {
            nmatch[itsim->second] = n;
            rvmatch[itsim->second] = itrec->second;
          }

          if (n > purity[itrec->second]) {
            purity[itrec->second] = n;
            svmatch[itrec->second] = itsim->second;
          }
        }

        // match by weight
        if (wt > wnmatch[itrec->second]) {
          wnmatch[itrec->second] = wt;
        }

        if (wt > wpurity[itrec->second]) {
          wpurity[itrec->second] = wt;
        }

      }  // end of reco vertex loop
    }
  } else {
    cout << " z[cm]       rec -->   ";
    cout.precision(4);
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      cout << setw(8) << fixed << itrec->first;
    }
    cout << endl;

    cout << "                        ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      cout << setw(7) << fixed << vtxs(itrec->second).tracksSize();
      if (itrec->second == 0) {
        cout << "*";
      } else {
        cout << " ";
      }
    }
    cout << "   rec tracks" << endl;

    if(tracking_truth_available_){
    
    // truthMatchedVertexTracks[irecvtx]= number of rec tracks in that vertex for which we have truth matched simtracks
    // (not necessarily all from the same simvertex)
    map<unsigned int, int> truthMatchedVertexTracks;
    
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      truthMatchedVertexTracks[itrec->second] =
	getTruthMatchedVertexTracks_obsolete(vtxs(itrec->second).recovertex(),min_trk_in_vtx_weight_).size();  // FIXME replace consistently
      cout << "DDDDD " << truthMatchedVertexTracks[itrec->second] << " =?=" << vtxs(itrec->second).countTruthMatchedTracks(min_trk_in_vtx_weight_);
    }
    
    cout << "                        ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
	continue;
      cout << setw(7) << fixed << vtxs[itrec->second].sumwnt
	   << " ";  // may contain tracks that have weight < 0.5 after the vertex fit
    }
    cout << "   truth matched " << endl;
    
    cout << "sim ------- trk  prim ----" << endl;
    
    for (vector<pair<double, unsigned int>>::iterator itsim = zsimv.begin(); itsim != zsimv.end(); itsim++) {
      // itsim->first = z of the simvx
      // itsim->second= index of the simvtx
      if ((itsim->first < zmin) || (itsim->first > zmax))
        continue;
      SimEvent* ev = &(simEvt[itsim->second]);
      rvmatch[itsim->second] = NOT_MATCHED_VTX_REC;
      
      cout.precision(4);
      if (itsim->second == 0) {
        cout << setw(8) << fixed << ev->z << ")*" << setw(5) << ev->rtk.size() << setw(5) << ev->rtkprim.size() << "  | ";
      } else {
        cout << setw(8) << fixed << ev->z << ") " << setw(5) << ev->rtk.size() << setw(5) << ev->rtkprim.size() << "  | ";
      }

      nmatch[itsim->second] = 0;   // highest number of truth-matched tracks of ev found in a recvtx
      wnmatch[itsim->second] = 0;  // highest sum of weights of truth-matched tracks of ev found in a recvtx
      double matchpurity = 0;      //,matchwpurity=0;//,matchpurity2=0;

      // compare this sim vertex to all recvertices:
      for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
        // itrec->first  = z coordinate of the recvtx
        // itrec->second = index of the recvtx
        unsigned int irecvtx = itrec->second;
        const auto v = &(vtxs(irecvtx));

        // count tracks found in both, sim and rec
        double n = 0, wt = 0;
	for( auto tv : v->tracks ){
          if (ev->hasTrack(tv)) {
            n++;
            wt += v->trackWeight(tv);
          }
        }

        if ((itrec->first >= zmin) && (itrec->first <= zmax)) {
          if (n > 0) {
            cout << setw(7) << int(n) << " ";
          } else {
            cout << "        ";
          }
        }

        // consider for matching if reasonably close in z
        double deltaz = fabs(v->z() - itsim->first);
        if ((deltaz < (5 * v->zError())) && (deltaz < 0.5)) {
          // match by number of tracks
          if (n > nmatch[itsim->second]) {
            nmatch[itsim->second] = n;
            rvmatch[itsim->second] = itrec->second;
            //matchpurity2=matchpurity;
            matchpurity = n / truthMatchedVertexTracks[itrec->second];
            //matchwpurity=wt/truthMatchedVertexTracks[itrec->second];
          }

          if (n > purity[itrec->second]) {
            purity[itrec->second] = n;
            svmatch[itrec->second] = itsim->second;
          }
        }

        // match by weight
        if (wt > wnmatch[itrec->second]) {
          wnmatch[itrec->second] = wt;
        }

        if (wt > wpurity[itrec->second]) {
          wpurity[itrec->second] = wt;
        }

      }  // end of reco vertex loop

      // summary of this sim vertex:
      cout << "  | " << setw(1) << ev->nwosmatch << "," << setw(1) << ev->nwntmatch;
      cout << " | ";
      if (nmatch[itsim->second] > 0) {
        if (matchpurity >= 0.5) {
          cout << "found  ";
        } else if (matchpurity >= 0.3) {
          cout << "ugly   ";
        } else {
          cout << "merged ";
        }
        cout << endl;
      } else {
        // sim vertex not matched to any rec vertex
        if (ev->rtk.size() == 0) {
          cout << "invisible" << endl;
        } else if (ev->rtk.size() == 1) {
          cout << "single track " << endl;
        } else if (ev->rec != NOT_MATCHED_VTX_REC) {
          cout << "poor,  quality " << ev->matchQuality << "  for zrec " << vtxs(ev->rec).z() << endl;
        } else {
          cout << "lost " << endl;
        }
      }
    }
    cout << "---------------------------" << endl;

    //  the purity of the reconstructed vertex
    cout << "               purity   ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      cout << setw(7) << fixed << purity[itrec->second] / truthMatchedVertexTracks[itrec->second];
      if (itrec->second == 0) {
        cout << "*";
      } else {
        cout << " ";
      }  // flag the tagged vertex
    }
    cout << endl;

    //  test
    cout << "               purity*  ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      if (vtxs[itrec->second].sumwnt > 0) {
        cout << setw(7) << fixed << vtxs[itrec->second].maxwnt / float(vtxs[itrec->second].sumwnt);
      } else {
        cout << "   -   ";
      }
      if (itrec->second == 0) {
        cout << "*";
      } else {
        cout << " ";
      }  // flag the tagged vertex
    }
    cout << endl;

    // rec vertex classification: the good, the bad, and the ugly
    cout << "                     |   ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;

      if ((svmatch[itrec->second] != NOT_MATCHED_VTX_SIM) && (rvmatch[svmatch[itrec->second]] == itrec->second)) {
        if ((purity[itrec->second] / truthMatchedVertexTracks[itrec->second]) >= 0.5) {
          cout << "   ok   ";
        } else {
          cout << "  ugly  ";
        }
      } else {
        if (vtxs[itrec->second].split_from() >= 0) {
          cout << " split  ";
        } else {
          cout << "  junk  ";
        }
      }
    }
    cout << endl;

  } // !tracking_truth_available_

  
    // wos matching // FIXME this is duplicated, can we re-use tpmatch here?
    cout << "        matched vtx z| ";
    for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
      if ((itrec->first < zmin) || (itrec->first > zmax))
        continue;
      if (vtxs[itrec->second].maxwos > 0) {
        cout << setw(8) << fixed << simEvt[vtxs[itrec->second].wosmatch].z;
      } else if (vtxs[itrec->second].maxwnt > 0) {
        cout << setw(8) << fixed << simEvt[vtxs[itrec->second].wntmatch].z;
      } else {
        cout << "    -    ";
      }
    }
    cout << endl;
    cout << "---------------------------" << endl;

  }  // end of the sim vs rec table dump

  // print recvertices that were not successfully matched
  cout << "list of junk vertices" << endl;
  for (vector<pair<double, unsigned int>>::iterator itrec = zrecv.begin(); itrec != zrecv.end(); itrec++) {
    if ((itrec->first < zmin) || (itrec->first > zmax))
      continue;

    unsigned int idxrec = itrec->second;

    if ((svmatch[idxrec] == NOT_MATCHED_VTX_SIM) || (rvmatch[svmatch[idxrec]] != idxrec)) {
      cout << "zrec=" << setw(9) << fixed << setprecision(4) << itrec->first << "  ";
      if (vtxs[idxrec].maxwos > 0) {
        cout << "maxwos = " << setw(10) << fixed << setprecision(1) << vtxs[idxrec].maxwos << " / " << setw(10)
             << fixed << setprecision(1) << vtxs[idxrec].sumwos << "   zwosmatch =" << setw(9) << fixed
             << setprecision(4) << simEvt[vtxs[idxrec].wosmatch].z << "  ";
      }
      if (vtxs[idxrec].maxwnt > 0) {
        cout << "maxnt = " << setw(3) << setprecision(0) << vtxs[idxrec].maxwnt << " /" << setw(3) << setprecision(0)
             << vtxs[idxrec].sumwnt << "  z=" << setw(8) << fixed << setprecision(4)
             << simEvt[vtxs[idxrec].wntmatch].z << "  ";
      }
      cout << "  fake=" << vtxs[idxrec].is_fake();
      cout << "  qual=" << vtxs[idxrec].matchQuality;
      cout << "  split=" << vtxs[idxrec].split_from();
      cout << endl;
    }
  }

  if (!fillHistograms)
    return;

}
/***************************************************************************************/








/***************************************************************************************/
void TestAnalyzer::signalvtxmatch(MVertexCollection & vertexes, std::vector<SimEvent> & simEvt)
// if we have MC info for the signal vertex, match it to the reco vertices
// as usual, assume that the first vertex in simEvt is the signal vertex
/***************************************************************************************/
{
  if (!MC_ || simEvt.size() < 1){
    return;
  }

  auto signal = simEvt.at(0);

  vector<double> score(vertexes.size());

  for( auto v : vertexes){
    unsigned int iv = v.index();
    double dz = signal.z - v.z();
    if( fabs(dz) > 1.0){
      score[iv] = 0;
      continue;
    }
    
    score[iv] = exp(-0.5 * pow(dz / v.zError(),2));
    

  }

}
/***************************************************************************************/






/***************************************************************************************/
 void TestAnalyzer::tpmatch(MVertexCollection & vtxs,
					std::vector<SimEvent>& simEvt,
					Tracks& tracks)
  // collects truth-matched track information for matching rec and sim vertices
  // two different matches are used to determine the dominant sim vertex in a recvertex
  // wos match  = sum of "weight-over-sigma**2"
  // ntrk match = number of tracks (implicit with weight > min_trk_in_vtx_weight_), weighted with pt (cut-off at 1Gev)
  // this information is filled into both, the MVertex objects and the SimEvents
  // a rec-to-sim match is made by identifying the the dominant sim vertex
  //   (wosmatch and wntmatch)
  // note: this is not a one-to-one match, a simvertex can dominate multiple recvtxs
  // simEvt.nwosmatch  and simEvt.nwtmatch keep count the number of recvertices dominated by a simvertex
/***************************************************************************************/
{

  // turn on debugging messages for one rec vertex and one sim vertex
  bool DEBUG = false;
  unsigned int ivdebug = 0;  // rec
  unsigned int ievdebug = 0; // sim
  if(DEBUG){ cout << "tpmatch  run:event = " << run_ << ":" << event_ << "  ivdebug=" << ivdebug << " ievdebug=" << ievdebug << endl;}


  // clear old sim->rec pointers, only valid for one vertex collection
  for (auto & ev : simEvt){
    ev.clear_matching_info();
  }

  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
    auto & v = vtxs[iv];
    v.clear_matching_info();
    
    if (v.isRecoFake()) continue;

    unsigned int iev = 0;  // simevent index
    for (auto & ev : simEvt){
      double evwos = 0;       // weight over sigma**2  of simEvt ev in the current recvtx
      double evwnt = 0;       // weighted number of tracks
      unsigned int evnt = 0;  // number of tracks from simEvt ev in the current recvtx

      assert(iev == ev.index);

      for (auto tk : v.tracks){

	double wt = v.trackWeight(tk);
	if(wt < min_trk_in_vtx_weight_) continue;  

        if (tk->is_matched() && (tk->_simEvt == &ev)) {
          double dz2_beam = pow(vertexBeamSpot_.BeamWidthX() * cos(tk->phi()) / tan(tk->theta()), 2) +
	    pow(vertexBeamSpot_.BeamWidthY() * sin(tk->phi()) / tan(tk->theta()), 2);
          double dz2 = pow(tk->dzError(), 2) + dz2_beam + pow(0.0020, 2); // added 20 um, some tracks have crazy small resolutions
          double wos = v.trackWeight(tk) / dz2;
          double wnt = v.trackWeight(tk) * min(tk->pt(), 1.0);  // (truncated-)pt-weighted track count (downweights pt < 1 GeV tracks)
          v.add_truthmatched_track(iev, wos, wnt);   // fill track(wos)  rec vtx <- sim vtx, increments sumnt and sumwos
          ev.addTrack(iv, wos, wnt);                 // fill track(wos)  sim vtx -> rec vtx
          evwos += wos;
          evwnt += wnt;
          evnt++;
	  if (DEBUG && (iev==ievdebug)){
	    cout << "tpmatch iev=" << iev << "  iv=" << iv << "    wos=" << wos << "   evwos= " << evwos << "   dz2= "  << dz2 << "  w=" <<  v.trackWeight(tk) << endl;
	  }
	  
        }
      }
      if(DEBUG && (iv == ivdebug) && (iev == ievdebug)){
	cout << "tpmatch rec " << iv << "  sim " << iev << "  evwos = " << evwos<< " v.maxwos=" << v.maxwos << "   evnt=" << evnt << endl;
      }
      // require 2 tracks for a wos-match
      if ((evwos > 0) && (evwos > v.maxwos) && (evnt > 1)) {
        v.wosmatch = iev;
        v.maxwos = evwos;
        v.maxwosnt = evnt;
      }

      // weighted track counting match, require at least one track
      if ((evnt > 0) && (evwnt > v.maxwnt)) {
        v.wntmatch = iev;
        v.maxwnt = evwnt;
      }

      iev++;
      
    } // end of simevent loop

    // now wosmatch  holds the index of the dominant simvertex for this recvtx and maxwos the wos value
    if(DEBUG && (iv==ivdebug)){
      cout << "tpmatch: recvtx " << ivdebug << ", v.maxwos=" << v.maxwos << "   v.wosmatch="  << v.wosmatch << endl;
    }

    if (v.maxwos > 0) {
      simEvt.at(v.wosmatch).nwosmatch++;  // count the recvertices dominated by a simvertex
      simEvt.at(v.wosmatch).wos_dominated_recv.push_back(iv);
      assert(iv < vtxs.size());
    }

    if (v.maxwnt > 0) {
      simEvt.at(v.wntmatch).nwntmatch++;  // count the recvertices dominated by a simvertex
    }

  }  // end of vertex loop


  // for comparison with simTrack like matching (signal vertex-only)
  double sumsigwntfrac = 0;
  for(auto & v : vtxs){
    if (v.wos.find(0) == v.wos.end()){
      v.sigwosfrac = 0;
    }else{
      assert(v.wos[0] == simEvt.at(0).wos[v.index()]);
      v.sigwosfrac = v.wos[0] / simEvt.at(0).sumwos;
    }

    if (v.wnt.find(0) == v.wnt.end()){
      v.sigwntfrac = 0;
    }else{
      assert(v.wos[0] == simEvt.at(0).wos[v.index()]);
      v.sigwntfrac = v.wnt[0] / simEvt.at(0).sumwnt;
      sumsigwntfrac += v.sigwntfrac; 
    }
  }

}//tpmatch
/*******************************************************************************************************************************/



/*******************************************************************************************************************************/
void TestAnalyzer::wos_match(MVertexCollection& vtxs,
                               std::vector<SimEvent>& simEvt,
                               Tracks& tracks)
  //make rec <-> sim associations using the information filled by tpmatch
  // in contrast to tpmatch, this tries to make a one-to-one match
  //
  // first make only associations when a recvertex is dominated by a simvertex that is not
  // already assigned to another rec vertex
  // dominated = either by wos or (weighted) nt
  // notes:
  //  * unless a rec vertex consists entirely of fake tracks, it must be dominated
  //    by a simvertex
  //  * a rec-vertex that does not get sim vertex in this way is either a split-vertex
  //    in the sense that the dominating sim vertex also dominates another rec vertex
  //    or it consists entirely of unmatched tracks
  // * a merged vertex would be a simvertex that shares (enough? the majority of its?) tracks with a recvertex
  //   but another sim vertex dominates that vertex
  //
  // finally match rec vertices that are not necessarily dominated by a sim vertex,
  // but the dominating sim-vertex has been matched already and there is a simvertex that
  // really wants to have that rec-vertex. This corresponds to a "small" vertex close to "big" one
  // where the small vertex is contaminated a lot by the big one, but still recognizable
  //
/*******************************************************************************************************************************/
{
  bool DEBUG = false;

  // reset
  if(DEBUG) {cout << "wos_match  vtxs.size() = " << vtxs.size()  << endl;}
  for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
    vtxs[iv].sim = NOT_MATCHED_VTX_SIM;
    vtxs[iv].matchQuality = 0;
  }
  for (unsigned int iev = 0; iev < simEvt.size(); iev++){
    simEvt[iev].rec = NOT_MATCHED_VTX_REC;
    simEvt[iev].matchQuality = 0;
  }


  if(DEBUG){
    cout << "DEBUG: wos_match nwosmatch[0] = " << simEvt.at(0).nwosmatch << endl;
  }

  // select a rec vertex for simEvt if that rec vertex is dominated by the sim vertex
  // start with assigning sim vertices that dominate exactly one rec vertex (rank 1)
  // then repeat with sim vertices that dominate more recvertices (rank 2, 3, ...)
  // when two or more rec vertices are dominated by the same sim vertex,
  // assign the rec vertex that got more from that sim (higher wos)
  for (unsigned int rank = 1; rank < 8; rank++)
    {
      for (unsigned int iev = 0; iev < simEvt.size(); iev++)
	{
	  assert(simEvt.at(iev).wos_dominated_recv.size() == simEvt.at(iev).nwosmatch);
	  if (simEvt.at(iev).nwosmatch == 0) continue;     // doesn't dominate anything
	  if (simEvt.at(iev).nwosmatch > rank) continue;   // less ambiguous first

	  // only continue for simEvts that have not already been matched 
	  if (simEvt.at(iev).rec != NOT_MATCHED_VTX_REC) continue;
	  
	  // select a rec vertex (index iv)
	  unsigned int iv = NOT_MATCHED_VTX_REC;
	  for (unsigned int k = 0; k < simEvt.at(iev).wos_dominated_recv.size(); k++)
	    {
	      unsigned int rec = simEvt.at(iev).wos_dominated_recv.at(k); // candidate (index in vtxs)
	      if (vtxs(rec).sim != NOT_MATCHED_VTX_SIM) continue; // already matched
	      if (fabs(simEvt.at(iev).z - vtxs(rec).z()) > zWosMatchMax_) continue;// insanely far away
	      if ( (iv == NOT_MATCHED_VTX_REC) || (simEvt.at(iev).wos.at(rec) > simEvt.at(iev).wos.at(iv)) )
		{
		  iv = rec;
		}
	    }
	  // if we have found a viable candidate, make the link
	  if (iv != NOT_MATCHED_VTX_REC)
	    {
	      vtxs.at(iv).sim = iev;
	      simEvt.at(iev).rec = iv;
	      vtxs.at(iv).matchQuality = rank;
	      simEvt.at(iev).matchQuality = rank;
	      if(DEBUG) {cout << "wos_match : match made  [" << iev << "]  <-> (" << iv << ")   rank=" << rank << endl;}
	    }
	}
    }
      

  // by now we have exhausted the rec vertices that are dominated by an unmatched simvertex
  // have we?
  for(unsigned int iv = 0; iv < vtxs.size(); iv++) {
    if ((vtxs.at(iv).sim == NOT_MATCHED_VTX_SIM) && (vtxs.at(iv).maxwos > 0)){
      if (simEvt.at(vtxs.at(iv).wosmatch).rec == NOT_MATCHED_VTX_REC){
	auto ev = simEvt.at(vtxs.at(iv).wosmatch);
	cout << "wos_match :  unmatched [" << vtxs.at(iv).wosmatch << "]   dominantes unmatched  ("<< iv <<")   ????"
	     << "  rec.maxwos="  <<  vtxs.at(iv).maxwos 
	     << "  rec.maxwosnt="  <<  vtxs.at(iv).maxwosnt
	     << "  rec.sumwos=" << vtxs.at(iv).sumwos 
	     << "  sim.nwosmatch = " << ev.nwosmatch
	     << "  zrec=" << vtxs.at(iv).z()
	     << "  zsim=" << ev.z
	     <<endl;
	reportEvent("wos_match :  unmatched dominates matched");
	}
      }
    }


  // give vertices a chance that have a lot of overlap, but are still recognizably
  // caused by a specific simvertex (without being classified as dominating)
  // like a small peak sitting on the flank of a larger nearby peak

  unsigned int ntry = 0;
  while (ntry++ < 10)
    {
      unsigned nmatch = 0;
      for (unsigned int sim = 0; sim < simEvt.size(); sim++)
	{
	  if ((simEvt.at(sim).rec != NOT_MATCHED_VTX_REC) || (simEvt.at(sim).wos.size() == 0))
	    continue;
	
	  // ok, single simvertex sim, who is your your favorite rec vertex?
	  unsigned int rec = NOT_MATCHED_VTX_REC;
	  for (auto rv : simEvt.at(sim).wos) {
	    if ((rec == NOT_MATCHED_VTX_REC) || (rv.second > simEvt.at(sim).wos.at(rec))) {
	      rec = rv.first;
	    }
	  }

	  if (rec == NOT_MATCHED_VTX_REC){
	    cout << "wos_match : last hope failed for (" << rec << ")" << endl;
	    // second chance if wos didn't work?
	    for (auto rv : simEvt.at(sim).wnt) {
	      if ((rec == NOT_MATCHED_VTX_REC) || (rv.second > simEvt.at(sim).wnt.at(rec))) {
		rec = rv.first;
	      }
	    }
	  }
	
	if (rec == NOT_MATCHED_VTX_REC)
	  continue;  // should not happen
	if (vtxs.at(rec).sim != NOT_MATCHED_VTX_SIM)
	  continue;  // already gone
	
	// do you, recvertex rec, take this simvertex sim as your lawful wedded truthmatch?
	unsigned int rec2sim = NOT_MATCHED_VTX_SIM;
	for (auto sv : vtxs.at(rec).wos) {
	  if (simEvt.at(sv.first).rec != NOT_MATCHED_VTX_REC)
	    continue;  // already used
	  if ((rec2sim == NOT_MATCHED_VTX_SIM) || (sv.second > vtxs.at(rec).wos.at(rec2sim))) {
	    rec2sim = sv.first;
	  }
	}
	
	if (sim == rec2sim) {
	  // I do
	  vtxs.at(rec).sim = sim;
	  vtxs.at(rec).matchQuality = 8;
	  simEvt.at(sim).rec = rec;
	  simEvt.at(sim).matchQuality = 8;
	  //cout << "wos_match : late match  [" << sim << "]   <--> (" << rec << ")" <<  "    try=" << ntry << endl;
	  nmatch++;
	}
      }  //sim loop
      // cout << "wos_match late matches " << nmatch <<  "  in try " << ntry << endl;
      if (nmatch == 0) {
	break;
      }
    }  // ntry
  
  if (DEBUG) {
    for (unsigned int iv = 0; iv < vtxs.size(); iv++) {
      if (iv == 0) {
        dumpThisEvent_ = true;
        cout << "wos_match    recvtx = " << setw(4) << iv << "   (sim=" << vtxs[iv].sim << ")" << endl;
        cout << " splitfrom= " << vtxs[iv].split_from() << endl;
        cout << " maxwos = " << vtxs[iv].maxwos << endl;
        cout << " sumwos = " << vtxs[iv].sumwos << endl;

        for (auto w = vtxs[iv].wos.begin(); w != vtxs[iv].wos.end(); w++) {
          cout << "matching (wos)  simevent = " << setw(4) << w->first << "  zsim = " << setw(8) << fixed << setprecision(4)
               << simEvt.at(w->first).z << "  wos=" << setw(10) << w->second << endl;
        }

        for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
          for (auto rv = simEvt.at(iev).wos.begin(); rv != simEvt.at(iev).wos.end(); rv++) {
            if (rv->first == iv) {
              cout << "matching (wos)  simevent " << iev << "    zsim= " << setw(8) << fixed << setprecision(4)
                   << simEvt.at(iev).z << "    zrec= " << setw(8) << fixed << setprecision(4)
                   << vtxs.at((*rv).first).z() << "  nt=" << (*rv).second << endl;
            }
          }
        }

	for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
          for (auto rv = simEvt.at(iev).wnt.begin(); rv != simEvt.at(iev).wnt.end(); rv++) {
            if (rv->first == iv) {
              cout << "matching (wnt)  simevent " << iev << "    zsim= " << setw(8) << fixed << setprecision(4)
                   << simEvt.at(iev).z << "    zrec= " << setw(8) << fixed << setprecision(4)
                   << vtxs.at((*rv).first).z() << "  nt=" << (*rv).second << endl;
            }
          }
        }
      }
    }
  }  // <<<< debugging

  /*
  // last chance for vertices that are not matched by tracks
  // but are reasonably close and not splitters
  for(unsigned int iev=0; iev<simEvt.size(); iev++)
    {
      if( simEvt[iev].rec == NOT_MATCHED )
	{
	  for(unsigned int iv=0; iv<vtxs.size(); iv++)
	    {
	      
	      if ((vtxs[iv].sim == NOT_MATCHED) && (vtxs[iv].split_from() < 0))
		{
		  double dz = simEvt.at(iev).z-vtxs.at(iv).z();
		  if( ( fabs(dz) < 0.1 )
		      ||( fabs(dz) < (2*vtxs.at(iv).zError()) )
		      )
		  {
		    vtxs[iv].matchQuality = 10;
		    vtxs[iv].sim = iev;
		    simEvt.at(iev).rec = iv;
		    simEvt.at(iev).matchQuality = 10;
		  }
		}
	    }
	}
    }
  */


  // for convenience, store some values from the matched rec vertex with the sim vertex
  for (unsigned int iev = 0; iev < simEvt.size(); iev++) {
    int iv = simEvt.at(iev).rec;
    if (iv != NOT_MATCHED_VTX_REC) {
      simEvt.at(iev).ndof = vtxs.at(iv).ndof();
      simEvt.at(iev).zrec = vtxs.at(iv).z();
    } else {
      simEvt.at(iev).ndof = 0.;
      simEvt.at(iev).zrec = 1000.;
    }
  }


  // consistency check
  for(unsigned int iev = 0; iev < simEvt.size(); iev++) {
    auto rec= simEvt.at(iev).rec;
    assert((rec == NOT_MATCHED_VTX_REC) || ((rec<vtxs.size()) && (vtxs.at(rec).sim == iev)));
  }
  for(unsigned int iv = 0; iv < vtxs.size(); iv++) {
    auto sim = vtxs.at(iv).sim;
    assert((sim == NOT_MATCHED_VTX_SIM) || ((sim < simEvt.size()) && (simEvt.at(sim).rec == iv)));
  }

} // wos_match
/******************************************************************************/









/******************************************************************************/
void TestAnalyzer::analyzeVertexTrackAssociation(std::map<std::string, TH1*>& h, MVertexCollection& vtxs, Tracks& tracks, std::vector<SimEvent> const&, float const npu){
  //
  // Track-Vertex Association: Efficiency
  //
  // - denominator:
  //   - tracks matched to a TrackingParticle (TP-matched tracks)
  // - numerators:
  //   - TP-matched tracks assigned to a recoVertex
  //   - TP-matched tracks assigned to the recoVertex matched to the TrackingVertex of the matched-TrackingParticle
  //
  for(auto const& tk : tracks){
    // skip recoTracks not matched to a TrackingParticle
    if (not tk.is_matched()) continue;

    // the matched-TrackingParticle is assigned to the signal TrackingVertex (index == 0)
    auto const isSignalSimVtx = tk._simEvt->is_signal();

    // properties of the matched-TrackingParticle
    auto const tpPt = tk._tpr->pt();
    auto const tpEta = tk._tpr->eta();

    // denominator: TP-matched tracks
    Fill(h, "trkVtxAssocEffic_TPMatchedTracks_tpPt", tpPt, isSignalSimVtx);
    Fill(h, "trkVtxAssocEffic_TPMatchedTracks_tpEta", tpEta, isSignalSimVtx);

    // find the recoVertex to which this recoTrack is assigned
    // (max trackWeight above 0.5)
    int bestRecoVtxIdx = -1;
    auto maxTrkWgt = 0.f;
    for (size_t vtxIdx=0; vtxIdx<vtxs.size(); ++vtxIdx){
      if(not select(vtxs[vtxIdx])) continue;

      if(not vtxs[vtxIdx].has_track_key(tk.key())) continue;

      auto const wt = vtxs[vtxIdx].trackWeight(tk);
      if(wt < min_trk_in_vtx_weight_) continue;

      if(wt > maxTrkWgt or bestRecoVtxIdx == -1){
        maxTrkWgt = wt;
        bestRecoVtxIdx = vtxIdx;
      }
    }

    // numerator #1: TP-matched recoTracks assigned to a(ny) recoVertex
    auto const belongsToOneRecoVtx = (bestRecoVtxIdx >= 0);
    if(not belongsToOneRecoVtx) continue;

    Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpPt", tpPt, isSignalSimVtx);
    Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithRecoVtx_tpEta", tpEta, isSignalSimVtx);

    // numerator #2: correct track-vertex association
    // - the recoVertex to which the recoTrack is assigned corresponds to
    //   the recoVertex matched to the TrackingVertex from which the matched-TrackingParticle originates

    // index of the recoVertex (if any) matched to the TrackingVertex from which the matched-TrackingParticle originates
    auto const idxOfRecoVtxOfSimVtx = tk._simEvt->rec;

    auto const belongsToCorrectRecoVtx = (uint(bestRecoVtxIdx) == idxOfRecoVtxOfSimVtx);
    if(not belongsToCorrectRecoVtx) continue;

    Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpPt", tpPt, isSignalSimVtx);
    Fill(h, "trkVtxAssocEffic_TPMatchedTracksWithCorrectRecoVtx_tpEta", tpEta, isSignalSimVtx);
  }

  //
  // Track-Vertex Association: Purity
  //
  // - denominator:
  //   - recoTracks of a recoVertex
  // - numerator:
  //   - recoTracks of a recoVertex that are TP-matched,
  //     and have the correct track-vertex association
  //     (SimVertex of TP is matched to recoVertex of recoTrack)
  //
  for (auto const& v : vtxs){
    // skip recoVertexs not passing selection criteria
    if (not select(v)) continue;

    auto const isSignalVtx = v.is_signal();

    uint ntk = 0;
    uint ntk_woFakeVtxs = 0;
    uint ntk_matched = 0;
    for (auto tv : v.tracks){
      // restrict to recoTracks truly assigned to this recoVertex (trackWeight >= min)
      auto const wt = v.trackWeight(tv);
      if (wt < min_trk_in_vtx_weight_) continue;

      // denominator: recoTracks of a recoVertex
      ++ntk;

      auto const trkPt = tv->pt();
      auto const trkEta = tv->eta();

      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkPt", trkPt, isSignalVtx);
      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta", trkEta, isSignalVtx);
      if     (trkPt <  1.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt000to001", trkEta, isSignalVtx); }
      else if(trkPt <  3.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt001to003", trkEta, isSignalVtx); }
      else if(trkPt < 10.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt003to010", trkEta, isSignalVtx); }
      else                { Fill(h, "trkVtxAssocPurity_TracksOfRecoVtx_trkEta_trkPt010toINF", trkEta, isSignalVtx); }

      // restrict to recoVertexs that are matched to a SimVertex (i.e. skip fake recoVertexs)
      auto const recoVtxHasSimMatch = (v.sim != NOT_MATCHED_VTX_SIM);
      if(not recoVtxHasSimMatch) continue;

      ++ntk_woFakeVtxs;

      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkPt", trkPt, isSignalVtx);
      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta", trkEta, isSignalVtx);
      if     (trkPt <  1.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt000to001", trkEta, isSignalVtx); }
      else if(trkPt <  3.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt001to003", trkEta, isSignalVtx); }
      else if(trkPt < 10.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt003to010", trkEta, isSignalVtx); }
      else                { Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxMatchedToSimVtx_trkEta_trkPt010toINF", trkEta, isSignalVtx); }

      // index of the SimVertex from which the TrackingParticle matched to this recoTrack originates
      // (for TP-unmatched recoTracks, this equals NOT_MATCHED_TK_SIM)
      auto const tk_sim = tracks.simevent_index_from_key(tv->key());

      // recoTrack is TP-matched and its recoVertex is matched to the SimVertex of the recoTrack's TP
      auto const hasCorrectTrkVtxAsso = (tk_sim != NOT_MATCHED_TK_SIM and tk_sim == v.sim);
      if(not hasCorrectTrkVtxAsso) continue;

      ++ntk_matched;

      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkPt", trkPt, isSignalVtx);
      Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta", trkEta, isSignalVtx);
      if     (trkPt <  1.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt000to001", trkEta, isSignalVtx); }
      else if(trkPt <  3.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt001to003", trkEta, isSignalVtx); }
      else if(trkPt < 10.){ Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt003to010", trkEta, isSignalVtx); }
      else                { Fill(h, "trkVtxAssocPurity_TracksOfRecoVtxWithCorrectSimVtxMatch_trkEta_trkPt010toINF", trkEta, isSignalVtx); }
    }

    if(ntk > 0){
      auto const purity = ntk_matched / float(ntk);
      Fill(h, "trkVtxAssocPurity_vs_vz", v.z(), purity, isSignalVtx);
      Fill(h, "trkVtxAssocPurity_vs_pu", npu, purity, isSignalVtx);
    }

    if(ntk_woFakeVtxs > 0){
      auto const purity = ntk_matched / float(ntk_woFakeVtxs);
      Fill(h, "trkVtxAssocPurityWithoutFakeRecoVtxs_vs_vz", v.z(), purity, isSignalVtx);
      Fill(h, "trkVtxAssocPurityWithoutFakeRecoVtxs_vs_pu", npu, purity, isSignalVtx);
    }
  }
}



///Helper function for analyzeVertexcollection TP 
   //first we define a function to do the calculation
    double distance_point_to_plane(const double point[3], const double plane_point[3], const double normal_vector[3]){
        // in order to make the formula easier we unpack the values from the vectors
        double x1 = point[0], y1 = point[1], z1 = point[2];

        double x0 = plane_point[0], y0 = plane_point[1], z0 = plane_point[2];
        double A = normal_vector[0], B = normal_vector[1], C = normal_vector[2];


        // now we calculate D
        double D = -(A * x0 + B * y0 + C * z0);
        // and the numerator | Ax1 + By1 + Cz1 + D |

        double numerator = fabs(A * x1 + B * y1 + C * z1 + D);
        // and denominator  square root ( A * A + B *B + C*C)
        double denominator = sqrt(A * A + B * B + C * C);
        if(denominator!= 0){
          return numerator / denominator;
        }else{
          std::cerr << " TrueE3DDistanceToPlane had a divide by zero error" << std::endl;

          return -0; // return unusal value to make
        }



    }





/***************************************************************************************/
void TestAnalyzer::analyzeVertexCollectionTP(std::map<std::string, TH1*>& h,
                                             MVertexCollection& vtxs,
                                             Tracks& tracks,
                                             vector<SimEvent>& simEvt,
                                             const string message)
// with track truthmatching (tp)  
/*********************************************************************************************/
{
    // Retrieve the histogram from the map
    TH1I* hist_reco_matches = dynamic_cast<TH1I*>(h["efficiency/SEEfficiency"]);
      if (!hist_reco_matches) {
        std::cerr << "Error: Histogram SEEfficiency not found!" << std::endl;
        return;
    }
   // if (!hist) {
    //    std::cerr << "Histogram not found!" << std::endl;
    //    return;
   // } shouldnt need that since fill should do that on its own
    std::cout << "I am in Analyzer function" << std::endl;

    double realCounter = 0, wrongCounter = 0;
    for(const auto& v : vtxs){
        if (!v.isRecoFake()) { // Check if the vertex is not a reconstruction fake
            if (v.is_real()) { // Check if the vertex is real
                realCounter++;
            } else {
                wrongCounter++;
            }
        }
    }

    std::cout << "Real and wrong reconstructed: " << realCounter << " " << wrongCounter << std::endl;

    // Fill the "Yes" bin (first bin, index 1) with realCounter
    hist_reco_matches->SetBinContent(1, realCounter);
    
    // Fill the "No" bin (second bin, index 2) with wrongCounter
    hist_reco_matches->SetBinContent(2, wrongCounter);


    std::cout << "Filling histogram: SEEfficiency" << std::endl;
    //Fill(h, "efficiency/SEEfficiency", realCounter, wrongCounter);


    // next histogram 
    TH2F* SERecoVsSimZPosition = dynamic_cast<TH2F*>(h["efficiency/SERecoVsSimZPosition"]);
      if (!SERecoVsSimZPosition) {
        std::cerr << "Error: Histogram SERecoVsSimZPosition not found!" << std::endl;
        return;
    }
    if (simEvt[0].is_matched()) {
    MVertex& matchedVtx = vtxs.at(simEvt[0].rec);
    SERecoVsSimZPosition->Fill(matchedVtx.z(), simEvt[0].z);
    }



    //new histogram 
    TH2F* PURecoVsSimZPosition = dynamic_cast<TH2F*>(h["efficiency/PURecoVsSimZPosition"]);
      if (!PURecoVsSimZPosition) {
        std::cerr << "Error: Histogram PURecoVsSimZPosition not found!" << std::endl;
        return;
    }

    for (size_t i = 0; i < simEvt.size();i++){
      if(!simEvt[i].is_signal()){
      if (simEvt[i].is_matched()) {  // Check if the simulated vertex is matched, if not we cannot take the distance

            unsigned int rec_index = simEvt[i].rec;  // Get the index of the matched reconstructed vertex
            double true_z = simEvt[i].z;
            double rec_z = vtxs.at(rec_index).z();
                PURecoVsSimZPosition->Fill(rec_z, true_z);
        }
      }

    }
      // new histogram, SE resolution

      TH1F *SEResolution = dynamic_cast<TH1F *>(h["efficiency/SEResolution"]);
    if (!SEResolution) {
        std::cerr << "Error: Histogram SEResolution not found!" << std::endl;
        return;
    }
      if (simEvt[0].is_matched()) {
      MVertex& matchedVtx = vtxs.at(simEvt[0].rec);
      SEResolution->Fill( simEvt[0].z -matchedVtx.z());
    }

    //new histogram PUResolution
  TH1F* PUResolution = dynamic_cast<TH1F*>(h["efficiency/PUResolution"]);
    if (!PUResolution) {
        std::cerr << "Error: Histogram PUResolution not found!" << std::endl;
        return;
    }
    for (size_t i = 0; i < simEvt.size(); ++i) {
    if (!simEvt[i].is_signal()) {  // Check if it is a pile-up event
        if (simEvt[i].is_matched()) {  // Check if the simulated vertex is matched, if not we cannot take the distance
            unsigned int rec_index = simEvt[i].rec;  // Get the index of the matched reconstructed vertex
            double true_z = simEvt[i].z;
            double rec_z = vtxs.at(rec_index).z();
            double delta_z = rec_z - true_z;

            // Fill the histogram with delta_z
           PUResolution->Fill(delta_z);
        }
    }
    }

        //definition of an 2H histogram


    
    // Initialize histograms for each category
    TH1F* PUConfusionMatrixCategorialC1 = dynamic_cast<TH1F*>(h["efficiency/reco_vs_true_z_position_hist_categorial_c1"]);
    if (!PUConfusionMatrixCategorialC1) {
        std::cerr << "Error: Histogram reco_vs_true_z_position_hist_categorial_c1 not found!" << std::endl;
        return;
    }

    TH1F* PUConfusionMatrixCategorialC2 = dynamic_cast<TH1F*>(h["efficiency/reco_vs_true_z_position_hist_categorial_c2"]);
    if (!PUConfusionMatrixCategorialC2) {
        std::cerr << "Error: Histogram reco_vs_true_z_position_hist_categorial_c2 not found!" << std::endl;
        return;
    }

    TH1F* PUConfusionMatrixCategorialC3 = dynamic_cast<TH1F*>(h["efficiency/reco_vs_true_z_position_hist_categorial_c3"]);
    if (!PUConfusionMatrixCategorialC3) {
        std::cerr << "Error: Histogram reco_vs_true_z_position_hist_categorial_c3 not found!" << std::endl;
        return;
        }
    // Map to track the count of reconstructed vertices for each simulated vertex
    std::unordered_map<unsigned int, int> simVertexToRecoCount;

    // Loop over all reconstructed vertices to count how many times each simulated vertex is matched
    for (unsigned int j = 0; j < vtxs.size(); ++j) {
        MVertex& vtx = vtxs.at(j);
        if (!vtx.isRecoFake() && vtx.sim != NOT_MATCHED_VTX_SIM) {
            ++simVertexToRecoCount[vtx.sim];// use a pre increment operator to increment value before returning it
           // cout << simVertexToRecoCount[vtx.sim] << endl; 
           /*
           be aware that this logic doesent quite work this way we have to look at the wos map (look for split_from
           in .H file or check teams post from saxermi1 on 9 of oktober)
           
           */
        }
    }

  // Loop over simulated vertices to categorize and fill histograms
for (size_t i = 0; i < simEvt.size(); i++) {
    // Check if the simulated event is matched to any reconstructed vertex
    if (simEvt[i].is_matched()) {
        // Attempt to find the count of reconstructed vertices associated with simulated vertex 'i'
        auto it = simVertexToRecoCount.find(i);
        if (it != simVertexToRecoCount.end()) {
            // Retrieved the count successfully
            int count = it->second;

            // Check if there is exactly one reconstructed vertex associated (Category 1)
            if (count == 1) {
                // Get the index of the associated reconstructed vertex
                unsigned int rec_index = simEvt[i].rec;
                // Retrieve the z-position of the reconstructed vertex
                double rec_z = vtxs.at(rec_index).z();
                PUConfusionMatrixCategorialC1->Fill(rec_z);
            } else if (simVertexToRecoCount[i] > 1) {
                // Category 2: Multiple reconstructed vertices per simulated vertex
                for (size_t j = 0; j < vtxs.size(); ++j) {
                    // Check if reconstructed vertex 'j' is associated with simulated vertex 'i'
                    if (vtxs.at(j).sim == i) {
                        // Retrieve the z-position of the reconstructed vertex
                        double rec_z = vtxs.at(j).z();
                        PUConfusionMatrixCategorialC2->Fill(rec_z);
                    }
                }
            }
            // Note: No else case; we use separate if statements for clarity and safety
        }
        // Optional: Handle the case where the simulated vertex 'i' is not found in the map
        // This means there are no reconstructed vertices associated with it
        else {
            // For example, you might want to log this event or perform additional processing
            // std::cout << "Simulated vertex " << i << " is matched but has no associated reconstructed vertices." << std::endl;
        }
    }
    // Optional: Handle the case where the simulated event is not matched
    // else {
    //     // For example, log or count unmatched simulated events
    // }
}


    // Category 3: Fake vertices (reconstructed vertices not matched to any simulated vertex)
    for (size_t i = 0; i < vtxs.size(); ++i) {
        MVertex& vtx = vtxs.at(i);
        if (vtx.isRecoFake() || vtx.sim == NOT_MATCHED_VTX_SIM) {
            double rec_z = vtx.z();
            PUConfusionMatrixCategorialC3->Fill(rec_z);
        }
    }
      // new histogram
      // this histogramm shows distance betweeen point in 3 d space and plane
      // we can later adapt this histogram to for example only show the distance between the border of the subspace and the fake vertices etc


       TH1F* TrueE3DDistanceToPlane = dynamic_cast<TH1F*>(h["efficiency/True_3D_point_to_plane_distance"]);
    if (!TrueE3DDistanceToPlane) {
        std::cerr << "Error: Histogram True_3D_point_to_plane_distance not found!" << std::endl;
        return;
        }
 
    // now we loop over the simulated vertices but first define the plane
     double plane_point[3] = {0.0, 0.0, 0.0};   // Point on the plane

    // Define the normal vector of the plane (45-degree slicing normal)
    double normal_vector[3] = {1.0 / sqrt(2), 0.0, -1.0 / sqrt(2)};

    for (size_t i = 0; i < simEvt.size(); ++i){
        if (simEvt[i].is_matched()) {

          double point[3] = {simEvt[i].x,
                             simEvt[i].y,
                             simEvt[i].z};
          double distance = distance_point_to_plane(point, plane_point, normal_vector);
          TrueE3DDistanceToPlane->Fill(distance);
        }
    }
  // histogram for SE track purity

  TH1F *SETracksPurity = dynamic_cast<TH1F *>(h["efficiency/SETracksPurity"]);
  if (!SETracksPurity) {
        std::cerr << "Error: Histogram SETracksPurity not found!" << std::endl;
        return;
    }

  if (simEvt[0].is_matched()) {
      MVertex& matchedVtx = vtxs.at(simEvt[0].rec);
      unsigned int numMatchedTracks = 0;
      unsigned int numTracks = matchedVtx.tracks.size();

        // Loop through the reconstructed tracks in the matched vertex
        for (auto tv : matchedVtx.tracks) {
            // Check if the track is matched to a simulated event
            unsigned int tk_sim = tracks.simevent_index_from_key(tv->key());
            assert(tv->_matched == tk_sim); // Ensure the track has the right matching

            // Check if the track is correctly assigned to the signal vertex
            bool correctly_assigned = (tk_sim == matchedVtx.sim);
            if (correctly_assigned) {
                numMatchedTracks++;
            }
        }

        // Calculate the purity as the fraction of correctly matched tracks
        float purity = (numTracks > 0) ? static_cast<float>(numMatchedTracks) / numTracks : 0;
        purity = purity * 100;

        // Fill the histogram with the calculated purity
        SETracksPurity->Fill(purity);
    } else {
        std::cerr << "No matched reconstructed vertex found!" << std::endl;
    }

    // histogram for SE track efficiency

    TH1F *SETracksEfficiency = dynamic_cast<TH1F *>(h["efficiency/SETracksEfficiency"]);
    if (!SETracksEfficiency) {
        std::cerr << "Error: Histogram SETracksEfficiency not found!" << std::endl;
        return;
    }

    // Ensure we have a signal vertex and it is matched
    if (simEvt.size() > 0 && simEvt[0].is_matched()) {
        
        unsigned int numSimTracks = simEvt[0].rtk.size(); // Number of simulated tracks
        unsigned int numMatchedTracks = 0; // Count of matched simulated tracks

        // Loop through the simulated tracks for the signal vertex
        for (auto& simTrack : simEvt[0].rtk) {
            // Check if this simulated track has a matching reconstructed track
            if (simTrack.is_matched()) {
                numMatchedTracks++;
            }
        }

        // Calculate the efficiency as the fraction of simulated tracks that are matched
        float efficiency = (numSimTracks > 0) ? static_cast<float>(numMatchedTracks) / numSimTracks : 0;
        efficiency = efficiency * 100;
        // Fill the histogram with the calculated efficiency
        SETracksEfficiency->Fill(efficiency);

    } else {
        std::cerr << "No signal vertex or no matched reconstructed vertex found!" << std::endl;
    }


    // histogram for PU track purity

    TH1F *PUTracksPurity = dynamic_cast<TH1F *>(h["efficiency/PUTracksPurity"]);
    if (!PUTracksPurity) {
          std::cerr << "Error: Histogram PUTracksPurity not found!" << std::endl;
          return;
      }

    for (size_t i = 1; i < simEvt.size(); i++) {
      if (simEvt[i].is_matched()) {
        MVertex& matchedVtx = vtxs.at(simEvt[i].rec);
        unsigned int numMatchedTracks = 0;
        unsigned int numTracks = matchedVtx.tracks.size();

          // Loop through the reconstructed tracks in the matched vertex
          for (auto tv : matchedVtx.tracks) {
              // Check if the track is matched to a simulated event
              unsigned int tk_sim = tracks.simevent_index_from_key(tv->key());
              assert(tv->_matched == tk_sim); // Ensure the track has the right matching

              // Check if the track is correctly assigned to the signal vertex
              bool correctly_assigned = (tk_sim == matchedVtx.sim);
              if (correctly_assigned) {
                  numMatchedTracks++;
              }
          }

          // Calculate the purity as the fraction of correctly matched tracks
          float purity = (numTracks > 0) ? static_cast<float>(numMatchedTracks) / numTracks : 0;
          purity = purity * 100;

          // Fill the histogram with the calculated purity
          PUTracksPurity->Fill(purity);
      } else {
          std::cerr << "No matched reconstructed vertex found!" << std::endl;
      }
    }

    // histogram for PU track efficiency

    TH1F *PUTracksEfficiency = dynamic_cast<TH1F *>(h["efficiency/PUTracksEfficiency"]);
    if (!PUTracksEfficiency) {
        std::cerr << "Error: Histogram PUTracksEfficiency not found!" << std::endl;
        return;
    }

    for (size_t i = 1; i < simEvt.size(); i++) {
      if (simEvt[i].is_matched()){
        unsigned int numSimTracks = simEvt[i].rtk.size(); // Number of simulated tracks
        unsigned int numMatchedTracks = 0; // Count of matched simulated tracks
        for (auto& simTrack : simEvt[i].rtk) {
              // Check if this simulated track has a matching reconstructed track
              if (simTrack.is_matched()) {
                  numMatchedTracks++;
              }
          }

        // Calculate the efficiency as the fraction of simulated tracks that are matched
        float efficiency = (numSimTracks > 0) ? static_cast<float>(numMatchedTracks) / numSimTracks : 0;
        efficiency = efficiency * 100;

        // Fill the histogram with the calculated efficiency
        PUTracksEfficiency->Fill(efficiency);

      } else {
          std::cerr << "No simulated vertex or no matched reconstructed vertex found!" << std::endl;
      }
    }
    

}

/*********************************************************************************************/






/*********************************************************************************************/
void TestAnalyzer::analyzeVertexCollectionZmatching(std::map<std::string, TH1*>& h,
								MVertexCollection& vtxs,
								std::vector<SimEvent>& simEvts,
								const std::string message,
								const double zwindow_sigmas
								)
{
  if (verbose_) {
    cout << "analyzeVertexCollectionZmatching, simEvts= " << simEvts.size() << " vtxs= " << vtxs.size()<<  " window size = " <<zwindow_sigmas  <<  endl;
  }

  unsigned int nsim = simEvts.size();
  unsigned int nvtx = vtxs.size();
  
  std::vector<unsigned int> nsimmatch(nvtx);  // count sim vertices within <zwindow_sigmas> sigma of the reconstructed position
  std::vector<unsigned int> nrecmatch(nsim);  // count rec vertices with reconstructed position within  <zwindow_sigmas> sigma 
  // repeat with the additional requirement of 2 common truth matched tracks
  std::vector<unsigned int> nsimmatch_tp(nvtx);
  std::vector<unsigned int> nrecmatch_tp(nsim);
  // and once more with a modified z-matching
  std::vector<unsigned int> nsimmatch_c(nvtx);
  std::vector<unsigned int> nrecmatch_c(nsim);


  for(unsigned int k=0; k < nvtx; k++){

    const auto v = vtxs.at(k);
    if (v.isRecoFake()) continue;

    unsigned int nearest_sim = NOT_MATCHED_VTX_REC;

    for(unsigned int i=0; i < nsim; i++){

      if (std::abs(simEvts.at(i).z - v.z()) < (zwindow_sigmas * v.zError())){
	//std::cout << "zmatching " << message << "  [" << i << "]" << simEvts.at(i).z << "   " << v.z() << "+/-" << v.zError() << " (" << k << ")" << endl;
	nrecmatch[i]++;
	nsimmatch[k]++;
	//nmatchall++;
	if(nearest_sim == NOT_MATCHED_VTX_REC){
	  nearest_sim = i;
	}else{
	  if(std::abs(simEvts.at(i).z - v.z()) < std::abs(v.z()-simEvts.at(nearest_sim).z)){
	    nearest_sim = i;
	  }
	}
      }

      // same but requiring two matched tracks
      if ((std::abs(simEvts.at(i).z - v.z()) < (zwindow_sigmas * v.zError())) && (simEvts.at(i).countVertexTracks(v, 0.2) > 1)){
	nrecmatch_tp[i]++;
	nsimmatch_tp[k]++;
      }

      // modified z-matching
      double dzmax = std::min(0.1, std::max(0.0100, zwindow_sigmas * v.zError()));
      if (std::abs(simEvts.at(i).z - v.z()) < dzmax){
	nrecmatch_c[i]++;
	nsimmatch_c[k]++;
      }
      
    }// sim events
  }// rec vtxs



  // fill histograms
  // 
  for(unsigned int i=0; i < nsim; i++){
    Fill(h, "zmatcheffvspu", nsim, nrecmatch[i] > 0 ? 1. : 0.);
    Fill(h, "zmatchambigvspu", nsim, nrecmatch[i] > 1 ? 1. : 0.);
    Fill(h, "zmatchnrecmatch", nrecmatch[i]);

    Fill(h, "zcmatcheffvspu", nsim, nrecmatch_c[i] > 0 ? 1. : 0.);
    Fill(h, "zcmatchambigvspu", nsim, nrecmatch_c[i] > 1 ? 1. : 0.);
    Fill(h, "zcmatchnrecmatch", nrecmatch_c[i]);

    Fill(h, "ztpmatcheffvspu", nsim, nrecmatch_tp[i] > 0 ? 1. : 0.);
    Fill(h, "ztpmatchambigvspu", nsim, nrecmatch_tp[i] > 1 ? 1. : 0.);
    Fill(h, "ztpmatchnrecmatch", nrecmatch_tp[i]);
  }

  // fake here means : no sim vertex within <zwindow_sigmas> x sigma_z
  for(unsigned int k=0; k < nvtx; k++){
    Fill(h, "zmatchfakevspu", nsim, nsimmatch[k] == 0 ? 1. : 0);
    Fill(h, "zcmatchfakevspu", nsim, nsimmatch_c[k] == 0 ? 1. : 0);
    Fill(h, "ztpmatchfakevspu", nsim, nsimmatch_tp[k] == 0 ? 1. : 0);
  }

  for(unsigned int k=0; k < nvtx; k++){
    Fill(h, "zmatchnsimmatch", nsimmatch[k]);
    Fill(h, "zcmatchnsimmatch", nsimmatch_c[k]);
    Fill(h, "ztpmatchnsimmatch", nsimmatch_tp[k]);
  }

  


    
  // go one step beyond looking at z-windows :  maximally greedy matching
  std::vector<std::pair<unsigned int, unsigned int>> recsim;    // based on z-distance only
  std::vector<std::pair<unsigned int, unsigned int>> recsim_c;    // based on z-distance, truncated
  std::vector<std::pair<unsigned int, unsigned int>> recsim_tp;  // additionally require at least one truth matched track
  // and the same for selected vertices
  std::vector<std::pair<unsigned int, unsigned int>> recselsim;    // based on z-distance only
  std::vector<std::pair<unsigned int, unsigned int>> recselsim_c;    // based on z-distance, truncated
  std::vector<std::pair<unsigned int, unsigned int>> recselsim_tp;  // additionally require at least one truth matched track

  unsigned int nvtxrec = 0, nvtxsel = 0;
  for(unsigned int k = 0; k < nvtx; k++){
    const auto v = vtxs.at(k);
    if( v.isRecoFake() ) continue;
    
    nvtxrec ++;
    if(select(v)) nvtxsel++;


    for(unsigned int i = 0; i < nsim; i++){

      // z-distance alone
      if (std::abs(simEvts.at(i).z - v.z()) < (zwindow_sigmas * v.zError())){
	recsim.emplace_back(k,i);
	if(select(v)) recselsim.emplace_back(k,i);
      }

      // truncated z-distance, allow at least 100 um, do not allow more than 1 mm
      double dzmax = std::min(0.1, std::max(0.0100, zwindow_sigmas * v.zError()));
      if (std::abs(simEvts.at(i).z - v.z()) < dzmax){
	recsim_c.emplace_back(k,i);
	if(select(v)) recselsim_c.emplace_back(k,i);
      }

      // z-distance + at least two tracks with weight > 0.5
      if ((std::abs(simEvts.at(i).z - v.z()) < (zwindow_sigmas * v.zError())) && (simEvts.at(i).countVertexTracks(v, 0.5) > 1)){
	recsim_tp.emplace_back(k,i);
	if(select(v)) recselsim_tp.emplace_back(k,i);
      }
    }
  }

  if( nsim > 0){
    int max_match = 0, max_match_c = 0,max_match_tp = 0;
    FFA(nvtxrec, nsim, recsim, max_match);
    Fill(h, "FFAzmatcheffvspu", nsim, float(max_match) / nsim);

    FFA(nvtxrec, nsim, recsim_c, max_match_c);
    Fill(h, "FFAzcmatcheffvspu", nsim, float(max_match_c) / nsim);

    FFA(nvtxrec, nsim, recsim_tp, max_match_tp);
    Fill(h, "FFAztpmatcheffvspu", nsim, float(max_match_tp) / nsim);

    if(nvtxrec > 0){
      Fill(h, "FFAzmatchfakevspu", nsim, float(nvtxrec - max_match) / nvtxrec);
      Fill(h, "FFAzcmatchfakevspu", nsim, float(nvtxrec - max_match_c) / nvtxrec);
      Fill(h, "FFAztpmatchfakevspu", nsim, float(nvtxrec - max_match_tp) / nvtxrec);
      //should the denominator contain only vertices with tp matchted tracks?
    }

    // include selection
    max_match = 0; 
    max_match_c = 0;
    max_match_tp = 0;
    FFA(nvtxsel, nsim, recselsim, max_match);
    Fill(h, "FFAzmatchseleffvspu", nsim, float(max_match) / nsim);

    FFA(nvtxsel, nsim, recselsim_c, max_match_c);
    Fill(h, "FFAzcmatchseleffvspu", nsim, float(max_match_c) / nsim);

    FFA(nvtxsel, nsim, recselsim_tp, max_match_tp);
    Fill(h, "FFAztpmatchseleffvspu", nsim, float(max_match_tp) / nsim);

    if(nvtxsel > 0){
      Fill(h, "FFAzmatchselfakevspu", nsim, float(nvtxsel - max_match) / nvtxsel);
      Fill(h, "FFAzcmatchselfakevspu", nsim, float(nvtxsel - max_match_c) / nvtxsel);
      Fill(h, "FFAztpmatchselfakevspu", nsim, float(nvtxsel - max_match_tp) / nvtxsel);
    }

  }
  
}
/*****************************analyzeVertexCollectionZmatching**********************************/







/***************************************************************************************/
void TestAnalyzer::analyzeVertexRecoCPUTime(std::map<std::string, TH1*>& h,
                                                        const reco::VertexCollection* recVtxs,
                                                        const std::string message)
/***************************************************************************************/
{
  double tclu = 0;
  double tfit = 0;
  double ttime = 0;
  int nsel = 0;
  bool found_timing_info = false;


  for (reco::VertexCollection::const_iterator v = recVtxs->begin(); v != recVtxs->end(); ++v) {
    if (v->isFake()) {
      tclu = v->covariance(iX, iX);
      tfit = v->covariance(iY, iY);
      ttime = v->covariance(iZ, iZ);
      found_timing_info = true;
      //std::cout << "extracted timing info " << message << "  tclu = "  << tclu << " ms   tfit =" << tfit << " ms" << endl;
    } else {
      if (select(*v)) {
        nsel++;
      }
    }
  }

  if (found_timing_info) {
    Fill(h, "cputime/tcluvsSimPU", simPU_, tclu);
    Fill(h, "cputime/tfitvsSimPU", simPU_, tfit);
    Fill(h, "cputime/ttimevsSimPU", simPU_, ttime);
    Fill(h, "cputime/tcluvsLPU", lumiPU_, tclu);
    Fill(h, "cputime/tfitvsLPU", lumiPU_, tfit);
    Fill(h, "cputime/ttimevsLPU", lumiPU_, ttime);
    Fill(h, "cputime/tcluvsnsel", nsel, tclu);
    Fill(h, "cputime/tfitvsnsel", nsel, tfit);
    Fill(h, "cputime/ttimevsnsel", nsel, ttime);
  }
}
/***************************************************************************************/





 

//define this as a plug-in
DEFINE_FWK_MODULE(TestAnalyzer);

