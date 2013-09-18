// -*- C++ -*-
//
// Package:   hitEffAnalyzer
// Class:      hitEffAnalyzer
// 
/**\class hitEffAnalyzer hitEffAnalyzer.cc hitEff/hitEffAnalyzer/src/hitEffAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Jessica Danielle Brinson,42 2-032,+41227662377,
//         Created:  Tue Aug 27 11:00:02 CEST 2013
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <math.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <map>
#include <set>
#include <vector>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"  

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/SiStripDetId/interface/TIBDetId.h"  
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"  
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"  
#include "DataFormats/SiStripDetId/interface/TECDetId.h"  

#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include <TrackingTools/MeasurementDet/interface/MeasurementDet.h>

#include "RecoTracker/TrackProducer/plugins/TrackProducer.h"

using namespace reco;
using namespace std;
using namespace edm;

//
// class declaration
//
class TrackAssociatorBase;
class TrackAssociatorByHits; 

class hitEffAnalyzer : public edm::EDAnalyzer {
  
public:
  hitEffAnalyzer(const edm::ParameterSet& conf);
  virtual ~hitEffAnalyzer();

  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  int getLayerHit(unsigned int rawId); 
  int getSubDet(unsigned int rawId);
  double getLocalPosXNorm(PSimHit& simHit,const MeasurementDet* measDet);
  double getLocalPosYNorm(PSimHit& simHit,const MeasurementDet* measDet);
private:
  TrackAssociatorBase * associatorByHits;

  // ----------member data ---------------------------
  edm::InputTag tracksTag, tpTag, simtracksTag;
  edm::Service<TFileService> fs;

 
  struct RecHitInfo {
    // info about a rec hit
    uint detRawId;  
    double energyLoss;
    double locPosXNorm;
    double locPosYNorm;
    int subdet;  // subdetector, TIB=3, TID=4, TOB=5, TEC=6, to match definition in http://cmslxr.fnal.gov/lxr/source/DataFormats/SiStripDetId/interface/SiStripDetId.h
    int layer;  // or wheel, or disk
    double eta;
    double phi;  
    GlobalPoint globalPt;  
    
    RecHitInfo() {  // constructor
      detRawId = -99;  
      energyLoss = -99;
      locPosXNorm = -99;
      locPosYNorm = -99;
      subdet = -99;  
      layer = -99;
      eta = -99;
      phi = -99;
    }
    
  };

  struct TrackAndHits {
    double pt;
    double eta;
    double phi;
    vector<RecHitInfo> recHits;
  };  

  //
  //declare histograms
  //
  TH1D* hNTrks;
  TH1D* hDistXMatchedRecHit;
  TH1D* hDistYMatchedRecHit;
  TH1D* hDistXNoMatchedRecHit;
  TH1D* hDistYNoMatchedRecHit;

};

//
// constructors and destructor
//

hitEffAnalyzer::hitEffAnalyzer(edm::ParameterSet const& conf) 
{
  //now do what ever initialization is needed

  tracksTag    = conf.getParameter< edm::InputTag >("tracksTag");
  tpTag        = conf.getParameter< edm::InputTag >("tpTag");
  simtracksTag = conf.getParameter< edm::InputTag >("simtracksTag");
  //
  //book histograms
  //
  hNTrks                = fs->make<TH1D>("numTrks"              ,";Num. trks"                     ,11 , -0.5, 10.5);
  hDistXMatchedRecHit   = fs->make<TH1D>("hDistXMatchedRecHit"  ,";x distance to sensor edge (cm)",100, -1.5, 1.5); 
  hDistYMatchedRecHit   = fs->make<TH1D>("hDistYMatchedRecHit"  ,";y distance to sensor edge (cm)",100, -1.5, 1.5); 
  hDistXNoMatchedRecHit = fs->make<TH1D>("hDistXNoMatchedRecHit",";x distance to sensor edge (cm)",100, -1.5, 1.5); 
  hDistYNoMatchedRecHit = fs->make<TH1D>("hDistYNoMatchedRecHit",";y distance to sensor edge (cm)",100, -1.5, 1.5); 

}
  
hitEffAnalyzer::~hitEffAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}
  
// ------------ method called for each event  ------------
void hitEffAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& setup)
{

  Handle<View<Track> > trackCollectionH;
  event.getByLabel(tracksTag,trackCollectionH);
  const View<Track>  tC = *(trackCollectionH.product()); 
 
  Handle<SimTrackContainer> simTrackCollection;
  event.getByLabel(simtracksTag, simTrackCollection);
  const SimTrackContainer simTC = *(simTrackCollection.product());
 
  Handle<TrackingParticleCollection>  TPCollectionH ;
  event.getByLabel(tpTag,TPCollectionH);
  const TrackingParticleCollection tPC   = *(TPCollectionH.product());

  Handle<TrajTrackAssociationCollection> trajTrackAssociationHandle;
  event.getByLabel("TrackRefitter", trajTrackAssociationHandle);  
  const TrajTrackAssociationCollection& TrajToTrackMap = *trajTrackAssociationHandle.product();
 
  Handle<TrackCollection> tracks;
  event.getByLabel(tracksTag,tracks);

  ESHandle<TrackAssociatorBase> theHitsAssociator;  
  setup.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",theHitsAssociator);
  associatorByHits = (TrackAssociatorBase *) theHitsAssociator.product();

  ESHandle<MeasurementTracker> theMeasTk;
  setup.get<CkfComponentsRecord>().get("", theMeasTk); 


  // ********************************************************
  // get info about recHits and store in struct
  // ********************************************************
  vector<TrackAndHits> trkRecHits;
  vector<RecHitInfo> recHits;
  
  double minTrkPt = 5;  
  int j=0;
  TrackCollection::const_iterator itTrack = tracks->begin();
  for(TrajTrackAssociationCollection::const_iterator cit=TrajToTrackMap.begin(); 
      cit!=TrajToTrackMap.end() && itTrack!=tracks->end(); 
      cit++, j++, itTrack++){
    const reco::TrackRef track = cit->val;  
    const edm::Ref<std::vector<Trajectory> > traj = cit->key;
    
    if (itTrack->pt() < minTrkPt) continue;  
    
    cout << "***********************************************************" << endl;  
    
    const vector<TrajectoryMeasurement> & measurements = traj->measurements();
    
    TrackAndHits newTrk;
    newTrk.pt  = itTrack->pt();
    newTrk.eta = itTrack->eta();
    newTrk.phi = itTrack->phi();
    
    // Loop over rechits
    for(vector<TrajectoryMeasurement>::const_iterator it = measurements.begin(); it!=measurements.end(); it++){
      TrajectoryStateOnSurface trajState=it->updatedState();
      if( !trajState.isValid()) continue;      
      const TrackingRecHit * recHit=(*it->recHit()).hit();
      if (!recHit) continue;

      RecHitInfo newhit;
      newhit.detRawId = recHit->geographicalId().rawId();   
      newhit.globalPt = trajState.globalPosition();  
      newhit.layer = getLayerHit(newhit.detRawId);  
      newhit.subdet = getSubDet(newhit.detRawId);  
      newhit.eta = newTrk.eta;
      newhit.phi = newTrk.phi;
      newTrk.recHits.push_back(newhit);        
    }
    trkRecHits.push_back(newTrk); 
  }

  // ********************************************************
  //associate reco'd tracks to sim tracks
  // ********************************************************
  cout << "                      ****************** Reco To Sim ****************** " << endl;
  cout << "-- Associator by hits --" << endl;  
  reco::RecoToSimCollection p = 
    associatorByHits->associateRecoToSim (trackCollectionH,TPCollectionH,&event );

  if (tC.size() != 2) { 
    cout << "Warning: found number of tracks in collection:  " << tC.size()
	 << "; skipping this event." << endl;
    return;  // only consider events with exactly two reco'd tracks  
  }
  
  hNTrks->Fill(tC.size());  

  for(View<Track>::size_type i=0; i<tC.size(); ++i) {
    RefToBase<Track> track(trackCollectionH, i);

    try{ 
      std::vector<std::pair<TrackingParticleRef, double> > tp = p[track];
      cout << "Reco Track pT: "  << setw(6) << track->pt() 
	   <<  " matched to " << tp.size() << " MC Tracks" << std::endl;

      for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin(); 
	   it != tp.end(); ++it) {
	TrackingParticleRef tpr = it->first;
	vector<PSimHit> simHits = tpr->trackPSimHit();  

	
	int nHitPix = 0;
	for (uint ihit=0; ihit<simHits.size(); ihit++) {  
	  DetId detid = DetId(simHits.at(ihit).detUnitId());  
	  if (detid.det() != DetId::Tracker) continue;  
	  if (detid.subdetId()<3) { 		  
	    nHitPix++;  
	    continue;  // only consider strip hits
	  }
	  if (detid.subdetId()>6) continue; 

	  bool hasMatchedRecHit = false;  	    
	  for(trackingRecHit_iterator iter = track->recHitsBegin(); 
	      iter != track->recHitsEnd(); iter++) {
	    if (!(*iter)->isValid()) continue;	      
	    if (simHits.at(ihit).detUnitId() == (*iter)->rawId()) hasMatchedRecHit = true;  
	  }  

	  //
	  //Fill hists for simHits matched/unmatched to recHit
	  //
	  if (hasMatchedRecHit) {
	    hDistXMatchedRecHit->Fill( getLocalPosXNorm(simHits.at(ihit), theMeasTk->idToDet(detid)) );  
	    hDistYMatchedRecHit->Fill( getLocalPosYNorm(simHits.at(ihit), theMeasTk->idToDet(detid)) );  
	  } else{ 
	    hDistXNoMatchedRecHit->Fill( getLocalPosXNorm(simHits.at(ihit), theMeasTk->idToDet(detid)) );  
	    hDistYNoMatchedRecHit->Fill( getLocalPosYNorm(simHits.at(ihit), theMeasTk->idToDet(detid)) );  
	  } 

	} //end for (uint ihit=0; ihit<simHits.size(); ihit++) {	  
      } //end  for (std::vector<std::pair<TrackingParticleRef, double> >::const_iterator it = tp.begin();
    }

    catch (Exception event) {
      cout << "->   Track pT: " 
	   << setprecision(2) << setw(6) << track->pt() 
	   <<  " matched to 0  MC Tracks" << endl;
    }    
 
  } //end for(View<Track>::size_type i=0; i<tC.size(); ++i) {

}
// ------------ method called once each job just before starting event loop  ------------
void 
hitEffAnalyzer::beginJob() 
{
}

// ------------ member functions  ------------
int hitEffAnalyzer::getSubDet(unsigned int rawId) {
  // return the subdetector ID
  DetId detid = DetId(rawId);
  return detid.subdetId();  
}


int hitEffAnalyzer::getLayerHit(unsigned int rawId) {
  //rawId must be the raw detector id

  SiStripDetId strDetId = SiStripDetId(rawId);  

  int layer = -1;  
  if (strDetId.subDetector() == SiStripDetId::TIB) {  
    TIBDetId tibid = TIBDetId(rawId);  
    layer = tibid.layerNumber();  
  } else if (strDetId.subDetector() == SiStripDetId::TOB) {  
    TOBDetId tobid = TOBDetId(rawId);  
    layer = tobid.layerNumber();  
  } else if (strDetId.subDetector() == SiStripDetId::TID) {  
    TIDDetId tidid = TIDDetId(rawId);  
    layer = tidid.diskNumber();  
  } else if (strDetId.subDetector() == SiStripDetId::TEC) {  
    TECDetId tecid = TECDetId(rawId);  
    layer = tecid.wheelNumber();  
  } 

  return layer;  

}  

double hitEffAnalyzer::getLocalPosXNorm(PSimHit& simHit,const MeasurementDet* measDet) {
  LocalPoint HitLocalPos = simHit.localPosition();  
  double HalfWidth      = measDet->surface().bounds().width()  /2.0;  // in X direction  
  double HalfLength     = measDet->surface().bounds().length() /2.0;  // in Y direction        
  const TrapezoidalPlaneBounds* trapezoidalBounds( dynamic_cast<const TrapezoidalPlaneBounds*>(&(measDet->surface().bounds())));   
  const RectangularPlaneBounds* rectangularBounds( dynamic_cast<const RectangularPlaneBounds*>(&(measDet->surface().bounds())));   
  if (trapezoidalBounds) {  
    std::vector<float> const & parameters = (*trapezoidalBounds).parameters();
    HalfLength     = parameters[3];  
    double t       = (HalfLength + HitLocalPos.y()) / (2*HalfLength); 
    HalfWidth      = parameters[0] + (parameters[1]-parameters[0]) * t;  
    //      cout << "Warning:  using trapezoidal bounds." << endl;  
  } else if (!rectangularBounds) {
    cout << "Warning:  neither trapezoidal nor rectangular bounds found!" << endl;
    return -99.; 
  }       
  double locPosXNorm =  HitLocalPos.x() / HalfWidth;  
  return locPosXNorm; 	      

}

double hitEffAnalyzer::getLocalPosYNorm(PSimHit& simHit,const MeasurementDet* measDet) {
  LocalPoint HitLocalPos = simHit.localPosition();  
  double HalfWidth      = measDet->surface().bounds().width()  /2.0;  // in X direction  
  double HalfLength     = measDet->surface().bounds().length() /2.0;  // in Y direction        
  const TrapezoidalPlaneBounds* trapezoidalBounds( dynamic_cast<const TrapezoidalPlaneBounds*>(&(measDet->surface().bounds())));   
  const RectangularPlaneBounds* rectangularBounds( dynamic_cast<const RectangularPlaneBounds*>(&(measDet->surface().bounds())));   
  if (trapezoidalBounds) {  
    std::vector<float> const & parameters = (*trapezoidalBounds).parameters();
    HalfLength     = parameters[3];  
    double t       = (HalfLength + HitLocalPos.y()) / (2*HalfLength); 
    HalfWidth      = parameters[0] + (parameters[1]-parameters[0]) * t;  
    //      cout << "Warning:  using trapezoidal bounds." << endl;  
  } else if (!rectangularBounds) {
    cout << "Warning:  neither trapezoidal nor rectangular bounds found!" << endl;
    return -99.; 
  }       
  double locPosYNorm =  HitLocalPos.y() / HalfLength;  
  return locPosYNorm; 	      

}

// ------------ method called once each job just after ending the event loop  ------------
void 
hitEffAnalyzer::endJob() 
{
}

DEFINE_FWK_MODULE(hitEffAnalyzer);
