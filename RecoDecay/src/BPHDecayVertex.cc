/*
 *  See header file for a description of this class.
 *
 *  $Date: 2015-07-03 13:49:53 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHDecayVertex.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>

using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
BPHDecayVertex::BPHDecayVertex( const edm::EventSetup* es ):
 evSetup( es ) {
  setNotUpdated();
}


BPHDecayVertex::BPHDecayVertex( const BPHDecayVertex* ptr,
                                const edm::EventSetup* es ):
 evSetup( es ) {
  const vector<Component>& list = ptr->BPHDecayMomentum::componentList();
  int i;
  int n = list.size();
  for ( i = 0; i < n; ++i ) {
    const Component& component = list[i];
    searchMap[component.cand] = component.searchList;
  }
  const vector<BPHRecoConstCandPtr>& dComp = daughComp();
  n = dComp.size();
  while ( n-- ) {
    const map<const reco::Candidate*,string>& dMap = dComp[n]->searchMap;
    searchMap.insert( dMap.begin(), dMap.end() );
  }
  setNotUpdated();
}

//--------------
// Destructor --
//--------------
void BPHDecayVertex::add( const string& name,
                          const reco::Candidate* daug, 
                          const string& searchList,
                          double mass ) {
  BPHDecayMomentum::add( name, daug, mass );
  searchMap[daug] = searchList;
  return;
}


BPHDecayVertex::~BPHDecayVertex() {
}

//--------------
// Operations --
//--------------
bool BPHDecayVertex::isValidVertex() const {
  if ( oldVertex ) fitVertex();
  return validVertex;
}


const reco::Vertex& BPHDecayVertex::vertex() const {
  if ( oldVertex ) fitVertex();
  return fittedVertex;
}


const vector<const reco::Track*>& BPHDecayVertex::tracks() const {
  if ( oldTracks ) tTracks();
  return rTracks;
}


const reco::Track* BPHDecayVertex::getTrack(
                                   const reco::Candidate* cand ) const {
  if ( oldTracks ) tTracks();
  map<const reco::Candidate*,
      const reco::Track*>::const_iterator iter = tkMap.find( cand );
  map<const reco::Candidate*,
      const reco::Track*>::const_iterator iend = tkMap.end();
  if ( iter == iend ) iter = tkMap.find( originalReco( cand ) );
  return ( iter != iend ? iter->second : 0 );
}


const vector<reco::TransientTrack>& BPHDecayVertex::transientTracks() const {
  if ( oldTracks ) tTracks();
  return trTracks;
}


reco::TransientTrack* BPHDecayVertex::getTransientTrack(
                                      const reco::Candidate* cand ) const {
  if ( oldTracks ) tTracks();
  map<const reco::Candidate*,
            reco::TransientTrack*>::const_iterator iter = ttMap.find( cand );
  map<const reco::Candidate*,
            reco::TransientTrack*>::const_iterator iend = ttMap.end();
  if ( iter == iend ) iter = ttMap.find( originalReco( cand ) );
  return ( iter != iend ? iter->second : 0 );
}


void BPHDecayVertex::setNotUpdated() const {
  BPHDecayMomentum::setNotUpdated();
  oldTracks = oldVertex = true;
  validVertex = false;
  return;
}


void BPHDecayVertex::tTracks() const {
  oldTracks = false;
   rTracks.clear();
  trTracks.clear();
  tkMap.clear();
  ttMap.clear();
  edm::ESHandle<TransientTrackBuilder> ttB;
  evSetup->get<TransientTrackRecord>().get( "TransientTrackBuilder", ttB );
  const vector<const reco::Candidate*>& dL = daughFull();
  int n = dL.size();
  trTracks.reserve( n );
  validVertex = true;
  while ( n-- ) {
    const reco::Candidate* rp = originalReco( dL[n] );
    tkMap[rp] = 0;
    ttMap[rp] = 0;
    if ( !rp->charge() ) continue;
    const reco::Track* tp;
    const char* searchList = "cfhp";
    map<const reco::Candidate*,string>::const_iterator iter =
                                                       searchMap.find( rp );
    if ( iter != searchMap.end() ) searchList = iter->second.c_str();
    tp = BPHTrackReference::getTrack( *rp, searchList );
    if ( tp == 0 ) {
      edm::LogPrint( "DataNotFound" )
                  << "BPHDecayVertex::tTracks: "
                  << "no track for reco::(PF)Candidate";
      validVertex = false;
      continue;
    }
     rTracks.push_back( tp );
    trTracks.push_back( ttB->build( tp ) );
    reco::TransientTrack* ttp = &trTracks.back();
    tkMap[rp] =  tp;
    ttMap[rp] = ttp;
  }
  return;
}


void BPHDecayVertex::fitVertex() const {
  oldVertex = false;
  if ( oldTracks ) tTracks();
  if ( trTracks.size() < 2 ) return;
  KalmanVertexFitter kvf( true );
  TransientVertex tv = kvf.vertex( trTracks );
  fittedVertex = tv;
  return;
}

