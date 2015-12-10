#ifndef BPHTrackReference_H
#define BPHTrackReference_H
/** \class BPHTrackReference
 *
 *  Description: 
 *
 *
 *  $Date: 2015-11-19 10:59:18 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

//---------------
// C++ Headers --
//---------------
#include <iostream>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHTrackReference {

 public:

  typedef reco::PFCandidate candidate;

  /** Constructor
   */
  BPHTrackReference();

  /** Destructor
   */
  ~BPHTrackReference();

  /** Operations
   */
  /// 
  static const reco::Track* getTrack( const reco::Candidate& rc,
                            const char* modeList = "cbfpmige" ) {
    const char* mptr = modeList;
    char mode;
    const reco::Track* tkp = 0;
    while ( ( mode = *mptr++ ) ) {
      switch ( mode ) {
      case 'c': if ( ( tkp = getFromRC( rc ) ) != 0 ) return tkp; break;
      case 'b': if ( ( tkp = getFromBT( rc ) ) != 0 ) return tkp; break;
      case 'f': if ( ( tkp = getFromPF( rc ) ) != 0 ) return tkp; break;
      case 'p': if ( ( tkp = getFromPC( rc ) ) != 0 ) return tkp; break;
      case 'm': if ( ( tkp = getMuonPF( rc ) ) != 0 ) return tkp; break;
      case 'i': if ( ( tkp = getMuonIT( rc ) ) != 0 ) return tkp; break;
      case 'g': if ( ( tkp = getMuonGT( rc ) ) != 0 ) return tkp; break;
      case 'e': if ( ( tkp = getElecPF( rc ) ) != 0 ) return tkp; break;
      }
    }
    return 0;
  }

  static const reco::Track* getFromRC( const reco::Candidate& rc ) {
//    std::cout << "getFromRC" << std::endl;
    try {
      const reco::TrackRef& tkr = rc.get<reco::TrackRef>();
      if ( !tkr.isNull() ) return tkr.get();
    }
    catch ( edm::Exception e ) {
    }
    return 0;
  }
  static const reco::Track* getFromBT( const reco::Candidate& rc ) {
//    std::cout << "getFromBT" << std::endl;
//    try {
//      const reco::Track* trk = rc.bestTrack();
//      return trk;
//    }
//    catch ( edm::Exception e ) {
//    }
    return 0;
  }
  static const reco::Track* getFromPF( const reco::Candidate& rc ) {
//    std::cout << "getFromPF" << std::endl;
    const reco::PFCandidate* pf =
          dynamic_cast<const reco::PFCandidate*>( &rc );
    if ( pf == 0 ) return 0;
    try {
      const reco::TrackRef& tkr = pf->trackRef();
      if ( !tkr.isNull() ) return tkr.get();
    }
    catch ( edm::Exception e ) {
    }
    return 0;
  }
  static const reco::Track* getFromPC( const reco::Candidate& rc ) {
//    std::cout << "getFromPC" << std::endl;
//    const pat::PackedCandidate* pp =
//        dynamic_cast<const pat::PackedCandidate*>( &rc );
//    if ( pp != 0 ) return &pp->pseudoTrack();
    return 0;
  }
  static const reco::Track* getMuonPF( const reco::Candidate& rc ) {
//    std::cout << "getMuonPF" << std::endl;
    const pat::Muon* mu = dynamic_cast<const pat::Muon*>( &rc );
    if ( mu == 0 ) return 0;
    try {
      const reco::PFCandidateRef& pcr = mu->pfCandidateRef();
      if ( !pcr.isNull() ) {
        const reco::TrackRef& tkr = pcr->trackRef();
        if ( !tkr.isNull() ) return tkr.get();
      }
    }
    catch ( edm::Exception e ) {
    }
    return 0;
  }
  static const reco::Track* getMuonIT( const reco::Candidate& rc ) {
//    std::cout << "getMuonIT" << std::endl;
    const pat::Muon* mu = dynamic_cast<const pat::Muon*>( &rc );
    if ( mu == 0 ) return 0;
    if ( !mu->isTrackerMuon() ) return 0;
    try {
      const reco::TrackRef& mit = mu->innerTrack();
      if ( !mit.isNull() ) return mit.get();
    }
    catch ( edm::Exception e ) {
    }
    return 0;
  }
  static const reco::Track* getMuonGT( const reco::Candidate& rc ) {
//    std::cout << "getMuonGT" << std::endl;
    const pat::Muon* mu = dynamic_cast<const pat::Muon*>( &rc );
    if ( mu == 0 ) return 0;
    if ( !mu->isGlobalMuon() ) return 0;
    try {
      const reco::TrackRef& mgt = mu->globalTrack();
      if ( !mgt.isNull() ) return mgt.get();
    }
    catch ( edm::Exception e ) {
    }
    return 0;
  }
  static const reco::Track* getElecPF( const reco::Candidate& rc ) {
//    std::cout << "getElecPF" << std::endl;
    const pat::Electron* el = dynamic_cast<const pat::Electron*>( &rc );
    if ( el == 0 ) return 0;
    try {
      const reco::PFCandidateRef& pcr = el->pfCandidateRef();
      if ( !pcr.isNull() ) {
        const reco::TrackRef& tkr = pcr->trackRef();
        if ( !tkr.isNull() ) return tkr.get();
      }
    }
    catch ( edm::Exception e ) {
    }
    return 0;
  }

 private:

  // private copy and assigment constructors
  BPHTrackReference           ( const BPHTrackReference& x );
  BPHTrackReference& operator=( const BPHTrackReference& x );

};


#endif // BPHTrackReference_H
