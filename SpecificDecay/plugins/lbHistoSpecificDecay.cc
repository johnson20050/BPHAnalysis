#include "BPHAnalysis/SpecificDecay/plugins/lbHistoSpecificDecay.h"

#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"

#include <TH1.h>
#include <TFile.h>

#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define SET_LABEL(NAME,PSET) ( NAME = PSET.getParameter<string>( #NAME ) )
// SET_LABEL(xyz,ps);
// is equivalent to
// xyz = ps.getParameter<string>( "xyx" )

// useful classes {{{
class BPHUserData {
 public:
  template<class T>
  static const T* get( const pat::CompositeCandidate& cand,
                       const string& name ) {
    if ( cand.hasUserData( name ) ) return cand.userData<T>( name );
    return 0;
  }
  template<class T>
  static const T* getByRef( const pat::CompositeCandidate& cand,
                            const string& name ) {
    if ( cand.hasUserData( name ) ) {
      typedef edm::Ref< std::vector<T> > objRef;
      const objRef* ref = cand.userData<objRef>( name );
      if ( ref ==      0 ) return 0;
      if ( ref->isNull() ) return 0;
      return ref->get();
    }
    return 0;
  }
};


class BPHDaughters {
 public:
  static vector<const reco::Candidate*> get(
         const pat::CompositeCandidate& cand, float massMin, float massMax ) {
    int i;
    int n = cand.numberOfDaughters();
    vector<const reco::Candidate*> v;
    v.reserve( n );
    for ( i = 0; i < n; ++i ) {
      const reco::Candidate* dptr = cand.daughter( i );
      float mass = dptr->mass();
      if ( ( mass > massMin ) && ( mass < massMax ) ) v.push_back( dptr );
    }
    return v;
  }
};


class BPHSoftMuonSelect {
 public:
  BPHSoftMuonSelect( int   cutTrackerLayers = 5,
                     int   cutPixelLayers   = 0,
                     float maxDxy           = 0.3,
                     float maxDz            = 20.0,
                     bool  goodMuon         = true ,
                     bool  highPurity       = true ):
   cutTL ( cutTrackerLayers ),
   cutPL ( cutPixelLayers   ),
   maxXY ( maxDxy           ),
   maxZ  ( maxDz            ),
   gM    ( goodMuon         ),
   hP    ( highPurity       ) {}
  bool accept( const reco::Candidate& cand,
               const reco::Vertex* pv ) const {
    const pat::Muon* p = dynamic_cast<const pat::Muon*>( &cand );
    if ( p == 0 ) return false;
    if ( gM &&
        !muon::isGoodMuon( *p, muon::TMOneStationTight ) ) return false;
    if ( p->innerTrack()->hitPattern().trackerLayersWithMeasurement() <= cutTL )
         return false;
    if ( p->innerTrack()->hitPattern().  pixelLayersWithMeasurement() <= cutPL )
         return false;
    if ( hP &&
        !p->innerTrack()->quality( reco::TrackBase::highPurity ) )
         return false;
    if ( pv == 0 ) return true;
    const reco::Vertex::Point& pos = pv->position();
    if ( fabs( p->innerTrack()->dxy( pos ) ) >= maxXY )
         return false;
    if ( fabs( p->innerTrack()->dz ( pos ) ) >= maxZ )
         return false;
    return true;
  }
 private:
  const reco::Vertex* pv;
  int cutTL;
  int cutPL;
  float maxXY;
  float maxZ;
  bool gM;
  bool hP;
};


class BPHDaughterSelect: public lbHistoSpecificDecay::CandidateSelect {
 public:
  BPHDaughterSelect( float  ptMinLoose,
                     float  ptMinTight,
                     float etaMaxLoose,
                     float etaMaxTight,
                     const BPHSoftMuonSelect*
                           softMuonselector = 0 ): pLMin(  ptMinLoose ),
                                                   pTMin(  ptMinTight ),
                                                   eLMax( etaMaxLoose ),
                                                   eTMax( etaMaxTight ),
                                                   sms( softMuonselector ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pv = 0 ) const {
    const reco::Candidate* dptr0 = cand.daughter( 0 );
    const reco::Candidate* dptr1 = cand.daughter( 1 );
    if ( dptr0 == 0 ) return false;
    if ( dptr1 == 0 ) return false;
    float pt0 = dptr0->pt();
    float pt1 = dptr1->pt();
    if ( ( pt0 < pLMin ) || ( pt1 < pLMin ) ) return false;
    if ( ( pt0 < pTMin ) && ( pt1 < pTMin ) ) return false;
    float eta0 = fabs( dptr0->eta() );
    float eta1 = fabs( dptr1->eta() );
    if (   ( eLMax > 0 ) &&
         ( ( eta0 > eLMax ) || ( eta1 > eLMax ) ) ) return false;
    if (   ( eTMax > 0 ) && 
         ( ( eta0 > eTMax ) && ( eta1 > eTMax ) ) ) return false;
    if ( sms != 0 ) {
      const reco::Vertex* pvtx = BPHUserData::getByRef
           <reco::Vertex>( cand, "primaryVertex" );
      if ( pvtx == 0 ) return false;
      if ( !sms->accept( *dptr0, pvtx ) ) return false;
      if ( !sms->accept( *dptr1, pvtx ) ) return false;
    }    
    return true;
  }
 private:
  float pLMin;
  float pTMin;
  float eLMax;
  float eTMax;
  const BPHSoftMuonSelect* sms;
};


class BPHCompositeBasicSelect: public lbHistoSpecificDecay::CandidateSelect {
 public:
  BPHCompositeBasicSelect( float     massMin,
                           float     massMax,
                           float       ptMin = -1.0,
                           float      etaMax = -1.0,
                           float rapidityMax = -1.0 ): mMin(     massMin ),
                                                       mMax(     massMax ),
                                                       pMin(       ptMin ),
                                                       eMax(      etaMax ),
                                                       yMax( rapidityMax ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pv = 0 ) const {
    if ( ( ( mMin > 0 ) && ( mMax < 0 ) ) ||
         ( ( mMin < 0 ) && ( mMax > 0 ) ) ||
         ( ( mMin > 0 ) && ( mMax > 0 ) && ( mMin < mMax ) ) ) {
      float mass = cand.mass();
      if (   mass < mMin   ) return false;
      if ( ( mMax >    0 ) &&
           ( mass > mMax ) ) return false;
    }
    if (         cand.      pt()   < pMin   ) return false;
    if ( ( eMax                    >    0 ) &&
         ( fabs( cand.     eta() ) > eMax ) ) return false;
    if ( ( yMax                    >    0 ) &&
         ( fabs( cand.rapidity() ) > yMax ) ) return false;
    return true;
  }

 private:
  float mMin;
  float mMax;
  float pMin;
  float eMax;
  float yMax;
};


class BPHFittedBasicSelect: public lbHistoSpecificDecay::CandidateSelect {
 public:
  BPHFittedBasicSelect( float     massMin,
                        float     massMax,
                        float       ptMin = -1.0,
                        float      etaMax = -1.0,
                        float rapidityMax = -1.0 ): mMin(     massMin ),
                                                    mMax(     massMax ),
                                                    pMin(       ptMin ),
                                                    eMax(      etaMax ),
                                                    rMax( rapidityMax ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pv = 0 ) const {
    if ( !cand.hasUserFloat( "fitMass"     ) ) return false;
    float mass = cand.userFloat( "fitMass" );
    if ( ( ( mMin > 0 ) && ( mMax < 0 ) ) ||
         ( ( mMin < 0 ) && ( mMax > 0 ) ) ||
         ( ( mMin > 0 ) && ( mMax > 0 ) && ( mMin < mMax ) ) ) {
      if (   mass < mMin   ) return false;
      if ( ( mMax >    0 ) &&
           ( mass > mMax ) ) return false;
    }
    const Vector3DBase<float,GlobalTag>* fmom = BPHUserData::get
        < Vector3DBase<float,GlobalTag> >( cand, "fitMomentum" );
    if ( fmom == 0 ) return false;
    if ( pMin > 0 ) {
      if ( fmom->transverse() < pMin ) return false;
    }
    if ( eMax > 0 ) {
      if ( fabs( fmom->eta() ) > eMax ) return false;
    }
    if ( rMax > 0 ) {
      float x = fmom->x();
      float y = fmom->y();
      float z = fmom->z();
      float e = sqrt( ( x * x ) + ( y * y ) + ( z * z ) + ( mass * mass ) );
      float r = log( ( e + z ) / ( e - z ) ) / 2;
      if ( fabs( r ) > rMax ) return false;
    }
    return true;
  }

 private:
  float mMin;
  float mMax;
  float pMin;
  float eMax;
  float rMax;
};


class BPHCompositeVertexSelect: public lbHistoSpecificDecay::CandidateSelect {
 public:
  BPHCompositeVertexSelect( float probMin,
                            float  cosMin = -1.0,
                            float  sigMin = -1.0 ): pMin( probMin ),
                                                    cMin(  cosMin ),
                                                    sMin(  sigMin ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pvtx = 0 ) const {
    const reco::Vertex* svtx = BPHUserData::get
         <reco::Vertex>( cand, "vertex" );
    if ( svtx == 0 ) return false;
    if ( pvtx == 0 ) return false;
    if ( pMin > 0 ) {
      if ( ChiSquaredProbability( svtx->chi2(),
                                  svtx->ndof() ) < pMin ) return false;
    }
    if ( ( cMin > 0 ) || ( sMin > 0 ) ) {
      TVector3 disp( svtx->x() - pvtx->x(),
                     svtx->y() - pvtx->y(),
                     0 );
      TVector3 cmom( cand.px(), cand.py(), 0 );
      float cosAlpha = disp.Dot( cmom ) / ( disp.Perp() * cmom.Perp() );
      if ( cosAlpha < cMin ) return false;
      if ( sMin < 0 ) return true;
      float mass = cand.mass();
      AlgebraicVector3 vmom( cand.px(), cand.py(), 0 );
      VertexDistanceXY vdistXY;
      Measurement1D distXY = vdistXY.distance( *svtx, *pvtx );
      double ctauPV = distXY.value() * cosAlpha * mass / cmom.Perp();
      GlobalError sve = svtx->error();
      GlobalError pve = pvtx->error();
      AlgebraicSymMatrix33 vXYe = sve.matrix() + pve.matrix();
      double ctauErrPV = sqrt( ROOT::Math::Similarity( vmom, vXYe ) ) * mass /
                               cmom.Perp2();
      if ( ( ctauPV / ctauErrPV ) < sMin ) return false;
    }
    return true;
  }

 private:
  float pMin;
  float cMin;
  float sMin;
};


class BPHFittedVertexSelect: public lbHistoSpecificDecay::CandidateSelect {
 public:
  BPHFittedVertexSelect( float probMin,
                         float  cosMin = -1.0,
                         float  sigMin = -1.0 ): pMin( probMin ),
                                                 cMin(  cosMin ),
                                                 sMin(  sigMin ) {
  }
  bool accept( const pat::CompositeCandidate& cand,
               const reco::Vertex* pvtx ) const {
    const reco::Vertex* svtx = BPHUserData::get
         <reco::Vertex>( cand, "fitVertex" );
    if ( svtx == 0 ) return false;
    if ( pMin > 0 ) {
      if ( ChiSquaredProbability( svtx->chi2(),
                                  svtx->ndof() ) < pMin ) return false;
    }
    if ( ( cMin > 0 ) || ( sMin > 0 ) ) {
      TVector3 disp( svtx->x() - pvtx->x(),
                     svtx->y() - pvtx->y(),
                     0 );
      const Vector3DBase<float,GlobalTag>* fmom = BPHUserData::get
          < Vector3DBase<float,GlobalTag> >( cand, "fitMomentum" );
      if ( fmom == 0 ) return false;
      TVector3 cmom( fmom->x(), fmom->y(), 0 );
      float cosAlpha = disp.Dot( cmom ) / ( disp.Perp() * cmom.Perp() );
      if ( cosAlpha < cMin ) return false;
      if ( sMin < 0 ) return true;
      if ( !cand.hasUserFloat( "fitMass" ) ) return false;
      float mass = cand.userFloat( "fitMass" );
      AlgebraicVector3 vmom( fmom->x(), fmom->y(), 0 );
      VertexDistanceXY vdistXY;
      Measurement1D distXY = vdistXY.distance( *svtx, *pvtx );
      double ctauPV = distXY.value() * cosAlpha * mass / cmom.Perp();
      GlobalError sve = svtx->error();
      GlobalError pve = pvtx->error();
      AlgebraicSymMatrix33 vXYe = sve.matrix() + pve.matrix();
      double ctauErrPV = sqrt( ROOT::Math::Similarity( vmom, vXYe ) ) * mass /
                               cmom.Perp2();
      if ( ( ctauPV / ctauErrPV ) < sMin ) return false;
    }
    return true;
  }

 private:
  float pMin;
  float cMin;
  float sMin;
};
// useful classes end }}}

lbHistoSpecificDecay::lbHistoSpecificDecay( const edm::ParameterSet& ps ) {

  useOnia = ( SET_LABEL( oniaCandsLabel, ps ) != "" );
  useLam0 = ( SET_LABEL( lam0CandsLabel, ps ) != "" );
  useTkTk = ( SET_LABEL( tktkCandsLabel, ps ) != "" );
  useLbL0 = ( SET_LABEL( LbL0CandsLabel, ps ) != "" );
  useLbTk = ( SET_LABEL( LbTkCandsLabel, ps ) != "" );
  if ( useOnia ) consume< vector<pat::CompositeCandidate> >( oniaCandsToken,
                                                             oniaCandsLabel );
  if ( useLam0 ) consume< vector<pat::CompositeCandidate> >( lam0CandsToken,
                                                             lam0CandsLabel );
  if ( useTkTk ) consume< vector<pat::CompositeCandidate> >( tktkCandsToken,
                                                             tktkCandsLabel );
  if ( useLbL0 ) consume< vector<pat::CompositeCandidate> >( LbL0CandsToken,
                                                             LbL0CandsLabel );
  if ( useLbTk ) consume< vector<pat::CompositeCandidate> >( LbTkCandsToken,
                                                             LbTkCandsLabel );

  static const BPHSoftMuonSelect sms;
}


lbHistoSpecificDecay::~lbHistoSpecificDecay() {

}


void lbHistoSpecificDecay::fillDescriptions(
                            edm::ConfigurationDescriptions& descriptions ) {
   edm::ParameterSetDescription desc;
   desc.add<string>( "oniaCandsLabel", "" );
   desc.add<string>( "lam0CandsLabel", "" );
   desc.add<string>( "tktkCandsLabel", "" );
   desc.add<string>( "LbL0CandsLabel", "" );
   desc.add<string>( "LbTkCandsLabel", "" );
   descriptions.add( "process.lbHistoSpecificDecay", desc );
   return;
}


void lbHistoSpecificDecay::beginJob() {
  createHisto( "massLam0"   ,  40, 1.10, 1.20 ); // Phi  mass
  createHisto( "massTkTk"   , 150, 1.00, 2.50 ); // Phi  mass
  createHisto( "massJPsi"   ,  35, 2.95, 3.30 ); // JPsi mass
  createHisto( "massPsi2"   ,  60, 3.40, 4.00 ); // Psi2 mass
  createHisto( "massUps123" , 125, 8.50, 11.0 ); // Ups  mass
  createHisto( "massLbL0"   , 100, 5.00, 6.00 ); // LbL0   mass
  createHisto( "massLbTk"   , 100, 5.00, 6.00 ); // LbTk   mass
  createHisto( "mfitLbL0"   , 100, 5.00, 6.00 ); // LbL0   mass, with constraint
  createHisto( "mfitLbTk"   , 100, 5.00, 6.00 ); // LbTk   mass, with constraint
  createHisto( "massLbL0JPsi",  35, 2.95, 3.30 ); // JPsi mass in LbL0 decay
  createHisto( "massLbTkJPsi",  35, 2.95, 3.30 ); // JPsi mass in LbTk decay
  createHisto( "massLbL0Lam0",  40, 1.10, 1.20 ); // Phi  mass in Bs decay
  createHisto( "massLbL0TkTk", 150, 1.00, 2.50 ); // Kx0  mass in LbTk decay

  createHisto( "massFull"   , 200, 2.00, 12.0 ); // Full onia mass

  return;
}

void lbHistoSpecificDecay::analyze( const edm::Event& ev,
                                     const edm::EventSetup& es ) {

  // get magnetic field
  edm::ESHandle<MagneticField> magneticField;
  es.get<IdealMagneticFieldRecord>().get( magneticField );

  // get object collections
  // collections are got through "BPHTokenWrapper" interface to allow
  // uniform access in different CMSSW versions

  //////////// quarkonia ////////////

  edm::Handle< vector<pat::CompositeCandidate> > oniaCands;
  int iqo;
  int nqo = 0;
  if ( useOnia ) {
    oniaCandsToken.get( ev, oniaCands );
    nqo = oniaCands->size();
    std::cout << "nqo = " << nqo << std::endl;
  }

  for ( iqo = 0; iqo < nqo; ++ iqo ) {
    LogTrace( "DataDump" )
           << "*********** quarkonium " << iqo << "/" << nqo << " ***********";
    const pat::CompositeCandidate& cand = oniaCands->at( iqo );
    //if ( !oniaVertexSelect->accept( cand,
    //                                BPHUserData::getByRef<reco::Vertex>( cand,
    //                                "primaryVertex" ) ) ) continue;
    //if ( !oniaDaughterSelect->accept( cand ) ) continue;
    fillHisto( "Full", cand );
    fillHisto( "Phi"   , cand );
    fillHisto( "JPsi"  , cand );
    fillHisto( "Psi2"  , cand );
    fillHisto( "Ups123", cand );
  }

  //////////// LbL0 ////////////

  edm::Handle< vector<pat::CompositeCandidate> > lbLbCands;
  int ilbLb;
  int nlbLb = 0;
  if ( useLbL0 ) {
    LbL0CandsToken.get( ev, lbLbCands );
    nlbLb = lbLbCands->size();
    std::cout << "nlbLb = " << nlbLb  << std::endl;
  }

  for ( ilbLb = 0; ilbLb < nlbLb; ++ ilbLb ) {
    LogTrace( "DataDump" )
           << "*********** LbL0 " << ilbLb << "/" << nlbLb << " ***********";
    const pat::CompositeCandidate& cand = lbLbCands->at( ilbLb );
    const pat::CompositeCandidate* jPsi = BPHUserData::getByRef
         <pat::CompositeCandidate>( cand, "refToJPsi" );
    LogTrace( "DataDump" )
           << "JPsi: " << jPsi;
    if ( jPsi == 0 ) continue;
    const pat::CompositeCandidate* lam0 = BPHUserData::getByRef
         <pat::CompositeCandidate>( cand, "refToLam0" );
    LogTrace( "DataDump" )
           << "Lam0: " << lam0;
    if ( lam0 == 0 ) continue;
    fillHisto( "LbL0"    ,  cand );
    fillHisto( "LbL0JPsi", *jPsi );
    fillHisto( "LbL0Lam0" , *lam0  );
  }

  //////////// LbTk ////////////

  edm::Handle< vector<pat::CompositeCandidate> > lbTkCands;
  int ilbTk;
  int nlbTk = 0;
  if ( useLbTk ) {
    LbTkCandsToken.get( ev, lbTkCands );
    nlbTk = lbTkCands->size();
    std::cout << "nlbTk = " << nlbTk << std::endl;
  }

  for ( ilbTk = 0; ilbTk < nlbTk; ++ ilbTk ) {
    LogTrace( "DataDump" )
           << "*********** LbTk " << ilbTk << "/" << nlbTk << " ***********";
    const pat::CompositeCandidate& cand = lbTkCands->at( ilbTk );
    const pat::CompositeCandidate* jPsi = BPHUserData::getByRef
         <pat::CompositeCandidate>( cand, "refToJPsi" );
    LogTrace( "DataDump" )
           << "JPsi: " << jPsi;
    if ( jPsi == 0 ) continue;
    fillHisto( "LbTk"    ,  cand );
    fillHisto( "LbTkJPsi", *jPsi );
  }
  //////////// TkTk ////////////

  edm::Handle< vector<pat::CompositeCandidate> > tktkCands;
  int itktk;
  int ntktk = 0;
  if ( useTkTk ) {
    tktkCandsToken.get( ev, tktkCands );
    ntktk = tktkCands->size();
    std::cout << "ntktk = " << ntktk << std::endl;
  }

  for ( itktk = 0; itktk < ntktk; ++ itktk ) {
    LogTrace( "DataDump" )
           << "*********** TkTk " << itktk << "/" << ntktk << " ***********";
    const pat::CompositeCandidate& cand = tktkCands->at( itktk );
    fillHisto( "TkTk"    ,  cand );
  }
  //////////// Lam0 ////////////

  edm::Handle< vector<pat::CompositeCandidate> > lam0Cands;
  int ilam0;
  int nlam0 = 0;
  if ( useLam0 ) {
    lam0CandsToken.get( ev, lam0Cands );
    nlam0 = lam0Cands->size();
    std::cout << "nlam0 = " << nlam0 << std::endl;
  }

  for ( ilam0 = 0; ilam0 < nlam0; ++ ilam0 ) {
    LogTrace( "DataDump" )
           << "*********** Lam0 " << ilam0 << "/" << nlam0 << " ***********";
    const pat::CompositeCandidate& cand = lam0Cands->at( ilam0 );
    fillHisto( "Lam0"    ,  cand );
  }
  return;

}


void lbHistoSpecificDecay::endJob() {
  return;
}


void lbHistoSpecificDecay::fillHisto( const string& name,
                                       const pat::CompositeCandidate& cand ) {
  float mass = ( cand.hasUserFloat( "fitMass" ) ?
                 cand.   userFloat( "fitMass" ) : -1 );
  fillHisto( "mass" + name, cand.mass() );
  fillHisto( "mfit" + name, mass );
  return;
}


void lbHistoSpecificDecay::fillHisto( const string& name, float x ) {
  map<string,TH1F*>::iterator iter = histoMap.find( name );
  map<string,TH1F*>::iterator iend = histoMap.end();
  if ( iter == iend ) return;
  iter->second->Fill( x );
  return;
}


void lbHistoSpecificDecay::createHisto( const string& name,
                                         int nbin, float hmin, float hmax ) {
  histoMap[name] = fs->make<TH1F>( name.c_str(), name.c_str(),
                                   nbin, hmin, hmax );
  return;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( lbHistoSpecificDecay );
