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

// select {{{
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
// select }}}

lbHistoSpecificDecay::lbHistoSpecificDecay( const edm::ParameterSet& ps ) {

  useOnia = ( SET_LABEL( oniaCandsLabel, ps ) != "" );
  useLbToJPsiLam0   = ( SET_LABEL(   LbToJPsiLam0CandsLabel, ps ) != "" );
  useLbToJPsiTkTk   = ( SET_LABEL(   LbToJPsiTkTkCandsLabel, ps ) != "" );
  if ( useOnia           ) consume< vector<pat::CompositeCandidate> >(           oniaCandsToken,
                                                                                 oniaCandsLabel );
  if ( useLbToJPsiLam0   ) consume< vector<pat::CompositeCandidate> >(   LbToJPsiLam0CandsToken,
                                                                         LbToJPsiLam0CandsLabel );
  if ( useLbToJPsiTkTk   ) consume< vector<pat::CompositeCandidate> >(   LbToJPsiTkTkCandsToken,
                                                                         LbToJPsiTkTkCandsLabel );

}


lbHistoSpecificDecay::~lbHistoSpecificDecay() {
}


void lbHistoSpecificDecay::fillDescriptions(
                            edm::ConfigurationDescriptions& descriptions ) {
   edm::ParameterSetDescription desc;
   desc.add<string>( "oniaCandsLabel", "" );
   desc.add<string>(   "LbToJPsiLam0CandsLabel", "" );
   desc.add<string>(   "LbToJPsiTkTkCandsLabel", "" );
   descriptions.add( "process.lbHistoSpecificDecay", desc );
   return;
}


void lbHistoSpecificDecay::beginJob() {

          createHisto( "mass_LbToJPsiLam0"        , 50, 5., 6.);
          createHisto( "mass_LbToJPsiLam0_refit"  , 50, 5., 6.);
          createHisto( "mass_JPsiinLbLam0"        , 35,2.9,3.3);
          createHisto( "mass_Lam0inLbLam0"        , 50,1.0,1.2);
          createHisto( "mass_LbToJPsiTkTk"        , 50, 5., 6.);
          createHisto( "mass_LbToJPsiTkTk_refit"  , 50, 5., 6.);
          createHisto( "mass_JPsiinLbTkTk"        , 35,2.9,3.3);

  return;
}

void lbHistoSpecificDecay::analyze( const edm::Event& ev,
                                     const edm::EventSetup& es ) {


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
  }

  for ( iqo = 0; iqo < nqo; ++ iqo ) {
    LogTrace( "DataDump" )
           << "*********** quarkonium " << iqo << "/" << nqo << " ***********";
    const pat::CompositeCandidate& cand = oniaCands->at( iqo );
    fillHisto( "Full", cand );
    fillHisto( "Phi"   , cand );
    fillHisto( "JPsi"  , cand );
    fillHisto( "Psi2"  , cand );
    fillHisto( "Ups123", cand );
  }
    
  //////// LbToJPsiLam0 //////
  edm::Handle< vector<pat::CompositeCandidate> > lbToJPsiLam0Cands;
  int ilb1;
  int nlb1;
  if ( useLbToJPsiLam0 )
  {
      LbToJPsiLam0CandsToken.get( ev, lbToJPsiLam0Cands );
      nlb1 = lbToJPsiLam0Cands->size();
  }
  for ( ilb1 = 0; ilb1 < nlb1; ++nlb1 )
  {
      const pat::CompositeCandidate& cand = lbToJPsiLam0Cands->at( ilb1 );
      const pat::CompositeCandidate* jpsi = BPHUserData::getByRef< pat::CompositeCandidate>( cand, "refToJPsi" );
      const pat::CompositeCandidate* lam0 = BPHUserData::getByRef< pat::CompositeCandidate>( cand, "refToLam0" );
      
      fillHisto( "mass_LbToJPsiLam0"        , cand.mass() );
      if ( cand.hasUserFloat( "fitMass" ) )
          fillHisto( "mass_LbToJPsiLam0_refit"  , cand.userFloat("fitMass") );
      if ( jpsi )
          fillHisto( "mass_JPsiinLbLam0"        , jpsi->mass() );
      if ( lam0 )
          fillHisto( "mass_Lam0inLbLam0"        , lam0->mass() );
  }


  //////// LbToJPsiTkTk ///////
  edm::Handle< vector<pat::CompositeCandidate> > lbToJPsiTkTkCands;
  int ilb2;
  int nlb2 =0;
  if ( useLbToJPsiTkTk )
  {
      LbToJPsiTkTkCandsToken.get( ev, lbToJPsiTkTkCands );
      nlb2 = lbToJPsiTkTkCands->size();
  }
  for ( ilb2 = 0;ilb2 < nlb2; ++ilb2 )
  {
      const pat::CompositeCandidate& cand = lbToJPsiTkTkCands->at( ilb2 );
      const pat::CompositeCandidate* jpsi = BPHUserData::getByRef< pat::CompositeCandidate>( cand, "refToJPsi" );
      const reco::Candidate* kptr = BPHDaughters::get( cand, 0.49, 0.51 ).front();
      const reco::Candidate* pptr = BPHDaughters::get( cand, 0.99, 1.01 ).front();
      if ( kptr->charge() * pptr->charge() == 1 ) printf( "error------------Lb is not neutral1!\n" );
      if ( kptr->charge() * pptr->charge() == 0 ) printf( "error------------Lb is not neutral2!\n" );
        
      fillHisto( "mass_LbToJPsiTkTk"        , cand.mass() );
      if ( cand.hasUserFloat( "fitMass" ) )
          fillHisto( "mass_LbToJPsiTkTk_refit"  , cand.userFloat("fitMass") );
      if ( jpsi )
          fillHisto( "mass_JPsiinLbTkTk"        , jpsi->mass() );
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
