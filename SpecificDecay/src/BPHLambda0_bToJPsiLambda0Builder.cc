/*
 *  See header file for a description of this class.
 *
 *  $Date: 2015-07-24 11:29:20 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHLambda0_bToJPsiLambda0Builder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHChi2Select.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassFitSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"

//---------------
// C++ Headers --
//---------------
using namespace std;

//-------------------
// Initializations --
//-------------------


//----------------
// Constructors --
//----------------
BPHLambda0_bToJPsiLambda0Builder::BPHLambda0_bToJPsiLambda0Builder( const edm::EventSetup& es,
    const std::vector<BPHPlusMinusConstCandPtr>& JpsiCollection,
    const std::vector<BPHPlusMinusConstCandPtr>& Lam0Collection ):
  jPsiName( "JPsi" ),
  lam0Name( "Lam0" ),
  evSetup( &es ),
  jpsiCollection( &JpsiCollection ),
  lam0Collection( &Lam0Collection ) {
  jpsiSel = new BPHMassSelect   ( 2.80, 3.40 );
  mlam0Sel= new BPHMassSelect   ( 1.05, 1.20 );
  massSel = new BPHMassSelect   ( 4.50, 6.50 );
  chi2Sel = new BPHChi2Select   ( 0.02 );
  ParticleMass jpsiMass( 3.0916 ); 
  jpsiConstr = new TwoTrackMassKinematicConstraint( jpsiMass );
  mFitSel = new BPHMassFitSelect( lam0Name,
                                  jpsiConstr,
                                  5.00, 6.00 );
  massConstr = true;
  minPDiff = 1.0e-4;
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHLambda0_bToJPsiLambda0Builder::~BPHLambda0_bToJPsiLambda0Builder() {
  delete jpsiConstr;
  delete  jpsiSel;
  delete mlam0Sel;
  delete  massSel;
  delete  chi2Sel;
  delete  mFitSel;
}

//--------------
// Operations --
//--------------
vector<BPHRecoConstCandPtr> BPHLambda0_bToJPsiLambda0Builder::build() {

  if ( updated ) return lbList;

  BPHRecoBuilder bLb( *evSetup );
  bLb.setMinPDiffererence( minPDiff );
  bLb.add( jPsiName, *jpsiCollection );
  bLb.add( lam0Name, *lam0Collection );
  bLb.filter( jPsiName, *jpsiSel  );
  bLb.filter( lam0Name, *mlam0Sel );

  bLb.filter( *massSel );
  bLb.filter( *chi2Sel );
  if ( massConstr ) bLb.filter( *mFitSel );

  lbList = BPHRecoCandidate::build( bLb );
//
//  Apply kinematic constraint on the JPsi mass.
//  The operation is already performed when apply the mass selection,
//  so it's not repeated. The following code is left as example
//  for similar operations
//
//  int iLb;
//  int nLb = ( massConstr ? lbList.size() : 0 );
//  for ( iLb = 0; iLb < nLb; ++iLb ) {
//    BPHRecoCandidate* cptr( const_cast<BPHRecoCandidate*>(
//                            lbList[iLb].get() ) );
//    BPHRecoConstCandPtr jpsi = cptr->getComp( jPsiName );
//    double jMass = jpsi->constrMass();
//    if ( jMass < 0 ) continue;
//    double sigma = jpsi->constrSigma();
//    cptr->kinematicTree( virtualParticleName, jpsiMass, jpsisigma );
//  }
  updated = true;

  return lbList;

}

/// set cuts
void BPHLambda0_bToJPsiLambda0Builder::setJPsiMassMin( double m ) {
  updated = false;
  jpsiSel->setMassMin( m );
  return;
}


void BPHLambda0_bToJPsiLambda0Builder::setJPsiMassMax( double m ) {
  updated = false;
  jpsiSel->setMassMax( m );
  return;
}


void BPHLambda0_bToJPsiLambda0Builder::setLam0MassMin( double m ) {
  updated = false;
  mlam0Sel->setMassMin( m );
  return;
}


void BPHLambda0_bToJPsiLambda0Builder::setLam0MassMax( double m ) {
  updated = false;
  mlam0Sel->setMassMax( m );
  return;
}


void BPHLambda0_bToJPsiLambda0Builder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHLambda0_bToJPsiLambda0Builder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHLambda0_bToJPsiLambda0Builder::setProbMin( double p ) {
  updated = false;
  chi2Sel->setProbMin( p );
  return;
}


void BPHLambda0_bToJPsiLambda0Builder::setMassFitMin( double m ) {
  updated = false;
  mFitSel->setMassMin( m );
  return;
}


void BPHLambda0_bToJPsiLambda0Builder::setMassFitMax( double m ) {
  updated = false;
  mFitSel->setMassMax( m );
  return;
}


void BPHLambda0_bToJPsiLambda0Builder::setConstr( bool flag ) {
  updated = false;
  massConstr = flag;
  return;
}

/// get current cuts
double BPHLambda0_bToJPsiLambda0Builder::getJPsiMassMin() const {
  return jpsiSel->getMassMin();
}


double BPHLambda0_bToJPsiLambda0Builder::getJPsiMassMax() const {
  return jpsiSel->getMassMax();
}


double BPHLambda0_bToJPsiLambda0Builder::getLam0MassMin() const {
  return mlam0Sel->getMassMin();
}


double BPHLambda0_bToJPsiLambda0Builder::getLam0MassMax() const {
  return mlam0Sel->getMassMax();
}


double BPHLambda0_bToJPsiLambda0Builder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHLambda0_bToJPsiLambda0Builder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHLambda0_bToJPsiLambda0Builder::getProbMin() const {
  return chi2Sel->getProbMin();
}


double BPHLambda0_bToJPsiLambda0Builder::getMassFitMin() const {
  return mFitSel->getMassMin();
}


double BPHLambda0_bToJPsiLambda0Builder::getMassFitMax() const {
  return mFitSel->getMassMax();
}


bool BPHLambda0_bToJPsiLambda0Builder::getConstr() const {
  return massConstr;
}

