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
#include "BPHAnalysis/SpecificDecay/interface/BPHLambda0_bToJPsiTkTkBuilder.h"

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
BPHLambda0_bToJPsiTkTkBuilder::BPHLambda0_bToJPsiTkTkBuilder( const edm::EventSetup& es,
    const std::vector<BPHPlusMinusConstCandPtr>&  JpsiCollection,
    const std::vector<BPHPlusMinusConstCandPtr>&  TkTkCollection):
  jPsiName( "JPsi" ),
  tktkName( "TkTk" ),
  evSetup( &es ),
   jpsiCollection( & JpsiCollection ),
   tktkCollection( & TkTkCollection ),
  _massDiff( 0.0148 ) {
  jpsiSel = new BPHMassSelect   ( 2.80, 3.40 );
  tktkSel = new BPHMassSelect   ( 0.  , 0.   );
  massSel = new BPHMassSelect   ( 4.50, 7.00 );
  chi2Sel = new BPHChi2Select   ( 0.02 );
  ParticleMass jpsiMass( 3.0916 ); 
  jpsiConstr = new TwoTrackMassKinematicConstraint( jpsiMass );
  mFitSel = new BPHMassFitSelect( tktkName,
                                  jpsiConstr,
                                  5.00, 6.00 );
  massConstr = true;
  minPDiff = 1.0e-4;
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHLambda0_bToJPsiTkTkBuilder::~BPHLambda0_bToJPsiTkTkBuilder() {
  delete jpsiConstr;
  delete  jpsiSel;
  delete  massSel;
  delete  chi2Sel;
  delete  mFitSel;
}

//--------------
// Operations --
//--------------
vector<BPHRecoConstCandPtr> BPHLambda0_bToJPsiTkTkBuilder::build() {

  if ( updated ) return lbList;

  BPHRecoBuilder bLb( *evSetup );
  bLb.setMinPDiffererence( minPDiff );
  bLb.add( jPsiName, *jpsiCollection );
  bLb.add( tktkName, *tktkCollection );
  bLb.filter( jPsiName, *jpsiSel  );
  if ( this->getTkTkMassMin() > 0.0001 && this->getTkTkMassMax() > 0.0001 )
      bLb.filter( tktkName, *tktkSel );

  bLb.filter( *massSel );
  bLb.filter( *chi2Sel );
  //if ( massConstr ) bLb.filter( *mFitSel );
  bLb.filter( *mFitSel );

  lbList = BPHRecoCandidate::build( bLb );

  updated = true;

  return lbList;

}

/// set cuts
void BPHLambda0_bToJPsiTkTkBuilder::setJPsiMassMin( double m ) {
  updated = false;
  jpsiSel->setMassMin( m );
  return;
}


void BPHLambda0_bToJPsiTkTkBuilder::setJPsiMassMax( double m ) {
  updated = false;
  jpsiSel->setMassMax( m );
  return;
}


void BPHLambda0_bToJPsiTkTkBuilder::setTkTkMassMin( double m ) {
  updated = false;
  tktkSel->setMassMin( m );
  return;
}


void BPHLambda0_bToJPsiTkTkBuilder::setTkTkMassMax( double m ) {
  updated = false;
  tktkSel->setMassMax( m );
  return;
}


void BPHLambda0_bToJPsiTkTkBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHLambda0_bToJPsiTkTkBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHLambda0_bToJPsiTkTkBuilder::setProbMin( double p ) {
  updated = false;
  chi2Sel->setProbMin( p );
  return;
}


void BPHLambda0_bToJPsiTkTkBuilder::setMassFitMin( double m ) {
  updated = false;
  mFitSel->setMassMin( m );
  return;
}


void BPHLambda0_bToJPsiTkTkBuilder::setMassFitMax( double m ) {
  updated = false;
  mFitSel->setMassMax( m );
  return;
}


void BPHLambda0_bToJPsiTkTkBuilder::setConstr( bool flag ) {
  updated = false;
  massConstr = flag;
  return;
}

/// get current cuts
double BPHLambda0_bToJPsiTkTkBuilder::getJPsiMassMin() const {
  return jpsiSel->getMassMin();
}


double BPHLambda0_bToJPsiTkTkBuilder::getJPsiMassMax() const {
  return jpsiSel->getMassMax();
}


double BPHLambda0_bToJPsiTkTkBuilder::getTkTkMassMin() const {
  return tktkSel->getMassMax();
}


double BPHLambda0_bToJPsiTkTkBuilder::getTkTkMassMax() const {
  return tktkSel->getMassMax();
}


double BPHLambda0_bToJPsiTkTkBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHLambda0_bToJPsiTkTkBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHLambda0_bToJPsiTkTkBuilder::getProbMin() const {
  return chi2Sel->getProbMin();
}


double BPHLambda0_bToJPsiTkTkBuilder::getMassFitMin() const {
  return mFitSel->getMassMin();
}


double BPHLambda0_bToJPsiTkTkBuilder::getMassFitMax() const {
  return mFitSel->getMassMax();
}


bool BPHLambda0_bToJPsiTkTkBuilder::getConstr() const {
  return massConstr;
}

