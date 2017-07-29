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
#include "BPHAnalysis/SpecificDecay/interface/BPHPentaQToJPsiPBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleNeutralVeto.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticlePtSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleEtaSelect.h"
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
BPHPentaQToJPsiPBuilder::BPHPentaQToJPsiPBuilder( const edm::EventSetup& es,
    const std::vector<BPHPlusMinusConstCandPtr>&   JpsiCollection,
    const BPHRecoBuilder::BPHGenericCollection*  ProtonCollection ):
    jPsiName(   "JPsi" ),
  protonName( "Proton" ),
  evSetup( &es ),
  jCollection( &  JpsiCollection ),
  pCollection(  ProtonCollection )  {
  jpsiSel = new BPHMassSelect       ( 2.80, 3.40 );
// knVeto = new BPHParticleNeutralVeto;
    ptSel = new BPHParticlePtSelect (  0.7 );
   etaSel = new BPHParticleEtaSelect( 10.0 );
  massSel = new BPHMassSelect       ( 3.50, 5.50 );
  chi2Sel = new BPHChi2Select       ( 0.02 );
  mFitSel = new BPHMassFitSelect    ( jPsiName,
                                      BPHParticleMasses::Psi1Mass,
                                      BPHParticleMasses::Psi1MWidth,
                                      5.00, 6.00 );
  massConstr = true;
  minPDiff = 1.0e-4;
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHPentaQToJPsiPBuilder::~BPHPentaQToJPsiPBuilder() {
  delete jpsiSel;
//delete  knVeto;
  delete   ptSel;
  delete  etaSel;
  delete massSel;
  delete chi2Sel;
  delete mFitSel;
}

//--------------
// Operations --
//--------------
vector<BPHRecoConstCandPtr> BPHPentaQToJPsiPBuilder::build() {

  if ( updated ) return pentaQList;

  BPHRecoBuilder bPentaQ( *evSetup );
  bPentaQ.setMinPDiffererence( minPDiff );
  bPentaQ.add(   jPsiName, *jCollection );
  bPentaQ.add( protonName,  pCollection, BPHParticleMasses::protonMass,
                                         BPHParticleMasses::protonMSigma );
  bPentaQ.filter(   jPsiName, *jpsiSel );
//bPentaQ.filter(   kaonName, * knVeto );
  bPentaQ.filter( protonName, *  ptSel );
  bPentaQ.filter( protonName, * etaSel );

  bPentaQ.filter( *massSel );
  bPentaQ.filter( *chi2Sel );
  if ( massConstr ) bPentaQ.filter( *mFitSel );

  pentaQList = BPHRecoCandidate::build( bPentaQ );
//
//  Apply kinematic constraint on the JPsi mass.
//  The operation is already performed when apply the mass selection,
//  so it's not repeated. The following code is left as example
//  for similar operations
//
//  int iPentaQ;
//  int nPentaQ = ( massConstr ? pentaQList.size() : 0 );
//  for ( iPentaQ = 0; iPentaQ < nPentaQ; ++iPentaQ ) {
//    BPHRecoCandidate* cptr( const_cast<BPHRecoCandidate*>(
//                            pentaQList[iPentaQ].get() ) );
//    BPHRecoConstCandPtr jpsi = cptr->getComp( jPsiName );
//    double jMass = jpsi->constrMass();
//    if ( jMass < 0 ) continue;
//    double sigma = jpsi->constrSigma();
//    cptr->kinematicTree( jPsiName, jMass, sigma );
//  }
  updated = true;
  return pentaQList;

}

/// set cuts
void BPHPentaQToJPsiPBuilder::setJPsiMassMin( double m ) {
  updated = false;
  jpsiSel->setMassMin( m );
  return;
}


void BPHPentaQToJPsiPBuilder::setJPsiMassMax( double m ) {
  updated = false;
  jpsiSel->setMassMax( m );
  return;
}


void BPHPentaQToJPsiPBuilder::setPPtMin( double pt ) {
  updated = false;
  ptSel->setPtMin( pt );
  return;
}


void BPHPentaQToJPsiPBuilder::setPEtaMax( double eta ) {
  updated = false;
  etaSel->setEtaMax( eta );
  return;
}


void BPHPentaQToJPsiPBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHPentaQToJPsiPBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHPentaQToJPsiPBuilder::setProbMin( double p ) {
  updated = false;
  chi2Sel->setProbMin( p );
  return;
}


void BPHPentaQToJPsiPBuilder::setMassFitMin( double m ) {
  updated = false;
  mFitSel->setMassMin( m );
  return;
}


void BPHPentaQToJPsiPBuilder::setMassFitMax( double m ) {
  updated = false;
  mFitSel->setMassMax( m );
  return;
}


void BPHPentaQToJPsiPBuilder::setConstr( bool flag ) {
  updated = false;
  massConstr = flag;
  return;
}

/// get current cuts
double BPHPentaQToJPsiPBuilder::getJPsiMassMin() const {
  return jpsiSel->getMassMin();
}


double BPHPentaQToJPsiPBuilder::getJPsiMassMax() const {
  return jpsiSel->getMassMax();
}


double BPHPentaQToJPsiPBuilder::getPPtMin() const {
  return ptSel->getPtMin();
}


double BPHPentaQToJPsiPBuilder::getPEtaMax() const {
  return etaSel->getEtaMax();
}


double BPHPentaQToJPsiPBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHPentaQToJPsiPBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHPentaQToJPsiPBuilder::getProbMin() const {
  return chi2Sel->getProbMin();
}


double BPHPentaQToJPsiPBuilder::getMassFitMin() const {
  return mFitSel->getMassMin();
}


double BPHPentaQToJPsiPBuilder::getMassFitMax() const {
  return mFitSel->getMassMax();
}


bool BPHPentaQToJPsiPBuilder::getConstr() const {
  return massConstr;
}

