/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHKshortToPiPiBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticlePtSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleEtaSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSymSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHChi2Select.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassFitSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

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
BPHKshortToPiPiBuilder::BPHKshortToPiPiBuilder(
               const edm::EventSetup& es,
               const BPHRecoBuilder::BPHGenericCollection* posCollection,
               const BPHRecoBuilder::BPHGenericCollection* negCollection ):
  pName( "PiPos" ),
  nName( "PiNeg" ),
  evSetup( &es ),
  pCollection( posCollection ),
  nCollection( negCollection ) {
    ptSel = new BPHParticlePtSelect (  0.7 );
   etaSel = new BPHParticleEtaSelect( 10.0 );
  massSel = new BPHMassSelect( 0.3 , 1.10 );
  chi2Sel = new BPHChi2Select( 0.0 );
  mFitSel = new BPHMassFitSelect( "", -1.0, -1.0, 0.3, 1.0 );
  massConstr = false;
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHKshortToPiPiBuilder::~BPHKshortToPiPiBuilder() {
  delete   ptSel;
  delete  etaSel;
  delete massSel;
  delete mFitSel;
  delete chi2Sel;
}

//--------------
// Operations --
//--------------
vector<BPHPlusMinusConstCandPtr> BPHKshortToPiPiBuilder::build() {

  if ( updated ) return ksList;

  BPHRecoBuilder bKs( *evSetup );
  bKs.add( pName, pCollection, BPHParticleMasses::pionMass,
                               BPHParticleMasses::pionMSigma );
  bKs.add( nName, nCollection, BPHParticleMasses::pionMass,
                               BPHParticleMasses::pionMSigma );
  bKs.filter( pName, *ptSel );
  bKs.filter( nName, *ptSel );
  bKs.filter( pName, *etaSel );
  bKs.filter( nName, *etaSel );

  if ( massConstr ) bKs.filter( *mFitSel );

  ksList = BPHPlusMinusCandidate::build( bKs, pName, nName );

  updated = true;
  return ksList;

}

/// set cuts
void BPHKshortToPiPiBuilder::setPtMin( double pt ) {
  updated = false;
  ptSel->setPtMin( pt );
  return;
}


void BPHKshortToPiPiBuilder::setEtaMax( double eta ) {
  updated = false;
  etaSel->setEtaMax( eta );
  return;
}


void BPHKshortToPiPiBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHKshortToPiPiBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHKshortToPiPiBuilder::setProbMin( double p ) {
  updated = false;
  chi2Sel->setProbMin( p );
  return;
}

void BPHKshortToPiPiBuilder::setConstr( double mass, double sigma ) {
  updated = false;
  cMass  = mass;
  cSigma = sigma;
  massConstr = true;

  return;
}

/// get current cuts
double BPHKshortToPiPiBuilder::getPtMin() const {
  return ptSel->getPtMin();
}


double BPHKshortToPiPiBuilder::getEtaMax() const {
  return etaSel->getEtaMax();
}


double BPHKshortToPiPiBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHKshortToPiPiBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHKshortToPiPiBuilder::getProbMin() const {
  return chi2Sel->getProbMin();
}


double BPHKshortToPiPiBuilder::getMassFitMin() const {
  return mFitSel->getMassMin();
}


double BPHKshortToPiPiBuilder::getMassFitMax() const {
  return mFitSel->getMassMax();
}
double BPHKshortToPiPiBuilder::getConstrMass() const {
  return cMass;
}


double BPHKshortToPiPiBuilder::getConstrSigma() const {
  return cSigma;
}

