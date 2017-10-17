/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHLambda0ToPPiBuilder.h"

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
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassFitSelect.h"
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
BPHLambda0ToPPiBuilder::BPHLambda0ToPPiBuilder(
               const edm::EventSetup& es,
               const BPHRecoBuilder::BPHGenericCollection* kaonCollection,
               const BPHRecoBuilder::BPHGenericCollection* pionCollection ):
  kaonName( "Proton" ),
  pionName( "Pion" ),
  evSetup( &es ),
  kCollection( kaonCollection ),
  pCollection( pionCollection ) {
    ptSel = new BPHParticlePtSelect (  0.7 );
   etaSel = new BPHParticleEtaSelect( 10.0 );
  massSel = new BPHMassSelect( 1.00, 1.20 );
  chi2Sel = new BPHChi2Select( 0.0 );
  mFitSel = new BPHMassFitSelect( "", 1.115683, -1.0, 1.00, 1.20 );
  massConstr = false;

  updated = false;
}

//--------------
// Destructor --
//--------------
BPHLambda0ToPPiBuilder::~BPHLambda0ToPPiBuilder() {
  delete   ptSel;
  delete  etaSel;
  delete massSel;
  delete chi2Sel;
}

//--------------
// Operations --
//--------------
vector<BPHPlusMinusConstCandPtr> BPHLambda0ToPPiBuilder::build() {

  if ( updated ) return kx0List;

  BPHRecoBuilder bKx0( *evSetup );
  bKx0.add( kaonName, kCollection, BPHParticleMasses::protonMass,
                                   BPHParticleMasses::protonMSigma );
  bKx0.add( pionName, pCollection, BPHParticleMasses::pionMass,
                                   BPHParticleMasses::pionMSigma );
  bKx0.filter( kaonName, *ptSel );
  bKx0.filter( pionName, *ptSel );
  bKx0.filter( kaonName, *etaSel );
  bKx0.filter( pionName, *etaSel );
  BPHMassSymSelect mTmpSel( kaonName, pionName, massSel );
  bKx0.filter( mTmpSel );

  vector<BPHPlusMinusConstCandPtr>
  tmpList = BPHPlusMinusCandidate::build( bKx0, kaonName, pionName );

  int ikx;
  int nkx = tmpList.size();
  kx0List.clear();
  kx0List.reserve( nkx );
  BPHPlusMinusConstCandPtr pxt( 0 );
  for ( ikx = 0; ikx < nkx; ++ikx ) {
    BPHPlusMinusConstCandPtr& px0 = tmpList[ikx];
    BPHPlusMinusCandidatePtr  pxb( new BPHPlusMinusCandidate( evSetup ) );
    const
    BPHPlusMinusCandidate* kx0 = px0.get();
    BPHPlusMinusCandidate* kxb = pxb.get();
    kxb->add( pionName, kx0->originalReco( kx0->getDaug( kaonName ) ),
              BPHParticleMasses::pionMass );
    kxb->add( kaonName, kx0->originalReco( kx0->getDaug( pionName ) ),
              BPHParticleMasses::protonMass );
    if ( fabs( kx0->composite().mass() - BPHParticleMasses::lambda0Mass ) <
         fabs( kxb->composite().mass() - BPHParticleMasses::lambda0Mass ) )
         pxt = px0;
    else pxt = pxb;
    if ( !massSel->accept( *pxt ) ) continue;
    if ( !chi2Sel->accept( *pxt ) ) continue;

    kx0List.push_back( pxt );
  }

  updated = true;

  return kx0List;

}

/// set cuts
void BPHLambda0ToPPiBuilder::setPtMin( double pt ) {
  updated = false;
  ptSel->setPtMin( pt );
  return;
}


void BPHLambda0ToPPiBuilder::setEtaMax( double eta ) {
  updated = false;
  etaSel->setEtaMax( eta );
  return;
}


void BPHLambda0ToPPiBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHLambda0ToPPiBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHLambda0ToPPiBuilder::setProbMin( double p ) {
  updated = false;
  chi2Sel->setProbMin( p );
  return;
}


void BPHLambda0ToPPiBuilder::setConstr( double mass, double sigma ) {
  updated = false;
  cMass  = mass;
  cSigma = sigma;
  massConstr = true;

  return;
}

/// get current cuts
double BPHLambda0ToPPiBuilder::getPtMin() const {
  return ptSel->getPtMin();
}


double BPHLambda0ToPPiBuilder::getEtaMax() const {
  return etaSel->getEtaMax();
}


double BPHLambda0ToPPiBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHLambda0ToPPiBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHLambda0ToPPiBuilder::getProbMin() const {
  return chi2Sel->getProbMin();
}


double BPHLambda0ToPPiBuilder::getConstrMass() const {
  return cMass;
}


double BPHLambda0ToPPiBuilder::getConstrSigma() const {
  return cSigma;
}

