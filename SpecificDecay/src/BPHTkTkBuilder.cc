/*
 *  See header file for a description of this class.
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//-----------------------
// This Class' Header --
//-----------------------
#include "BPHAnalysis/SpecificDecay/interface/BPHTkTkBuilder.h"

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
BPHTkTkBuilder::BPHTkTkBuilder(
               const edm::EventSetup& es,
               const BPHRecoBuilder::BPHGenericCollection* ptkCollection,
               const BPHRecoBuilder::BPHGenericCollection* ntkCollection,
               const std::string ptkName, const double ptkMass, const double ptkSigma,
               const std::string ntkName, const double ntkMass, const double ntkSigma):
  evSetup( &es ),
  pCollection( ptkCollection ),
  nCollection( ntkCollection ),
  pName( ptkName ), pMass( ptkMass ), pSigma( ptkSigma ), 
  nName( ntkName ), nMass( ntkMass ), nSigma( ntkSigma ) 
{
    ptSel = new BPHParticlePtSelect (  0.8 );
   etaSel = new BPHParticleEtaSelect(  2.5 );
  massSel = new BPHMassSelect( 1.00, 2.20 );
  chi2Sel = new BPHChi2Select( 0.00005 );
  updated = false;
}

//--------------
// Destructor --
//--------------
BPHTkTkBuilder::~BPHTkTkBuilder() {
  delete   ptSel;
  delete  etaSel;
  delete massSel;
  delete chi2Sel;
}

//--------------
// Operations --
//--------------
vector<BPHPlusMinusConstCandPtr> BPHTkTkBuilder::build() {

  if ( updated ) return tktkList;

  BPHRecoBuilder bTkTk( *evSetup );
  bTkTk.add( pName, pCollection, pMass, pSigma );
  bTkTk.add( nName, nCollection, nMass, nSigma );
  bTkTk.filter( pName, * ptSel );
  bTkTk.filter( nName, * ptSel );
  bTkTk.filter( pName, *etaSel );
  bTkTk.filter( nName, *etaSel );
  //BPHMassSymSelect mTmpSel( pName,   nName, massSel );
  //bTkTk.filter( mTmpSel );

  vector<BPHPlusMinusConstCandPtr>
  tmpList = BPHPlusMinusCandidate::build( bTkTk, pName,   nName );

  return tmpList;
  int ikx;
  int nkx = tmpList.size();
  tktkList.clear();
  tktkList.reserve( nkx );
  BPHPlusMinusConstCandPtr pxt( 0 );
  for ( ikx = 0; ikx < nkx; ++ikx ) {
    BPHPlusMinusConstCandPtr& pxt = tmpList[ikx];
    if ( !massSel->accept( *pxt ) ) continue;
    if ( !chi2Sel->accept( *pxt ) ) continue;

    pxt->updateMom();
    pxt->composite();
    tktkList.push_back( pxt );
  }

  updated = true;
  return tktkList;

}

/// set cuts
void BPHTkTkBuilder::setPtMin( double pt )  {
  updated = false;
  ptSel->setPtMin( pt );
  return;
}


void BPHTkTkBuilder::setEtaMax( double eta ) {
  updated = false;
  etaSel->setEtaMax( eta );
  return;
}


void BPHTkTkBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHTkTkBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHTkTkBuilder::setProbMin( double p ) {
  updated = false;
  chi2Sel->setProbMin( p );
  return;
}


void BPHTkTkBuilder::setConstr( double mass, double sigma ) {
  updated = false;
  cMass  = mass;
  cSigma = sigma;
  return;
}

/// get current cuts
double BPHTkTkBuilder::getPtMin()  const {
  return ptSel->getPtMin();
}


double BPHTkTkBuilder::getEtaMax() const {
  return etaSel->getEtaMax();
}


double BPHTkTkBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHTkTkBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHTkTkBuilder::getProbMin() const {
  return chi2Sel->getProbMin();
}


double BPHTkTkBuilder::getConstrMass() const {
  return cMass;
}


double BPHTkTkBuilder::getConstrSigma() const {
  return cSigma;
}
