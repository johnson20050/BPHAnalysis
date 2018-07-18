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
    const std::vector<BPHPlusMinusConstCandPtr>& pTkTkCollection,
    const std::vector<BPHPlusMinusConstCandPtr>& Collection ):
  jPsiName( "JPsi" ),
  tktkName( "TkTk" ),
  evSetup( &es ),
   jpsiCollection( & JpsiCollection ),
  ptktkCollection( &pTkTkCollection ),
  ntktkCollection( &Collection ),
  _massDiff( 0.0148 ) {
  jpsiSel = new BPHMassSelect   ( 2.80, 3.40 );
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

  //BPHRecoBuilder bLb( *evSetup );
  //bLb.setMinPDiffererence( minPDiff );
  //bLb.add( jPsiName, * jpsiCollection );
  //bLb.add( tktkName, *ptktkCollection );
  //bLb.filter( jPsiName, *jpsiSel  );

  //bLb.filter( *massSel );
  //bLb.filter( *chi2Sel );
  //if ( massConstr ) bLb.filter( *mFitSel );

  //std::vector<BPHRecoConstCandPtr> tmpList = BPHRecoCandidate::build( bLb );

  if ( ptktkCollection->size() != ntktkCollection->size() )
      printf("2 tktk list not the same!\n");
  int iJPsi;
  int nJPsi =  jpsiCollection->size();
  int iTkTk;
  int nTkTk = ptktkCollection->size();
  lbList.clear();
  lbList.reserve( nTkTk*nJPsi );
  BPHRecoConstCandPtr pxt( 0 );
  for ( iJPsi = 0; iJPsi < nJPsi; ++iJPsi )
  for ( iTkTk = 0; iTkTk < nTkTk; ++iTkTk ) {
    BPHRecoCandidatePtr px0( new BPHRecoCandidate( evSetup ) );
    BPHRecoCandidatePtr pxb( new BPHRecoCandidate( evSetup ) );
    BPHRecoCandidate* kx0 = px0.get();
    BPHRecoCandidate* kxb = pxb.get();

    //kx0->add( jPsiName, BPHRecoConstCandPtr( jpsiCollection->at(iJPsi).get()) );
    //kxb->add( jPsiName, BPHRecoConstCandPtr( jpsiCollection->at(iJPsi).get()) );
    kx0->add( jPsiName, jpsiCollection->at(iJPsi) );
    kxb->add( jPsiName, jpsiCollection->at(iJPsi) );
    //kx0->add( tktkName, BPHRecoConstCandPtr(ptktkCollection->at(iTkTk).get()) );
    //kxb->add( tktkName, BPHRecoConstCandPtr(ntktkCollection->at(iTkTk).get()) );
    kx0->add( tktkName, ptktkCollection->at(iTkTk) );
    kxb->add( tktkName, ntktkCollection->at(iTkTk) );

    if ( massConstr )
        if ( !mFitSel->accept( *px0 ) ) continue;
    if ( massConstr )
        if ( !mFitSel->accept( *pxb ) ) continue;
    if ( !kx0->isValidFit() ) continue;
    if ( !kxb->isValidFit() ) continue;
    if ( !massSel->accept( *px0 ) ) continue;
    if ( !massSel->accept( *pxb ) ) continue;
    if ( !chi2Sel->accept( *px0 ) ) continue;
    if ( !chi2Sel->accept( *pxb ) ) continue;
    float mass0 = kx0->currentParticle()->currentState().mass();
    float massb = kxb->currentParticle()->currentState().mass();
    // if particle and anti-particle are in Lam0 signal region, it cannot be distinguished, throwout it.
    if ( fabs( mass0 - BPHParticleMasses::Lambda0_bMass ) < _massDiff &&
         fabs( massb - BPHParticleMasses::Lambda0_bMass ) < _massDiff   )
        continue;
    else if
       ( fabs( mass0 - BPHParticleMasses::Lambda0_bMass ) <
         fabs( massb - BPHParticleMasses::Lambda0_bMass ) )
       {
           lbList.push_back( px0 );
           ltktkList.push_back( ptktkCollection->at(iTkTk) );
       }
    else
    {
        lbList.push_back( pxb );
        ltktkList.push_back( ntktkCollection->at(iTkTk) );
    }
  }
  updated = true;

  return lbList;

}
std::vector<BPHPlusMinusConstCandPtr> BPHLambda0_bToJPsiTkTkBuilder::getTkTkList()
{ if ( !updated ) build(); return ltktkList; }


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

