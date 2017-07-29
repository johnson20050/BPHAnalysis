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
#include "BPHAnalysis/SpecificDecay/interface/BPHLbToJPsiTkTkBuilder.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticlePtSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleEtaSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSymSelect_jpsi.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHChi2Select.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleNeutralVeto.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHCompChargeSelect.h"
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
BPHLbToJPsiTkTkBuilder::BPHLbToJPsiTkTkBuilder( const edm::EventSetup& es,
    const std::vector<BPHPlusMinusConstCandPtr>& jPsiCollection,
    const std::vector<BPHPlusMinusConstCandPtr>& TkTkCollection,
    const std::string ptkName, const double ptkMass, const double ptkSigma,
    const std::string ntkName, const double ntkMass, const double ntkSigma):
    evSetup( &es ),
    jpsiCollection( &jPsiCollection ),
    tktkCollection( &TkTkCollection ),
    jPsiName( "JPsi" ),
       pName( ptkName), pMass( ptkMass), pSigma(ptkSigma),
       nName( ntkName), nMass( ntkMass), nSigma(ntkSigma)
{
    jpsiSel = new BPHMassSelect           ( 2.80, 3.40 );
      ptSel = new BPHParticlePtSelect     ( 0.7 );
     etaSel = new BPHParticleEtaSelect    ( 2.5 ); 

    massSel = new BPHMassSelect           ( 4.50, 6.50 );
 massTmpSel = new BPHMassSymSelect_jpsi   ( pName, nName, massSel );
    chi2Sel = new BPHChi2Select           ( 0.02 );
    mFitSel = new BPHMassFitSelect( jPsiName,
                                    BPHParticleMasses::jPsiMass,
                                    BPHParticleMasses::jPsiMSigma,
                                    5.00, 6.00 );
    massConstr = false;
    minPDiff = 1.0e-4;
    updated = false;
}

//--------------
// Destructor --
//--------------
BPHLbToJPsiTkTkBuilder::~BPHLbToJPsiTkTkBuilder() {
  delete     jpsiSel;
  delete       ptSel;
  delete      etaSel;

  delete     massSel;
  delete  massTmpSel;
  delete     chi2Sel;
  delete     mFitSel;
}

//--------------
// Operations --
//--------------
vector<BPHRecoConstCandPtr> BPHLbToJPsiTkTkBuilder::build() {

    if ( updated ) return lbList;

    vector<BPHRecoConstCandPtr> tmpList;
    BPHRecoConstCandPtr pxt( 0 );
    vector<BPHPlusMinusConstCandPtr>::const_iterator jpsiIter = jpsiCollection->begin();
    vector<BPHPlusMinusConstCandPtr>::const_iterator jpsiIend = jpsiCollection->end  ();
    while ( jpsiIter != jpsiIend )
    {
        const BPHPlusMinusConstCandPtr njp = *jpsiIter;
        const BPHRecoCandidate* _jpsi = jpsiIter++->get();
        const reco::Candidate* muPos = _jpsi->originalReco(
                                       _jpsi->getDaug( "MuPos" ) );
        const reco::Candidate* muNeg = _jpsi->originalReco(
                                       _jpsi->getDaug( "MuNeg" ) );
        vector<BPHPlusMinusConstCandPtr>::const_iterator tktkIter = tktkCollection->begin();
        vector<BPHPlusMinusConstCandPtr>::const_iterator tktkIend = tktkCollection->end  ();
        while ( tktkIter != tktkIend )
        {
            const BPHRecoCandidate* _tktk = tktkIter++->get();
            const reco::Candidate* pTk = _tktk->originalReco(
                                         _tktk->getDaug( pName ) );
            const reco::Candidate* nTk = _tktk->originalReco(
                                         _tktk->getDaug( nName ) );
            if ( BPHRecoBuilder::sameTrack( muPos, pTk, 0.0005 ) ) continue;
            if ( BPHRecoBuilder::sameTrack( muPos, nTk, 0.0005 ) ) continue;
            if ( BPHRecoBuilder::sameTrack( muNeg, pTk, 0.0005 ) ) continue;
            if ( BPHRecoBuilder::sameTrack( muNeg, nTk, 0.0005 ) ) continue;
            if ( ! ptSel->accept( *pTk ) ) continue;
            if ( ! ptSel->accept( *nTk ) ) continue;
            if ( !etaSel->accept( *pTk ) ) continue;
            if ( !etaSel->accept( *nTk ) ) continue;

            BPHRecoCandidatePtr nLb1( new BPHRecoCandidate( evSetup ) );
            nLb1->add( jPsiName, njp );
            nLb1->add( pName, pTk, pMass, pSigma );
            nLb1->add( nName, nTk, nMass, nSigma );
            nLb1->kinematicTree( jPsiName,
                                 BPHParticleMasses::jPsiMass,
                                 BPHParticleMasses::jPsiMWidth );
            ////nLb1.setNotUpdated();
            //// change the order of tk1 and tk2 to find anti-particle
            //BPHRecoCandidatePtr nLb2( new BPHRecoCandidate( evSetup ) );
            //nLb2->add( jPsiName, njp );
            //nLb2->add( pName, nTk, pMass, pSigma );
            //nLb2->add( nName, pTk, nMass, nSigma );
            //nLb2->kinematicTree( jPsiName,
            //                     BPHParticleMasses::jPsiMass,
            //                     BPHParticleMasses::jPsiMWidth );
            ////nLb2.setNotUpdated();
            //if ( fabs( nLb1->composite().mass() - BPHParticleMasses::Lambda0_bMass ) <
            //     fabs( nLb2->composite().mass() - BPHParticleMasses::Lambda0_bMass ) )
            //    pxt = nLb1;
            //else 
            //    pxt = nLb2;
            pxt = nLb1;


            if ( !   massSel->accept( *pxt ) ) continue;
            //if ( !massTmpSel->accept( *pxt ) ) continue; // asdf disableed when I don't need to change ptk and ntk
            if ( !   chi2Sel->accept( *pxt ) ) continue;
            if ( massConstr ) if ( ! mFitSel->accept( *pxt ) ) continue;
            tmpList.push_back(pxt);
        }
    }


    updated = true;
    return tmpList;
}

/// set cuts
void BPHLbToJPsiTkTkBuilder::setJPsiMassMin( double m ) {
  updated = false;
  jpsiSel->setMassMin( m );
  return;
}


void BPHLbToJPsiTkTkBuilder::setJPsiMassMax( double m ) {
  updated = false;
  jpsiSel->setMassMax( m );
  return;
}


void BPHLbToJPsiTkTkBuilder::setPtMin      ( double m ) {
    updated = false;
    ptSel->setPtMin( m );
    return;
}

void BPHLbToJPsiTkTkBuilder::setEtaMax     ( double m ) {
    updated = false;
    etaSel->setEtaMax( m );
    return;
}

void BPHLbToJPsiTkTkBuilder::setMassMin( double m ) {
  updated = false;
  massSel->setMassMin( m );
  return;
}


void BPHLbToJPsiTkTkBuilder::setMassMax( double m ) {
  updated = false;
  massSel->setMassMax( m );
  return;
}


void BPHLbToJPsiTkTkBuilder::setProbMin( double p ) {
  updated = false;
  chi2Sel->setProbMin( p );
  return;
}

void BPHLbToJPsiTkTkBuilder::setMassFitMin( double m )
{
    updated = false;
    mFitSel->setMassMin( m );
    return;
}

void BPHLbToJPsiTkTkBuilder::setMassFitMax( double m )
{
    updated = false;
    mFitSel->setMassMax( m );
    return;
}


void BPHLbToJPsiTkTkBuilder::setConstr( double mass, double sigma ) {
  updated = false;
  massConstr = true;
  mFitSel->setFitConstraint( mFitSel->getConstrainedName(), mass, sigma );
  return;
}

void BPHLbToJPsiTkTkBuilder::setConstr( bool m ) {
  updated = false;
  massConstr = true;
  return;
}



/// get current cuts
double BPHLbToJPsiTkTkBuilder::getJPsiMassMin() const {
  return jpsiSel->getMassMin();
}


double BPHLbToJPsiTkTkBuilder::getJPsiMassMax() const {
  return jpsiSel->getMassMax();
}


double BPHLbToJPsiTkTkBuilder::getPtMin() const {
  return ptSel->getPtMin();
}


double BPHLbToJPsiTkTkBuilder::getEtaMax() const {
  return etaSel->getEtaMax();
}


double BPHLbToJPsiTkTkBuilder::getMassMin() const {
  return massSel->getMassMin();
}


double BPHLbToJPsiTkTkBuilder::getMassMax() const {
  return massSel->getMassMax();
}


double BPHLbToJPsiTkTkBuilder::getProbMin() const {
  return chi2Sel->getProbMin();
}


double BPHLbToJPsiTkTkBuilder::getMassFitMin() const {
    return mFitSel->getMassMin();
}

double BPHLbToJPsiTkTkBuilder::getMassFitMax() const {
    return mFitSel->getMassMax();
}

double BPHLbToJPsiTkTkBuilder::getMassFitMass() const {
  return mFitSel->getMass();
}


double BPHLbToJPsiTkTkBuilder::getMassFitSigma() const {
  return mFitSel->getSigma();
}


bool BPHLbToJPsiTkTkBuilder::getConstr() const {
  return massConstr;
}

