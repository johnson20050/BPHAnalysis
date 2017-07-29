#ifndef BPHAnalysis_SpecificDecay_BPHMassSymSelect_jpsi_h
#define BPHAnalysis_SpecificDecay_BPHMassSymSelect_jpsi_h
/** \class BPHMassSymSelect_jpsi
 *
 *  Description: 
 *     Class for candidate selection by invariant mass (at momentum sum level)
 *     allowing for decay product mass swap
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHMomentumSelect.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHDecayMomentum.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"

//---------------
// C++ Headers --
//---------------
#include <string>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHMassSymSelect_jpsi: public BPHMomentumSelect {

 public:

  /** Constructor
   */
  BPHMassSymSelect_jpsi( const std::string& np, const std::string& nn,
                    const BPHMassSelect* ms ): nPos( np ), nNeg( nn ),
                                               mSel( ms ) {}

  /** Destructor
   */
  virtual ~BPHMassSymSelect_jpsi() {}

  /** Operations
   */
  /// select particle
  virtual bool accept( const BPHDecayMomentum& cand ) const {

    if ( mSel->accept( cand ) ) return true;

    const
    reco::Candidate* pp = cand.getDaug( nPos );
    const
    reco::Candidate* np = cand.getDaug( nNeg );

    const 
    BPHRecoCandidate* cand_jpsi = cand.getComp( "JPsi" ).get();

    reco::Candidate* pc = cand.originalReco( pp )->clone();
    reco::Candidate* nc = cand.originalReco( np )->clone();

    pc->setMass( np->p4().mass() );
    nc->setMass( pp->p4().mass() );
    const reco::Candidate::LorentzVector  s4 = pc->p4() + nc->p4() + cand_jpsi->composite().p4();
    double mass = s4.mass();

    delete pc;
    delete nc;
    return ( ( mass > mSel->getMassMin() ) &&
             ( mass < mSel->getMassMax() ) );

  }

 private:

  // private copy and assigment constructors
  BPHMassSymSelect_jpsi           ( const BPHMassSymSelect_jpsi& x );
  BPHMassSymSelect_jpsi& operator=( const BPHMassSymSelect_jpsi& x );

  std::string nPos;
  std::string nNeg;
  const BPHMassSelect* mSel;

};


#endif

