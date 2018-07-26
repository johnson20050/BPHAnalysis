#ifndef BPHAnalysis_SpecificDecay_BPHCompChargeSelect_h
#define BPHAnalysis_SpecificDecay_BPHCompChargeSelect_h
/** \class BPHCompChargeSelect
 *
 *  Description: 
 *  put a limitaion on final composite particle. Force the particle need to take charge
 *
 *  \author : Lian-Sheng, Tsai. 2017/02/14
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHMomentumSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassCuts.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHDecayMomentum.h"

//---------------
// C++ Headers --
//---------------


//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHCompChargeSelect: public BPHMomentumSelect {

 public:

  /** Constructor
   */
  BPHCompChargeSelect( int c ): charge(c) {}

  /** Destructor
   */
  virtual ~BPHCompChargeSelect() {}

  /** Operations
   */
  /// select particle
  virtual bool accept( const BPHDecayMomentum& cand ) const {
      return ( cand.composite().charge() == charge );
  }

  int getCharge() const 
  { return charge; }

  void setCharge( int c )
  { charge = c; return; }



 private:

  // private copy and assigment constructors
  BPHCompChargeSelect           ( const BPHCompChargeSelect& x );
  BPHCompChargeSelect& operator=( const BPHCompChargeSelect& x );
  int charge;

};


#endif

