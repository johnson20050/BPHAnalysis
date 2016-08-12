#ifndef BPHParticleEtaSelect_H
#define BPHParticleEtaSelect_H
/** \class BPHParticleEtaSelect
 *
 *  Descrietaion: 
 *     Class for particle selection by eta
 *
 *
 *  $Date: 2016-05-03 14:57:34 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoSelect.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------


//---------------
// C++ Headers --
//---------------


//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHParticleEtaSelect: public BPHRecoSelect {

 public:

  /** Constructor
   */
  BPHParticleEtaSelect( double eta );

  /** Destructor
   */
  virtual ~BPHParticleEtaSelect();

  /** Operations
   */
  /// select particle
  virtual bool accept( const reco::Candidate& cand ) const;

  /// set eta max
  void setEtaMax( double eta );

  /// get current eta max
  double getEtaMax() const;

 private:

  // private copy and assigment constructors
  BPHParticleEtaSelect           ( const BPHParticleEtaSelect& x );
  BPHParticleEtaSelect& operator=( const BPHParticleEtaSelect& x );

  double etaMax;

};


#endif // BPHParticleEtaSelect_H

