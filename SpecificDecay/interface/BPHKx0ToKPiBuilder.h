#ifndef BPHKx0ToKPiBuilder_H
#define BPHKx0ToKPiBuilder_H
/** \class BPHKx0ToKPiBuilder
 *
 *  Description: 
 *     Class to build Phi candidates
 *
 *
 *  $Date: 2015-07-24 11:29:20 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------


//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"

#include "FWCore/Framework/interface/Event.h"

class BPHParticlePtSelect;
class BPHParticleEtaSelect;
class BPHChi2Select;
class BPHMassSelect;

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHKx0ToKPiBuilder {

 public:

  /** Constructor
   */
  BPHKx0ToKPiBuilder( const edm::EventSetup& es,
       const BPHRecoBuilder::BPHGenericCollection* kaonCollection,
       const BPHRecoBuilder::BPHGenericCollection* pionCollection );

  /** Destructor
   */
  virtual ~BPHKx0ToKPiBuilder();

  /** Operations
   */
  /// build Phi candidates
  std::vector<BPHPlusMinusConstCandPtr> build();

  /// set cuts
  void setPtMin  ( double pt  );
  void setEtaMax ( double eta );
  void setMassMin( double m   );
  void setMassMax( double m   );
  void setProbMin( double p   );
  void setConstr ( double mass, double sigma );

  /// get current cuts
  double getPtMin  () const;
  double getEtaMax () const;
  double getMassMin() const;
  double getMassMax() const;
  double getProbMin() const;
  double getConstrMass () const;
  double getConstrSigma() const;

 private:

  // private copy and assigment constructors
  BPHKx0ToKPiBuilder           ( const BPHKx0ToKPiBuilder& x );
  BPHKx0ToKPiBuilder& operator=( const BPHKx0ToKPiBuilder& x );

  std::string kaonName;
  std::string pionName;

  const edm::EventSetup* evSetup;
  const BPHRecoBuilder::BPHGenericCollection* kCollection;
  const BPHRecoBuilder::BPHGenericCollection* pCollection;

  BPHParticlePtSelect *  ptSel;
  BPHParticleEtaSelect* etaSel;
  BPHMassSelect* massSel;
  BPHChi2Select* chi2Sel;
  double cMass;
  double cSigma;
  bool updated;

  std::vector<BPHPlusMinusConstCandPtr> kx0List;

};


#endif // BPHKx0ToKPiBuilder_H

