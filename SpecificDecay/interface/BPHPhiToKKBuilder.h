#ifndef BPHPhiToKKBuilder_H
#define BPHPhiToKKBuilder_H
/** \class BPHPhiToKKBuilder
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

class BPHPhiToKKBuilder {

 public:

  /** Constructor
   */
  BPHPhiToKKBuilder( const edm::EventSetup& es,
       const BPHRecoBuilder::BPHGenericCollection* kPosCollection,
       const BPHRecoBuilder::BPHGenericCollection* kNegCollection );

  /** Destructor
   */
  virtual ~BPHPhiToKKBuilder();

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
  BPHPhiToKKBuilder           ( const BPHPhiToKKBuilder& x );
  BPHPhiToKKBuilder& operator=( const BPHPhiToKKBuilder& x );

  std::string kPosName;
  std::string kNegName;

  const edm::EventSetup* evSetup;
  const BPHRecoBuilder::BPHGenericCollection* posCollection;
  const BPHRecoBuilder::BPHGenericCollection* negCollection;

  BPHParticlePtSelect *  ptSel;
  BPHParticleEtaSelect* etaSel;
  BPHMassSelect* massSel;
  BPHChi2Select* chi2Sel;
  double cMass;
  double cSigma;
  bool updated;

  std::vector<BPHPlusMinusConstCandPtr> phiList;

};


#endif // BPHPhiToKKBuilder_H

