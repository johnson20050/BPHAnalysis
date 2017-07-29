#ifndef BPHAnalysis_SpecificDecay_BPHLambda0ToPPiBuilder_h
#define BPHAnalysis_SpecificDecay_BPHLambda0ToPPiBuilder_h
/** \class BPHLambda0ToPPiBuilder
 *
 *  Description: 
 *     Class to build K*0 to K+ pi- candidates
 *
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

class BPHLambda0ToPPiBuilder {

 public:

  /** Constructor
   */
  BPHLambda0ToPPiBuilder( const edm::EventSetup& es,
       const BPHRecoBuilder::BPHGenericCollection* kaonCollection,
       const BPHRecoBuilder::BPHGenericCollection* pionCollection );

  /** Destructor
   */
  virtual ~BPHLambda0ToPPiBuilder();

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
  BPHLambda0ToPPiBuilder           ( const BPHLambda0ToPPiBuilder& x );
  BPHLambda0ToPPiBuilder& operator=( const BPHLambda0ToPPiBuilder& x );

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


#endif

