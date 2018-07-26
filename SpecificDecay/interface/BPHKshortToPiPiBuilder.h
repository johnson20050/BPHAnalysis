#ifndef BPHAnalysis_SpecificDecay_BPHKshortToPiPiBuilder_h
#define BPHAnalysis_SpecificDecay_BPHKshortToPiPiBuilder_h
/** \class BPHKshortToPiPiBuilder
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
class BPHMassFitSelect;

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHKshortToPiPiBuilder {

 public:

  /** Constructor
   */
  BPHKshortToPiPiBuilder( const edm::EventSetup& es,
       const BPHRecoBuilder::BPHGenericCollection* posCollection,
       const BPHRecoBuilder::BPHGenericCollection* negCollection );

  /** Destructor
   */
  virtual ~BPHKshortToPiPiBuilder();

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
  double getMassFitMin () const;
  double getMassFitMax () const;
  double getConstrMass () const;
  double getConstrSigma() const;

 private:

  // private copy and assigment constructors
  BPHKshortToPiPiBuilder           ( const BPHKshortToPiPiBuilder& x );
  BPHKshortToPiPiBuilder& operator=( const BPHKshortToPiPiBuilder& x );

  std::string pName;
  std::string nName;

  const edm::EventSetup* evSetup;
  const BPHRecoBuilder::BPHGenericCollection* pCollection;
  const BPHRecoBuilder::BPHGenericCollection* nCollection;

  BPHParticlePtSelect *  ptSel;
  BPHParticleEtaSelect* etaSel;
  BPHMassSelect* massSel;
  BPHChi2Select* chi2Sel;
  BPHMassFitSelect*  mFitSel;
  double cMass;
  double cSigma;
  bool massConstr;
  bool updated;

  std::vector<BPHPlusMinusConstCandPtr> ksList;

};


#endif

