#ifndef BPHAnalysis_SpecificDecay_BPHTkTkBuilder_h
#define BPHAnalysis_SpecificDecay_BPHTkTkBuilder_h
/** \class BPHTkTkBuilder
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

class BPHTkTkBuilder {

 public:

  /** Constructor
   */
  BPHTkTkBuilder( const edm::EventSetup& es,
       const BPHRecoBuilder::BPHGenericCollection* ptkCollection,
       const BPHRecoBuilder::BPHGenericCollection* mtkCollection,
       const std::string ptkName, const double ptkMass, const double ptkSigma,
       const std::string ntkName, const double ntkMass, const double ntkSigma);

  /** Destructor
   */
  virtual ~BPHTkTkBuilder();

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
  void setpTkCut ( bool flag ) { ptkCut = flag; }
  void setnTkCut ( bool flag ) { ntkCut = flag; }

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
  BPHTkTkBuilder           ( const BPHTkTkBuilder& x );
  BPHTkTkBuilder& operator=( const BPHTkTkBuilder& x );


  const edm::EventSetup* evSetup;
  const BPHRecoBuilder::BPHGenericCollection* pCollection;
  const BPHRecoBuilder::BPHGenericCollection* nCollection;
  const std::string pName;
  const double pMass;
  const double pSigma;
  const std::string nName;
  const double nMass;
  const double nSigma;

  BPHParticlePtSelect * ptSel;
  BPHParticleEtaSelect* etaSel;
  BPHMassSelect* massSel;
  BPHChi2Select* chi2Sel;
  double cMass;
  double cSigma;
  bool updated;
  bool ptkCut, ntkCut;

  std::vector<BPHPlusMinusConstCandPtr> tktkList;

};


#endif

