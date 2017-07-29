#ifndef BPHLambda0_bToPentaQKBuilder_H
#define BPHLambda0_bToPentaQKBuilder_H
/** \class BPHLambda0_bToPentaQKBuilder
 *
 *  Description: 
 *     Class to build Lambda0_b to penta_quark K+- candidates
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

class BPHParticleNeutralVeto;
class BPHParticlePtSelect;
class BPHParticleEtaSelect;
class BPHMassSelect;
class BPHChi2Select;
class BPHMassFitSelect;

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHLambda0_bToPentaQKBuilder {

 public:

  /** Constructor
   */
  BPHLambda0_bToPentaQKBuilder( const edm::EventSetup& es,
      const std::vector<BPHPlusMinusConstCandPtr>& PentaQCollection,
      const BPHRecoBuilder::BPHGenericCollection*    kaonCollection ); 

  /** Destructor
   */
  virtual ~BPHLambda0_bToPentaQKBuilder();

  /** Operations
   */
  /// build Bu candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// set cuts
  void setKPtMin       ( double pt  );
  void setKEtaMax      ( double eta );
  void setPenQMassMin  ( double m   );
  void setPenQMassMax  ( double m   );
  void setMassMin      ( double m   );
  void setMassMax      ( double m   );
  void setProbMin      ( double p   );
  void setMassFitMin   ( double m   );
  void setMassFitMax   ( double m   );
  void setConstr       ( bool flag  );

  /// get current cuts
  double getKPtMin       () const;
  double getKEtaMax      () const;
  double getPenQMassMin  () const;
  double getPenQMassMax  () const;
  double getMassMin      () const;
  double getMassMax      () const;
  double getProbMin      () const;
  double getMassFitMin   () const;
  double getMassFitMax   () const;
  bool   getConstr       () const;

 private:

  // private copy and assigment constructors
  BPHLambda0_bToPentaQKBuilder           ( const BPHLambda0_bToPentaQKBuilder& x );
  BPHLambda0_bToPentaQKBuilder& operator=( const BPHLambda0_bToPentaQKBuilder& x );

  std::string   penQName;
  std::string   kaonName;

  const edm::EventSetup* evSetup;
  const BPHRecoBuilder::BPHGenericCollection*  pCollection;
  const BPHRecoBuilder::BPHGenericCollection*  kCollection;

  BPHParticlePtSelect   *     ptSel;
  BPHParticleEtaSelect  *    etaSel;
  BPHMassSelect         *  mPenQSel;

  BPHMassSelect         *   massSel;
  BPHChi2Select         *   chi2Sel;
  BPHMassFitSelect      *   mFitSel;

  bool massConstr;
  float minPDiff;
  bool updated;

  std::vector<BPHRecoConstCandPtr> lbList;

};


#endif // BPHLambda0_bToPentaQKBuilder_H

