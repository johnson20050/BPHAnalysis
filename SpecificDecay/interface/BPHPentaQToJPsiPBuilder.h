#ifndef BPHPentaQToJPsiPBuilder_H
#define BPHPentaQToJPsiPBuilder_H
/** \class BPHPentaQToJPsiPBuilder
 *
 *  Description: 
 *     Class to build B+- to JPsi K+- candidates
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

class BPHPentaQToJPsiPBuilder {

 public:

  /** Constructor
   */
  BPHPentaQToJPsiPBuilder( const edm::EventSetup& es,
      const std::vector<BPHPlusMinusConstCandPtr>&   JpsiCollection,
      const BPHRecoBuilder::BPHGenericCollection*  ProtonCollection );

  /** Destructor
   */
  virtual ~BPHPentaQToJPsiPBuilder();

  /** Operations
   */
  /// build Bu candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// set cuts
  void setPPtMin     ( double pt  );
  void setPEtaMax    ( double eta );
  void setJPsiMassMin( double m   );
  void setJPsiMassMax( double m   );
  void setMassMin    ( double m   );
  void setMassMax    ( double m   );
  void setProbMin    ( double p   );
  void setMassFitMin ( double m   );
  void setMassFitMax ( double m   );
  void setConstr     ( bool flag  );

  /// get current cuts
  double getPPtMin     () const;
  double getPEtaMax    () const;
  double getJPsiMassMin() const;
  double getJPsiMassMax() const;
  double getMassMin    () const;
  double getMassMax    () const;
  double getProbMin    () const;
  double getMassFitMin () const;
  double getMassFitMax () const;
  bool   getConstr     () const;

 private:

  // private copy and assigment constructors
  BPHPentaQToJPsiPBuilder           ( const BPHPentaQToJPsiPBuilder& x );
  BPHPentaQToJPsiPBuilder& operator=( const BPHPentaQToJPsiPBuilder& x );

  std::string   jPsiName;
  std::string protonName;

  const edm::EventSetup* evSetup;
  const std::vector<BPHPlusMinusConstCandPtr>* jCollection;
  const BPHRecoBuilder::BPHGenericCollection*  pCollection;

  BPHMassSelect         * jpsiSel;
//BPHParticleNeutralVeto*  knVeto;
  BPHParticlePtSelect   *   ptSel;
  BPHParticleEtaSelect  *  etaSel;

  BPHMassSelect         * massSel;
  BPHChi2Select         * chi2Sel;
  BPHMassFitSelect      * mFitSel;

  bool  massConstr;
  float minPDiff;
  bool  updated;

  std::vector<BPHRecoConstCandPtr> pentaQList;

};


#endif // BPHPentaQToJPsiPBuilder_H

