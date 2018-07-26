#ifndef BPHLambda0_bToJPsiLambda0Builder_H
#define BPHLambda0_bToJPsiLambda0Builder_H
/** \class BPHLambda0_bToJPsiLambda0Builder
 *
 *  Description: 
 *     Class to build Lb to Jpsi Lambda* candidates
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
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"

#include "FWCore/Framework/interface/Event.h"

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

class BPHLambda0_bToJPsiLambda0Builder {

 public:

  /** Constructor
   */
  BPHLambda0_bToJPsiLambda0Builder( const edm::EventSetup& es,
      const std::vector<BPHPlusMinusConstCandPtr>& JpsiCollection,
      const std::vector<BPHPlusMinusConstCandPtr>& Lam0Collection );

  /** Destructor
   */
  virtual ~BPHLambda0_bToJPsiLambda0Builder();

  /** Operations
   */
  /// build Bs candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// set cuts
  void setJPsiMassMin( double m  );
  void setJPsiMassMax( double m  );
  void setLam0MassMin( double m  );
  void setLam0MassMax( double m  );
  void setMassMin    ( double m  );
  void setMassMax    ( double m  );
  void setProbMin    ( double p  );
  void setMassFitMin ( double m  );
  void setMassFitMax ( double m  );
  void setConstr     ( bool flag );

  /// get current cuts
  double getJPsiMassMin() const;
  double getJPsiMassMax() const;
  double getLam0MassMin() const;
  double getLam0MassMax() const;
  double getMassMin    () const;
  double getMassMax    () const;
  double getProbMin    () const;
  double getMassFitMin () const;
  double getMassFitMax () const;
  bool   getConstr     () const;

 private:

  // private copy and assigment constructors
  BPHLambda0_bToJPsiLambda0Builder           ( const BPHLambda0_bToJPsiLambda0Builder& x );
  BPHLambda0_bToJPsiLambda0Builder& operator=( const BPHLambda0_bToJPsiLambda0Builder& x );

  std::string jPsiName;
  std::string lam0Name;

  const edm::EventSetup* evSetup;
  const std::vector<BPHPlusMinusConstCandPtr>* jpsiCollection;
  const std::vector<BPHPlusMinusConstCandPtr>* lam0Collection;

  BPHMassSelect   *  jpsiSel;
  BPHMassSelect   * mlam0Sel;

  BPHMassSelect   *  massSel;
  BPHChi2Select   *  chi2Sel;
  BPHMassFitSelect*  mFitSel;

  bool massConstr;
  float minPDiff;

  // if you have any change on the particle, the 'updated' will become false
  // Then when you build(), you can update the particle
  bool updated;

  std::vector<BPHRecoConstCandPtr> lbList;
  TwoTrackMassKinematicConstraint* jpsiConstr;

};


#endif // BPHLambda0_bToJPsiLambda0Builder_H

