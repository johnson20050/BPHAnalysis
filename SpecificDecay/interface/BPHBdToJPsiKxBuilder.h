#ifndef BPHBdToJPsiKxBuilder_H
#define BPHBdToJPsiKxBuilder_H
/** \class BPHBdToJPsiKxBuilder
 *
 *  Description: 
 *     Class to build Bs to JPsi Phi candidates
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

class BPHBdToJPsiKxBuilder {

 public:

  /** Constructor
   */
  BPHBdToJPsiKxBuilder( const edm::EventSetup& es,
      const std::vector<BPHPlusMinusConstCandPtr>& jpsiCollection,
      const std::vector<BPHPlusMinusConstCandPtr>&  kx0Collection );

  /** Destructor
   */
  virtual ~BPHBdToJPsiKxBuilder();

  /** Operations
   */
  /// build Bs candidates
  std::vector<BPHRecoConstCandPtr> build();

  /// set cuts
  void setJPsiMassMin( double m  );
  void setJPsiMassMax( double m  );
  void setKxMassMin  ( double m  );
  void setKxMassMax  ( double m  );
  void setMassMin    ( double m  );
  void setMassMax    ( double m  );
  void setProbMin    ( double p  );
  void setMassFitMin ( double m  );
  void setMassFitMax ( double m  );
  void setConstr     ( bool flag );

  /// get current cuts
  double getJPsiMassMin() const;
  double getJPsiMassMax() const;
  double getKxMassMin  () const;
  double getKxMassMax  () const;
  double getMassMin    () const;
  double getMassMax    () const;
  double getProbMin    () const;
  double getMassFitMin () const;
  double getMassFitMax () const;
  bool   getConstr     () const;

 private:

  // private copy and assigment constructors
  BPHBdToJPsiKxBuilder           ( const BPHBdToJPsiKxBuilder& x );
  BPHBdToJPsiKxBuilder& operator=( const BPHBdToJPsiKxBuilder& x );

  std::string jPsiName;
  std::string  kx0Name;

  const edm::EventSetup* evSetup;
  const std::vector<BPHPlusMinusConstCandPtr>* jCollection;
  const std::vector<BPHPlusMinusConstCandPtr>* kCollection;

  BPHMassSelect   * jpsiSel;
  BPHMassSelect   * mkx0Sel;

  BPHMassSelect   * massSel;
  BPHChi2Select   * chi2Sel;
  BPHMassFitSelect* mFitSel;

  bool massConstr;
  float minPDiff;
  bool updated;

  std::vector<BPHRecoConstCandPtr> bdList;

};


#endif // BPHBdToJPsiKxBuilder_H

