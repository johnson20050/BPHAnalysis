#ifndef BPHLbToJPsiTkTkBuilder_H
#define BPHLbToJPsiTkTkBuilder_H
/** \class BPHLbToJPsiTkTkBuilder
 *
 *  Description:
 *     Class to build Lb to Jpsi 2 tracks candidates
 *     assume one track is proton track, and the other is kaon track
 *
 *
 *  $Date: 2015-07-24 11:29:20 $
 *  $Revision: 1.1 $
 *  \author Paolo Ronchese INFN Padova
 *
 *  Update: Modifier Lian-Sheng, Tsai
 *          Mon Feb 13 21:27:07 CST 2017
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

class BPHParticlePtSelect;
class BPHParticleEtaSelect;
class BPHMassSelect;
class BPHMassSymSelect_jpsi;
class BPHChi2Select;
class BPHMassFitSelect;
class BPHParticleNeutralVeto;
class BPHCompChargeSelect;

//---------------
// C++ Headers --
//---------------
#include <string>
#include <vector>

//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHLbToJPsiTkTkBuilder
{

public:

    /** Constructor
     */
    BPHLbToJPsiTkTkBuilder( const edm::EventSetup& es,
                            const std::vector<BPHPlusMinusConstCandPtr>& jPsiCollection,
                            const std::vector<BPHPlusMinusConstCandPtr>& TkTkCollection,
                            const std::string ptkName, const double ptkMass, const double ptkSigma,
                            const std::string ntkName, const double ntkMass, const double ntkSigma);

    /** Destructor
     */
    virtual ~BPHLbToJPsiTkTkBuilder();

    /** Operations
     */
    /// build Bs candidates
    std::vector<BPHRecoConstCandPtr> build();

    /// set cuts
    void setJPsiMassMin ( double m  );
    void setJPsiMassMax ( double m  );
    void setPtMin       ( double m  );
    void setEtaMax      ( double m  );
    void setMassMin     ( double m  );
    void setMassMax     ( double m  );
    void setProbMin     ( double p  );
    void setMassFitMin  ( double m  );
    void setMassFitMax  ( double m  );
    void setConstr      ( double mass, double sigma );
    void setConstr      ( bool   m  );

    /// get current cuts
    double getJPsiMassMin () const;
    double getJPsiMassMax () const;
    double getPtMin       () const;
    double getEtaMax      () const;
    double getMassMin     () const;
    double getMassMax     () const;
    double getProbMin     () const;
    double getMassFitMin  () const;
    double getMassFitMax  () const;
    double getMassFitMass () const;
    double getMassFitSigma() const;
    bool   getConstr      () const;

private:

    // private copy and assigment constructors
    BPHLbToJPsiTkTkBuilder           ( const BPHLbToJPsiTkTkBuilder& x );
    BPHLbToJPsiTkTkBuilder& operator=( const BPHLbToJPsiTkTkBuilder& x );


    const edm::EventSetup* evSetup;
    const std::vector<BPHPlusMinusConstCandPtr>*  jpsiCollection;
    const std::vector<BPHPlusMinusConstCandPtr>*  tktkCollection;
    std::string jPsiName;
    const std::string pName;
    double pMass;
    double pSigma;
    std::string nName;
    double nMass;
    double nSigma;

    BPHMassSelect           *     jpsiSel;
    BPHParticlePtSelect     *       ptSel;
    BPHParticleEtaSelect    *      etaSel;
    BPHParticleNeutralVeto  *    nVetoSel;

    BPHMassSelect           *     massSel;
    BPHMassSymSelect_jpsi   *  massTmpSel;
    BPHChi2Select           *     chi2Sel;
    BPHMassFitSelect        *     mFitSel;
    BPHCompChargeSelect     *   chargeSel;


    bool massConstr;
    double cMass;
    double cSigma;
    float minPDiff;

    // if you have any change on the particle, the 'updated' will become false
    // Then when you build(), you can update the particle
    bool updated;

    std::vector<BPHRecoConstCandPtr> lbList;
  TwoTrackMassKinematicConstraint* jpsiConstr;

};


#endif // BPHLbToJPsiTkTkBuilder_H

