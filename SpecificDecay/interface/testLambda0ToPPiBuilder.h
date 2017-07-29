#ifndef testAnalysis_smallAnalyzer_testLambda0ToPPiBuilder_h
#define testAnalysis_smallAnalyzer_testLambda0ToPPiBuilder_h
/** \class testLambda0ToPPiBuilder
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

class testLambda0ToPPiBuilder
{

public:

    /** Constructor
     */
    testLambda0ToPPiBuilder( const edm::EventSetup& es,
                             const BPHRecoBuilder::BPHGenericCollection* protCollection,
                             const BPHRecoBuilder::BPHGenericCollection* pionCollection );

    /** Destructor
     */
    virtual ~testLambda0ToPPiBuilder();

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
    testLambda0ToPPiBuilder           ( const testLambda0ToPPiBuilder& x );
    testLambda0ToPPiBuilder& operator=( const testLambda0ToPPiBuilder& x );

    std::string protName;
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

    std::vector<BPHPlusMinusConstCandPtr> lam0List;

};


#endif

