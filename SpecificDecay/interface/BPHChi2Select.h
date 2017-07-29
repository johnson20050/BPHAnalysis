#ifndef BPHAnalysis_SpecificDecay_BPHChi2Select_h
#define BPHAnalysis_SpecificDecay_BPHChi2Select_h
/** \class BPHChi2Select
 *
 *  Description: 
 *     Class for candidate selection by chisquare (at vertex fit level)
 *
 *  \author Paolo Ronchese INFN Padova
 *
 */

//----------------------
// Base Class Headers --
//----------------------
#include "BPHAnalysis/RecoDecay/interface/BPHVertexSelect.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
#include "BPHAnalysis/RecoDecay/interface/BPHDecayVertex.h"
#include "TMath.h"
#include <vector>

//--------------------------------
// Jack's code used in accept() --
//--------------------------------
#include "RecoVertex/KinematicFitPrimitives/interface/VirtualKinematicParticleFactory.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"


//              ---------------------
//              -- Class Interface --
//              ---------------------

class BPHChi2Select: public BPHVertexSelect {

 public:

  /** Constructor
   */
  BPHChi2Select( double prob ): probMin( prob ) {}

  /** Destructor
   */
  virtual ~BPHChi2Select() {}

  /** Operations
   */
  /// select vertex
  virtual bool accept( const BPHDecayVertex& cand ) const {
    const reco::Vertex& v = cand.vertex();
    if ( v.isFake() ) return false;
    if ( !v.isValid() ) return false;
    if ( probMin == -1.0 ) return true;
    return ( TMath::Prob( v.chi2(), lround( v.ndof() ) ) > probMin );
  }

  // from Jack's code. added in Dec 1 2016
  // Specifically used in lambda0, because it has defferent vertex point.
  // e.g. Lambda0_b -> J/psi + lambda0
  // J/psi has very short lifetime, you can regard the vertex point of J/psi and Lambda0_b as the same.
  // But lambda0 has long life time(about 0.4cm). Thus you need to assign for different vertex point and recalculate the momentum
  virtual bool accept( const BPHRecoCandidate* cand ) const {
    // perform vertex fit
    std::vector<RefCountedKinematicParticle> kComp = cand->kinParticles();

    KinematicParticleVertexFitter vtxFitter;
    RefCountedKinematicTree compTree = vtxFitter.fit( kComp );
    if (compTree->isEmpty()) return false;
    compTree->movePointerToTheTop();
    const RefCountedKinematicParticle kPart = compTree->currentParticle();
    const RefCountedKinematicVertex   kVtx  = compTree->currentDecayVertex();
    const KinematicState& kState = kPart->currentState();
    if (!kState.isValid()) return false;

    if ( probMin == -1.0 ) return true;
    double vtxProb = TMath::Prob( kVtx->chiSquared(), kVtx->degreesOfFreedom() );
    return ( vtxProb > probMin );
  }

  /// set prob min
  void setProbMin( double p ) { probMin = p; return; }

  /// get current prob min
  double getProbMin() const { return probMin; }

 private:

  // private copy and assigment constructors
  BPHChi2Select           ( const BPHChi2Select& x );
  BPHChi2Select& operator=( const BPHChi2Select& x );

  double probMin;

};


#endif

