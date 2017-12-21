#include "BPHAnalysis/SpecificDecay/interface/usefulFuncs.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "TVector2.h"
#include <vector>
#include <utility> // std::pair

bool usefulFuncs::ccMatch( const pat::CompositeCandidate& p1, const pat::CompositeCandidate& p2, double& ref_minDeltaR2 )
{
    double deltaPT = fabs( p1.pt() - p2.pt() );
    if ( deltaPT/p1.pt() > 0.8 ) return false;
    double deltaR2 = pow( TVector2::Phi_mpi_pi( p1.phi() - p2.phi() ), 2 ) + pow( p1.eta() - p2.eta(), 2 );
    if ( ref_minDeltaR2 < deltaR2 ) return false;
    ref_minDeltaR2 = deltaR2;
    return true;
}
bool usefulFuncs::candidateMatch( const RefCountedKinematicParticle& refCand, const reco::Candidate* cand, double& ref_minDeltaR2 )
{
    const GlobalVector& refMom = refCand->currentState().kinematicParameters().momentum();
    double deltaPT = fabs( refMom.transverse() - cand->pt() );
    if ( deltaPT/refMom.transverse() > 0.8 ) return false;
    double deltaR2 = pow( TVector2::Phi_mpi_pi( refMom.phi() - cand->phi() ), 2 ) + pow( refMom.eta() - cand->eta(), 2 );
    if ( ref_minDeltaR2 < deltaR2 ) return false;
    ref_minDeltaR2 = deltaR2;
    return true;
}
usefulFuncs::recoParticleInfo::recoParticleInfo( const std::string& compName, const std::string& daugName, const reco::Candidate* inptr ) :
    cName( compName ),
    dName( daugName ),
    cptr( inptr ),
    refptr( nullptr ) {}
std::string usefulFuncs::recoParticleInfo::getFullName() const
{
    if ( this->isCompDaughter() )
        return cName+"/"+dName;
    return dName;
}
std::string usefulFuncs::recoParticleInfo::getCompName() const
{ return cName; }
std::string usefulFuncs::recoParticleInfo::getDaugName() const
{ return dName; }
int usefulFuncs::recoParticleInfo::getRecoCharge() const
{ return cptr->charge(); }
float usefulFuncs::recoParticleInfo::getRecoPt() const
{ return cptr->pt(); }
const reco::Candidate* usefulFuncs::recoParticleInfo::getRecoParticle() const
{ return cptr; }
GlobalVector usefulFuncs::recoParticleInfo::getRefitMom() const
{ return refptr->currentState().kinematicParameters().momentum(); }
int usefulFuncs::recoParticleInfo::getRefitCharge() const
{ return refptr->currentState().particleCharge(); }
void usefulFuncs::recoParticleInfo::setRefitParticle( RefCountedKinematicParticle& obj )
{ refptr.swap(obj); }
RefCountedKinematicParticle usefulFuncs::recoParticleInfo::getRefitParticle() const
{ return refptr; }
bool usefulFuncs::recoParticleInfo::isCompDaughter() const
{ return !cName.empty(); }
bool usefulFuncs::recoParticleInfo::checkSameCharge() const
{
    if ( cptr && refptr )
    {
        if ( this->getRecoCharge() == this->getRefitCharge() )
            return true;
    }
    else
        printf ("recoParticleInfo::cptr or refptr doesn't exit!\n");
    return false;
}
