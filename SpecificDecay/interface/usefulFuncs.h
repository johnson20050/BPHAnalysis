#ifndef __usefulFuncs_h__
#define __usefulFuncs_h__
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"

namespace usefulFuncs
{
    // check two pat::CompositeCandidate are the same one or not.
    // you need to loop this function to reach the target.

    // algorithm:
    // use pt preselection and find minimize delta R.
    //          ** delta R is defined as (delta phi)^2+(delta eta)^2
    // if delta R is smaller than input delta R, the program returns true.
    bool ccMatch( const pat::CompositeCandidate& p1, const pat::CompositeCandidate& p2, double& ref_minDeltaR2 );
    bool candidateMatch( const RefCountedKinematicParticle& refCand, const reco::Candidate* cand, double& _minDeltaR2 );

    struct recoParticleInfo // {{{
    {
        // this struct uses the map on particle name and reco::Candidate
        // and finally pair particle name and RefCountedKinematicParticle
        recoParticleInfo( const std::string& compName, const std::string& daugName, const reco::Candidate* inptr );
        std::string getFullName() const;
        std::string getCompName() const;
        std::string getDaugName() const;
        int getRecoCharge() const;
        float getRecoPt() const;
        const reco::Candidate* getRecoParticle() const;
        GlobalVector getRefitMom() const;
        int getRefitCharge() const;
        void setRefitParticle( RefCountedKinematicParticle& obj );
        RefCountedKinematicParticle getRefitParticle() const;
        bool isCompDaughter() const;
        bool checkSameCharge() const;

        // used to get jpsi and pv
        template<class T>
        static const T* getByRef( const pat::CompositeCandidate& cand, const std::string& name )
        {
            if ( cand.hasUserData( name ) )
            {
                typedef edm::Ref< std::vector<T> > objRef;
                const objRef* ref = cand.userData<objRef>( name );
                if ( ref ==      0 ) return nullptr;
                if ( ref->isNull() ) return nullptr;
                return ref->get();
            }
            return nullptr;
        }
        //private:
        std::string cName;
        std::string dName;
        const reco::Candidate* cptr;
        RefCountedKinematicParticle refptr;
    }; // end struct }}}

    reco::TrackRef getTrackRefFromRC( const reco::Candidate& rc );
    reco::TrackRef getTrackRefFromPF( const reco::Candidate& rc );
    reco::TrackRef getTrackRefFromGP( const reco::Candidate& rc );

}
#endif
