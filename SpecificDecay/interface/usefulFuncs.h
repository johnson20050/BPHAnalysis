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
// my added struct recoParticleInfo {{{
// this struct uses the map on particle name and reco::Candidate
// and finally pair particle name and RefCountedKinematicParticle
    struct recoParticleInfo
    {
        recoParticleInfo( const std::string& compName, const std::string& daugName, const reco::Candidate* inptr );// :
        //    cName( compName ),
        //    dName( daugName ),
        //    cptr( inptr ),
        //    refptr( nullptr ) {}
        std::string getFullName() const;
        //{
        //    if ( this->isCompDaughter() )
        //        return cName+"/"+dName;
        //    return dName;
        //}

        std::string getCompName() const;
        //{ return cName; }
        std::string getDaugName() const;
        //{ return dName; }
        int getRecoCharge() const;
        //{ return cptr->charge(); }
        float getRecoPt() const;
        //{ return cptr->pt(); }
        const reco::Candidate* getRecoParticle() const;
        //{ return cptr; }
        GlobalVector getRefitMom() const;
        //{ return refptr->currentState().kinematicParameters().momentum(); }
        int getRefitCharge() const;
        //{ return refptr->currentState().particleCharge(); }
        void setRefitParticle( RefCountedKinematicParticle& obj );
        //{ refptr.swap(obj); }
        RefCountedKinematicParticle getRefitParticle() const;
        //{ return refptr; }
        bool isCompDaughter() const;
        //{ return !cName.empty(); }
        bool checkSameCharge() const;
        //{
        //    if ( cptr && refptr )
        //    {
        //        if ( this->getRecoCharge() == this->getRefitCharge() )
        //            return true;
        //    }
        //    else
        //        printf ("recoParticleInfo::cptr or refptr doesn't exit!\n");
        //    return false;
        //}

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
    };
// my added recoParticleInfo end}}}
}
#endif
