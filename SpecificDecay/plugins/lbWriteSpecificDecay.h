#ifndef BPHAnalysis_SpecificDecay_lbWriteSpecificDecay_h
#define BPHAnalysis_SpecificDecay_lbWriteSpecificDecay_h

#include "BPHAnalysis/RecoDecay/interface/BPHAnalyzerTokenWrapper.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"

#include <string>
#include <vector>
#include <map>
#include <map>
#include <iostream>
#include <fstream>

class TH1F;
class BPHRecoCandidate;

class lbWriteSpecificDecay:
    public BPHAnalyzerWrapper<BPHModuleWrapper::one_producer>
{

public:

    explicit lbWriteSpecificDecay( const edm::ParameterSet& ps );
    virtual ~lbWriteSpecificDecay();

    static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );

    virtual void beginJob();
    virtual void produce( edm::Event& ev, const edm::EventSetup& es );
    virtual void fill( edm::Event& ev, const edm::EventSetup& es );
    virtual void endJob();

private:

    std::string pVertexLabel;
    std::string patMuonLabel;
    std::string ccCandsLabel;
    std::string pfCandsLabel;
    std::string pcCandsLabel;
    std::string gpCandsLabel;

    // token wrappers to allow running both on "old" and "new" CMSSW versions
    BPHTokenWrapper< std::vector<reco::Vertex>                 > pVertexToken;
    BPHTokenWrapper< pat::MuonCollection                       > patMuonToken;
    BPHTokenWrapper< std::vector<pat::CompositeCandidate     > > ccCandsToken;
    BPHTokenWrapper< std::vector<reco::PFCandidate           > > pfCandsToken;
    BPHTokenWrapper< std::vector<BPHTrackReference::candidate> > pcCandsToken;
    BPHTokenWrapper< std::vector<pat::GenericParticle        > > gpCandsToken;


    bool usePV;
    bool usePM;
    bool useCC;
    bool usePF;
    bool usePC;
    bool useGP;

// The label used in output product
    std::string     oniaName;
    std::string     Lam0Name;
    std::string     TkTkName;
    std::string LbToLam0Name;
    std::string LbToTkTkName;

    std::string PhiName;
    std::string BsName;

    enum recoType { Onia, Psi1, Psi2, Lam0, TkTk, LbToLam0, LbToTkTk, Phi, Bs };
    enum  parType { ptMin, etaMax,
                    mPsiMin, mPsiMax, mLam0Min, mLam0Max,
                    mPhiMin, mPhiMax,
                    massMin, massMax, probMin, mFitMin, mFitMax,
                    constrMass, constrSigma, constrMJPsi, writeCandidate,
                  };

    std::map<std::string,recoType> rMap;
    std::map<std::string, parType> pMap;
    std::map<std::string, parType> fMap;
    std::map< recoType, std::map<parType,double> >  parMap;


    bool recoOnia     ;
    bool recoLam0     ;
    bool recoTkTk     ;
    bool recoLbToLam0 ;
    bool recoLbToTkTk ;
    bool recoPhi;
    bool recoBs;

    bool writeOnia    ;
    bool writeLam0    ;
    bool writeTkTk    ;
    bool writeLbToLam0;
    bool writeLbToTkTk;
    bool writePhi;
    bool writeBs;

    bool writeVertex;
    bool writeMomentum;

    std::vector<BPHPlusMinusConstCandPtr> lFull;
    std::vector<BPHPlusMinusConstCandPtr> lJPsi;
    std::vector<BPHPlusMinusConstCandPtr> lLam0;
    std::vector<BPHPlusMinusConstCandPtr> lTkTk;
    std::vector<BPHRecoConstCandPtr>      lLbToLam0;
    std::vector<BPHRecoConstCandPtr>      lLbToTkTk;
    std::vector<BPHPlusMinusConstCandPtr> lPhi;
    std::vector<BPHRecoConstCandPtr>      lBs;

    std::map<const BPHRecoCandidate*,const BPHRecoCandidate*> jPsiOMap;
    typedef edm::Ref< std::vector<reco::Vertex> > vertex_ref;
    std::map<const BPHRecoCandidate*,vertex_ref> pvRefMap;
    typedef edm::Ref< pat::CompositeCandidateCollection > compcc_ref;
    std::map<const BPHRecoCandidate*,compcc_ref> ccRefMap;

    void setRecoParameters( const edm::ParameterSet& ps );
// my added struct recoParticleInfo {{{
// this struct uses the map on particle name and reco::Candidate
// and finally pair particle name and RefCountedKinematicParticle
struct recoParticleInfo
{
    recoParticleInfo( const std::string& compName, const std::string& daugName, const reco::Candidate* inptr ) :
        cName( compName ),
        dName( daugName ),
        cptr( inptr ),
        refptr( nullptr ) {}
    std::string getFullName() const
    { 
        if ( this->isCompDaughter() )
            return cName+"/"+dName; 
        return dName;
    }
    std::string getCompName() const
    { return cName; }
    std::string getDaugName() const
    { return dName; }
    int getRecoCharge() const
    { return cptr->charge(); }
    float getRecoPt() const
    { return cptr->pt(); }
    GlobalVector getRefMom() const
    { return refptr->currentState().kinematicParameters().momentum(); }
    int getRefCharge() const
    { return refptr->currentState().particleCharge(); }
    void setRefParticle( RefCountedKinematicParticle& obj )
    { refptr.swap(obj); }
    RefCountedKinematicParticle getRefParticle() const
    { return refptr; }
    bool isCompDaughter() const
    { return !cName.empty(); }
    bool checkSameCharge() const
    {
        if ( cptr && refptr ) 
        {
            if ( this->getRecoCharge() == this->getRefCharge() )
                return true;
        }
        else
            printf ("recoParticleInfo::cptr or refptr doesn't exit!\n");
        return false;
    }
    private:        
    std::string cName;
    std::string dName;
    const reco::Candidate* cptr;
    RefCountedKinematicParticle refptr;
};
// my added recoParticleInfo end}}}
    template <class T>
    edm::OrphanHandle<pat::CompositeCandidateCollection> write( edm::Event& ev,
            const std::vector<T>& list, const std::string& name, bool writeSecParticleInf )
    {
        pat::CompositeCandidateCollection* ccList =
            new pat::CompositeCandidateCollection;
        int i;
        int n = list.size();
        std::map<const BPHRecoCandidate*,
            const BPHRecoCandidate*>::const_iterator jpoIter;
        std::map<const BPHRecoCandidate*,
            const BPHRecoCandidate*>::const_iterator jpoIend = jPsiOMap.end();
        std::map<const BPHRecoCandidate*,vertex_ref>::const_iterator pvrIter;
        std::map<const BPHRecoCandidate*,vertex_ref>::const_iterator pvrIend =
            pvRefMap.end();
        std::map<const BPHRecoCandidate*,compcc_ref>::const_iterator ccrIter;
        std::map<const BPHRecoCandidate*,compcc_ref>::const_iterator ccrIend =
            ccRefMap.end();

        for ( i = 0; i < n; ++i )
        {
            const T& ptr = list[i];
            ccList->push_back( ptr->composite() );
            pat::CompositeCandidate& cc = ccList->back();
            if ( ( pvrIter = pvRefMap.find( ptr.get() ) ) != pvrIend )
                cc.addUserData ( "primaryVertex", pvrIter->second );
            const std::vector<std::string>& cNames = ptr->compNames();
            int j = 0;
            int m = cNames.size();
            while ( j < m )
            {
                const std::string& compName = cNames[j++];
                const BPHRecoCandidate* cptr = ptr->getComp( compName ).get();
                if ( ( ccrIter = ccRefMap.find( cptr ) ) == ccrIend )
                {
                    if ( ( jpoIter = jPsiOMap.find( cptr ) ) != jpoIend )
                        cptr = jpoIter->second;
                    else cptr = 0;
                }
                if ( ( ccrIter = ccRefMap.find( cptr ) ) != ccrIend )
                {
                    compcc_ref cref = ccrIter->second;
                    if ( cref.isNonnull() ) cc.addUserData ( "refTo" + compName, cref );
                }
            }
            const BPHPlusMinusCandidate* pmp =
                dynamic_cast<const BPHPlusMinusCandidate*>( ptr.get() );
            if ( pmp != 0 ) cc.addUserData( "cowboy", pmp->isCowboy() );
            if ( ptr->isEmpty() )
            {
                if ( writeVertex ) cc.addUserData( "vertex" , ptr->vertex() );
                continue;
            }
            if ( writeVertex ) cc.addUserData( "fitVertex",
                                                   reco::Vertex( *ptr->currentDecayVertex() ) );

            // store refit information {{{
            if ( ptr->isValidFit() ) 
            {
                const RefCountedKinematicParticle kinPart = ptr->currentParticle();
                const           KinematicState    kinStat = kinPart->currentState();
                cc.addUserFloat( "fitMass", kinStat.mass() );
                if ( writeMomentum )
                    cc.addUserData ( "fitMomentum",
                                     kinStat.kinematicParameters().momentum() );


                if ( writeSecParticleInf )
                {
                    std::vector<RefCountedKinematicParticle> daughters = ptr->kinematicTree()->finalStateParticles();
                    std::vector<RefCountedKinematicParticle> compDaugh;
                    const std::string& constrName = ptr->getConstrName();
                    std::vector<recoParticleInfo> myParticleList;
                
                    // find RefCountedKinematicParticle
                    if ( ptr->seckinematicTree()->isValid() )
                    {
                        // the last one is the constrained particle, it is virtual particle, need to be removed.
                        daughters.erase( --daughters.end() );
                        compDaugh = ptr->seckinematicTree()->finalStateParticles();
                    }


                    // all RefCountedKinematicParticle are found, check if it losts particle?
                    if ( (daughters.size()+compDaugh.size()) != ptr->daughFull().size() )
                        printf("-----bugs on RefCountedKinematicParticle in write() -----\n");
                    else
                    {
                        //get the daughterFull name stored in BPHRecoCandidate.
                        const std::vector<std::string>& cNames_ = ptr->compNames();
                        const std::vector<std::string>& dNames_ = ptr->daugNames();
                        std::map<const reco::Candidate*, std::string> componentList;
                        for ( const std::string& cName_ : cNames_ )
                        {
                            const BPHRecoCandidate* cand_ = ptr->getComp( cName_ ).get();
                            const std::vector<std::string>& dNames__ = cand_->daugNames();
                            for ( const std::string& dName__ : dNames__ )
                            {
                                const reco::Candidate* cand__ = cand_->getDaug( dName__ );
                                //componentList.insert( std::make_pair(cand__, (cName_+"/"+dName__)) );
                                myParticleList.emplace_back( cName_, dName__, cand__ );
                            }
                        }
                        for ( const std::string& dName_ : dNames_ )
                        {
                            const reco::Candidate* cand_ = ptr->getDaug( dName_ );
                            //componentList.insert( std::make_pair(cand_, dName_) );
                            myParticleList.emplace_back( "", dName_, cand_ );
                        }
                        // get daughFull name end

                        int totalCharge = 0;
                        std::map<const reco::Candidate*, std::string> componentList_ = componentList;
                        std::map< const RefCountedKinematicParticle, std::pair<const reco::Candidate*, std::string> > particleList;
    
                        //use pt and charge to match the reco::Candidate* and RefCountedKinematicParticle
                        //And the constrained name is used to separate virtual particle and other particle 
                        const std::string& constrName = ptr->getConstrName();
                        for ( recoParticleInfo& _oldParticle : myParticleList )
                        {
                            int charge_ = _oldParticle.getRecoCharge();
                            float ptDiff = 999.0;

                            std::vector<RefCountedKinematicParticle> lists_;
                            if ( _oldParticle.getCompName() == constrName )
                                lists_ = compDaugh;
                            else
                                lists_ = daughters;

                            for ( RefCountedKinematicParticle& _cDau : lists_ )
                            {
                                // use charge and pt to distinguish particle
                                if ( _oldParticle.getRecoCharge() != charge_ ) continue;
                                float dPt_ = _oldParticle.getRecoPt();
                                float rPt_ = _cDau->currentState().kinematicParameters().momentum().transverse();
                                if ( fabs( dPt_-rPt_ ) < ptDiff )
                                {
                                    ptDiff = fabs( dPt_-rPt_ );
                                    _oldParticle.setRefParticle( _cDau );
                                }
                            }
                            totalCharge += _oldParticle.getRecoCharge();
                        }

                        // particle name - RefCountedKinematicParticle pair is prepared.
                        if ( myParticleList.size() != ptr->daughFull().size() ) printf( "-----bug: particleList and daughters don't owns the same size! in write() \n" );
                        else if ( totalCharge != 0 ) printf( "-----bug: the built particle is not neutral! in write()\n" );
                        else
                        {
                            for ( const auto& particleContainer : myParticleList )
                            {
                                cc.addUserData ( particleContainer.getFullName()+".fitMom",
                                                 particleContainer.getRefMom() );
                                cc.addUserData ( particleContainer.getFullName()+".charge",
                                                 particleContainer.getRefCharge() );
                            }
                        }
                    }
                } // if writeSecParticleInf
            } // store refit information end }}}

        }
        typedef std::unique_ptr<pat::CompositeCandidateCollection> ccc_pointer;
        edm::OrphanHandle<pat::CompositeCandidateCollection> ccHandle =
            ev.put( ccc_pointer( ccList ), name );
        for ( i = 0; i < n; ++i )
        {
            const BPHRecoCandidate* ptr = list[i].get();
            edm::Ref<pat::CompositeCandidateCollection> ccRef( ccHandle, i );
            ccRefMap[ptr] = ccRef;
        }
        return ccHandle;
    }
};

#endif
