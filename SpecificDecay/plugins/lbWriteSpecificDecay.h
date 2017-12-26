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
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"

// use reco::DeDxData
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Ref.h"

#include <string>
#include <vector>
#include <map>
#include <map>
#include <iostream>
#include <fstream>

#include "BPHAnalysis/SpecificDecay/interface/usefulFuncs.h"

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

    std::string bsPointLabel;
    std::string pVertexLabel;
    std::string patMuonLabel;
    std::string ccCandsLabel;
    std::string pfCandsLabel;
    std::string pcCandsLabel;
    std::string gpCandsLabel;
    std::string dedxHrmLabel;
    std::string dedxPLHLabel;
    //std::string bsLabel;

    // token wrappers to allow running both on "old" and "new" CMSSW versions
    BPHTokenWrapper< reco::BeamSpot                            > bsPointToken;
    BPHTokenWrapper< std::vector<reco::Vertex>                 > pVertexToken;
    BPHTokenWrapper< pat::MuonCollection                       > patMuonToken;
    BPHTokenWrapper< std::vector<pat::CompositeCandidate     > > ccCandsToken;
    BPHTokenWrapper< std::vector<reco::PFCandidate           > > pfCandsToken;
    BPHTokenWrapper< std::vector<BPHTrackReference::candidate> > pcCandsToken;
    BPHTokenWrapper< std::vector<pat::GenericParticle        > > gpCandsToken;
    BPHTokenWrapper< edm::ValueMap<reco::DeDxData            > > dedxHrmToken;
    BPHTokenWrapper< edm::ValueMap<reco::DeDxData            > > dedxPLHToken;
    //BPHTokenWrapper< reco::BeamSpot > bsToken;


    bool useBS;
    bool usePV;
    bool usePM;
    bool useCC;
    bool usePF;
    bool usePC;
    bool useGP;
    bool useHrm;
    bool usePLH;

// The label used in output product
    std::string     oniaName;
    std::string     Lam0Name;
    std::string LbToLam0Name;


    enum recoType { Onia, Psi1, Psi2, Lam0, LbToLam0 };
    enum  parType { ptMin, etaMax,
                    mPsiMin, mPsiMax,
                    mLam0Min, mLam0Max,
                    massMin, massMax, probMin, mFitMin, mFitMax,
                    constrMass, constrSigma, constrMJPsi, writeCandidate,
                  };

    std::map<std::string,recoType> rMap;
    std::map<std::string, parType> pMap;
    std::map<std::string, parType> fMap;
    std::map< recoType, std::map<parType,double> >  parMap;


    bool recoOnia     ;
    bool recoLam0     ;
    bool recoLbToLam0 ;

    bool writeOnia    ;
    bool writeLam0    ;
    bool writeLbToLam0;

    bool writeVertex;
    bool writeMomentum;

    std::vector<BPHPlusMinusConstCandPtr> lFull;
    std::vector<BPHPlusMinusConstCandPtr> lJPsi;
    std::vector<BPHPlusMinusConstCandPtr> lLam0;
    std::vector<BPHRecoConstCandPtr>      lLbToLam0;

    // the connection between JPsi & lFull.
    // lFull is the mumu candidate which will be recorded in storage.
    // JPsi is the particle you regard it as jpsi candidate.
    std::map<const BPHRecoCandidate*,const BPHRecoCandidate*> jPsiOMap;
    // the connection between lFull & pvRef
    typedef edm::Ref< std::vector<reco::Vertex> > vertex_ref;
    std::map<const BPHRecoCandidate*,vertex_ref> pvRefMap;
    // after candidates are "put" in storage, the function 'write' will create Ref object.
    // and the
    typedef edm::Ref< pat::CompositeCandidateCollection > compcc_ref;
    std::map<const BPHRecoCandidate*,compcc_ref> ccRefMap;

    void setRecoParameters( const edm::ParameterSet& ps );
    template <class T>
    edm::OrphanHandle<pat::CompositeCandidateCollection> write( edm::Event& ev,
            const std::vector<T>& list, const std::string& name, bool writeDownThisEvent )
    {
        pat::CompositeCandidateCollection* ccList =
            new pat::CompositeCandidateCollection;
        int i;
        int n = list.size();
        //int listTag[ n ];
        //for ( int kk = 0; kk < n; ++kk )
        //    listTag[kk] = -1;

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

        edm::Handle<edm::ValueMap<reco::DeDxData>> dedxHrmHandle;
        edm::Handle<edm::ValueMap<reco::DeDxData>> dedxPLHHandle;
        edm::Handle< std::vector<pat::GenericParticle> > gpHandle;
        edm::Handle< std::vector<reco::PFCandidate> > pfHandle;
        edm::Handle< reco::BeamSpot > bsHandle;

        if ( usePLH )
            dedxPLHToken.get( ev, dedxPLHHandle );
        if ( useHrm )
            dedxHrmToken.get( ev, dedxHrmHandle );
        if ( usePLH || useHrm )
        {
            if ( usePF )
                pfCandsToken.get( ev, pfHandle );
            if ( useGP )
                gpCandsToken.get( ev, gpHandle );
        }
        if ( useBS ) // used to find Impact parameter for secondary particle
            if ( list.size() )
                if ( !list[0]->compNames().size() )
                    bsPointToken.get( ev, bsHandle );

        if ( writeDownThisEvent )
        {
            for ( i = 0; i < n; ++i )
            {
                const T& ptr = list[i];
                // bool-int to check if it is needed to write event.
                int dddPassed = 0;
                ccList->push_back( ptr->composite() );
                pat::CompositeCandidate& cc = ccList->back();
                //pat::CompositeCandidate cc = ptr->composite();
                if ( ( pvrIter = pvRefMap.find( ptr.get() ) ) != pvrIend )
                    cc.addUserData ( "primaryVertex", pvrIter->second );
                const std::vector<std::string>& cNames = ptr->compNames();
                int j = 0;
                int m = cNames.size();
                while ( j < m )
                {
                    const std::string& compName = cNames[j++];
                    if ( !ptr->getComp( compName ) ) continue;
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
                if ( !ptr->isValidFit() )  continue;

                const RefCountedKinematicParticle kinPart = ptr->currentParticle();
                const           KinematicState    kinStat = kinPart->currentState();
                cc.addUserFloat( "fitMass", kinStat.mass() );
                std::vector<usefulFuncs::recoParticleInfo> myParticleList;
                myParticleList.reserve( cc.numberOfDaughters()+1 );
                if ( writeMomentum )
                {
                    // store refit particle momentum
                    cc.addUserData ( "fitMomentum", kinStat.kinematicParameters().momentum() );

                    // store refit daughter momentum
                    std::vector<RefCountedKinematicParticle> daughters = ptr->kinematicTree()->finalStateParticles();
                    std::vector<RefCountedKinematicParticle> compDaugh;


                    // find RefCountedKinematicParticle

                    if ( ptr->seckinematicTree().get() )
                    //if ( ptr->seckinematicTree()->isValid() )
                    if ( ptr->compNames().size() )
                    {
                        // the last one is the constrained particle, it is virtual particle, need to be removed.
                        daughters.erase( --daughters.end() );
                        compDaugh = ptr->seckinematicTree()->finalStateParticles();
                    }


                    //if ( !(ptr->compNames().size()) )
                    //{
                    //    std::cout << "d : c : dF = " << daughters.size() <<","<< compDaugh.size() <<","<< ptr->daughFull().size() << std::endl;
                    //    for ( const auto& dN : ptr->daugNames() )
                    //        std::cout << "---" << dN;
                    //    std::vector<RefCountedKinematicParticle> cDaugh;
                    //    if ( ptr->seckinematicTree().get() )
                    //    {
                    //        cDaugh = ptr->seckinematicTree()->finalStateParticles();
                    //        for ( int i=0; i<2;++i )
                    //            if ( fabs( cDaugh[i]->currentState().mass() - daughters[i]->currentState().mass() ) > 0.000000000001 )
                    //                std::cout << "sec fit & first fit results are different!!\n";
                    //    }
                    //    std::cout << std::endl;
                    //}
                    //else
                    //{
                    //    std::cout << "d : c : dF = " << daughters.size() <<","<< compDaugh.size() <<","<< ptr->daughFull().size() << std::endl;
                    //    for ( const auto& dN : ptr->daugNames() )
                    //        std::cout << "---" << dN;
                    //    for ( const auto& cN : ptr->compNames() )
                    //        std::cout << "..." << cN;
                    //    std::cout << std::endl;
//
                    //}
                    // all RefCountedKinematicParticle are found, check if it losts particle?
                    if ( (daughters.size()+compDaugh.size()) != ptr->daughFull().size() ) continue;

                    //get the daughterFull name stored in BPHRecoCandidate. {{{
                    const std::vector<std::string>& cNames_ = ptr->compNames();
                    const std::vector<std::string>& dNames_ = ptr->daugNames();
                    std::map<const reco::Candidate*, std::string> componentList;
                    for ( const std::string& cName_ : cNames_ )
                    {
                        if ( !ptr->getComp(cName_) ) continue;
                        const BPHRecoCandidate* cand_ = ptr->getComp( cName_ ).get();
                        const std::vector<std::string>& dNames__ = cand_->daugNames();
                        for ( const std::string& dName__ : dNames__ )
                        {
                            const reco::Candidate* cand__ = cand_->getDaug( dName__ );
                            myParticleList.emplace_back( cName_, dName__, cand__ );
                        }
                    }
                    for ( const std::string& dName_ : dNames_ )
                    {
                        const reco::Candidate* cand_ = ptr->getDaug( dName_ );
                        myParticleList.emplace_back( "", dName_, cand_ );
                    }
                    // get daughFull name end }}}

                    int totalCharge = 0;
                    std::map<const reco::Candidate*, std::string> componentList_ = componentList;
                    std::map< const RefCountedKinematicParticle, std::pair<const reco::Candidate*, std::string> > particleList;

                    //use pt and charge to match the reco::Candidate* and RefCountedKinematicParticle
                    //And the constrained name is used to separate virtual particle and other particle
                    const std::string& constrName = ptr->getConstrName();
                    for ( usefulFuncs::recoParticleInfo& _oldParticle : myParticleList )
                    {

                        std::vector<RefCountedKinematicParticle> lists_;
                        if ( _oldParticle.getCompName() == constrName && !constrName.empty() )
                            lists_ = compDaugh;
                        else
                            lists_ = daughters;

                        int charge_ = _oldParticle.getRecoCharge();
                        double minDeltaR2 = 999.;
                        for ( RefCountedKinematicParticle& _cDau : lists_ )
                        {
                            // use charge and pt to distinguish particle
                            if ( _oldParticle.getRecoCharge() != charge_ ) continue;
                            if ( usefulFuncs::candidateMatch( _cDau, _oldParticle.getRecoParticle(), minDeltaR2 ) )
                                _oldParticle.setRefitParticle( _cDau );
                        }
                        totalCharge += _oldParticle.getRecoCharge();
                    }

                    // particle name - RefCountedKinematicParticle pair is prepared.
                    if ( myParticleList.size() != ptr->daughFull().size() )
                    {
                        printf( "-----bug: particleList and daughters don't owns the same size! in write() \n" );
                        continue;
                    }
                    if ( totalCharge != 0 )
                    {
                        printf( "-----bug: the built particle is not neutral! in write()\n" );
                        continue;
                    }

                    // add refit momentum of daughters

                    dddPassed = myParticleList.size();
                    for ( const usefulFuncs::recoParticleInfo& particleContainer : myParticleList )
                    {
                        if ( !particleContainer.refptr.get() ) continue;
                        cc.addUserData ( particleContainer.getFullName()+".fitMom",
                                         particleContainer.getRefitMom() );

                        std::unique_ptr<GlobalPoint> myReferencePoint(nullptr);
                        if ( ptr->compNames().size() ) // if it is Lambda0_b or Bs ( find PV )
                        {
                            if ( !ptr->getComp("JPsi") ) continue;
                            const BPHRecoCandidate* jpsi = ptr->getComp( "JPsi" ).get();
                            const reco::Vertex* _pv = nullptr;

                            // find the jpsi in lFull. And use the pvMap to find pv.
                            if (( jpoIter = jPsiOMap.find(jpsi) ) != jpoIend )
                                jpsi = jpoIter->second;
                            if (( pvrIter = pvRefMap.find(jpsi) ) != pvRefMap.end() )
                                _pv = pvrIter->second.get();
                            if ( !_pv ) continue;
                            if ( !_pv->isValid() ) continue;
                            std::unique_ptr<GlobalPoint> pvPoint( new GlobalPoint( _pv->x(), _pv->y(), _pv->z() ) );
                            myReferencePoint = std::move( pvPoint );
                        }
                        else if ( useBS )// if it is secondary candidate like Lam0 or Kshort or JPsi
                        {
                            if ( !bsHandle.isValid() ) continue;
                            std::unique_ptr<GlobalPoint> bsPoint( new GlobalPoint( bsHandle->x0(), bsHandle->y0(), bsHandle->z0() ) );
                            myReferencePoint = std::move( bsPoint );
                        }


                        const reco::TransientTrack& newTT = particleContainer.getRefitParticle()->refittedTransientTrack();
                        if ( !newTT.isValid() ) continue;
                        TrajectoryStateClosestToPoint traj = newTT.trajectoryStateClosestToPoint( *myReferencePoint );
                        float IPt ( traj.perigeeParameters().transverseImpactParameter() );
                        float IPt_err ( traj.perigeeError().transverseImpactParameterError() );


                        // !!! refit result preselection
                        //if ( IPt/IPt_err < 1. ) continue;
                        const GlobalPoint& pVtx = ptr->currentDecayVertex()->position();
                        const GlobalPoint& pPv  = *myReferencePoint.get();
                        GlobalVector distDiff( pVtx.x()-pPv.x(), pVtx.y()-pPv.y(), pVtx.z()-pPv.z() );
                        //GlobalVector momentum(

                        // !!! refit result preselection end

                        cc.addUserFloat ( particleContainer.getFullName()+".IPt", IPt );
                        cc.addUserFloat ( particleContainer.getFullName()+".IPt.Error", IPt_err );


                        reco::TrackRef trackID;
                        if ( trackID.isNull() ) trackID = usefulFuncs::getTrackRefFromGP( *particleContainer.getRecoParticle() );
                        if ( trackID.isNull() ) trackID = usefulFuncs::getTrackRefFromPF( *particleContainer.getRecoParticle() );
                        if ( trackID.isNull() ) trackID = usefulFuncs::getTrackRefFromRC( *particleContainer.getRecoParticle() );
                        if ( trackID.isNonnull() )
                        {
                                if ( useHrm )
                                {
                                    const edm::ValueMap<reco::DeDxData>& dedxHrmMap = *( dedxHrmHandle.product() );
                                    const reco::DeDxData& _dedxHrm = dedxHrmMap[ trackID ];
                                    if ( _dedxHrm.dEdx() )
                                        cc.addUserFloat( particleContainer.getFullName()+".dEdx.Harmonic", _dedxHrm.dEdx() );
                                }
                                if ( usePLH )
                                {
                                    const edm::ValueMap<reco::DeDxData>& dedxPLHMap = *( dedxPLHHandle.product() );
                                    const reco::DeDxData& _dedxPLH = dedxPLHMap[ trackID ];
                                    if ( _dedxPLH.dEdx() )
                                        cc.addUserFloat( particleContainer.getFullName()+".dEdx.pixelHrm", _dedxPLH.dEdx() );
                                }
                        }



                        --dddPassed;
                    } // end of particleContainer
                    //for ( const usefulFuncs::recoParticleInfo& particleContainer : myParticleList )
                    //    std::cout << ", " << particleContainer.getFullName();
                    //std::cout << std::endl;

                } // if writeMomentum end

                // store refit information end }}}

                // if all daughter information are stored
                //if ( dddPassed == 0 )
                    //ccList->push_back( cc );
                //listTag[i] = dddPassed;
            } // run over all candidate end
        } // if writeDownThisEvent
        typedef std::unique_ptr<pat::CompositeCandidateCollection> ccc_pointer;
        edm::OrphanHandle<pat::CompositeCandidateCollection> ccHandle =
            ev.put( ccc_pointer( ccList ), name );
        // after candidates are "put", start to obtain Ref object
        // ex : lFull, lJPsi, lLam0, lLbToJPsiLam0
        for ( i = 0; i < n; ++i )
        {
            //if ( listTag[i] != 0 ) continue;
            //if ( listTag[i] ) listTag[i] = listTag[i];
            const BPHRecoCandidate* ptr = list[i].get();
            edm::Ref<pat::CompositeCandidateCollection> ccRef( ccHandle, i );
            ccRefMap[ptr] = ccRef;
        }
        return ccHandle;
    }
};

#endif
