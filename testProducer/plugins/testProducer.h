#ifndef BPHAnalysis_testAnalyzer_testProducer_h
#define BPHAnalysis_testAnalyzer_testProducer_h

#include "BPHAnalysis/RecoDecay/interface/BPHAnalyzerTokenWrapper.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"

#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "TTree.h"
#include "BPHAnalysis/testAnalyzer/interface/format.h"
class TH1F;
class BPHRecoCandidate;

class testProducer:
    public BPHAnalyzerWrapper<BPHModuleWrapper::one_analyzer>
{

public:

    explicit testProducer( const edm::ParameterSet& ps );
    virtual ~testProducer();

    static void fillDescriptions( edm::ConfigurationDescriptions& descriptions );

    virtual void beginJob();
    //virtual void produce( edm::Event& ev, const edm::EventSetup& es );
    virtual void analyze( const edm::Event& ev, const edm::EventSetup& es );
    virtual void fill   ( const edm::Event& ev, const edm::EventSetup& es );
    virtual void endJob();

private:


// The label used in reading data
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

    enum recoType { Onia, Psi1, Psi2, Lam0, TkTk, LbToLam0, LbToTkTk };
    enum  parType { ptMin, etaMax,
                    mPsiMin, mPsiMax, mLam0Min, mLam0Max,
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

    bool writeOnia    ;
    bool writeLam0    ;
    bool writeTkTk    ;
    bool writeLbToLam0;
    bool writeLbToTkTk;

    bool writeVertex;
    bool writeMomentum;

    std::vector<BPHPlusMinusConstCandPtr> lFull;
    std::vector<BPHPlusMinusConstCandPtr> lJPsi;
    std::vector<BPHPlusMinusConstCandPtr> lLam0;
    std::vector<BPHPlusMinusConstCandPtr> lTkTk;
    std::vector<BPHRecoConstCandPtr>      lLbToLam0;
    std::vector<BPHRecoConstCandPtr>      lLbToTkTk;
    std::vector<BPHRecoConstCandPtr>      laLbToTkTk;

    std::map<const BPHRecoCandidate*,const BPHRecoCandidate*> jPsiOMap;
    typedef edm::Ref< std::vector<reco::Vertex> > vertex_ref;
    std::map<const BPHRecoCandidate*,vertex_ref> pvRefMap;
    typedef edm::Ref< pat::CompositeCandidateCollection > compcc_ref;
    std::map<const BPHRecoCandidate*,compcc_ref> ccRefMap;
    std::map<const BPHRecoCandidate*,const BPHRecoCandidate*> secVtxParticleMap;

    // my code added:{{{
    //

    edm::Service<TFileService> fs;
    TTree *jpsiTree, *lam0Tree, *lambTree, *lbv1Tree, *lbv2Tree;

    JpsiBranches       jpsiBr;
    Lam0Branches       lam0Br;
    LambToLam0Branches lambBr;
    LambToTkTkBranches lbv1Br;
    LambToTkTkBranches lbv2Br;
    UInt_t eventNo;


    struct SecRecoResult// {{{
    {
    public:
        SecRecoResult(const reco::Vertex& SecVtx, const double& SecMass, const GlobalVector* SecMom, int dTag) : VERTEX(SecVtx), MASS(SecMass), flag(true), debugTag(dTag)
        { MOMENTUM = new GlobalVector( *SecMom ); }
        virtual ~SecRecoResult()
        { delete MOMENTUM;  MOMENTUM = NULL; }
        SecRecoResult( const SecRecoResult& in ) : VERTEX(in.VERTEX), MASS(in.MASS), flag(in.flag), debugTag(in.debugTag)
        { MOMENTUM = new GlobalVector( *(in.MOMENTUM) ); }
        SecRecoResult( bool fff, int dTag ) : VERTEX(), MASS(-999.), flag(fff), debugTag(dTag)
        { MOMENTUM = NULL; }
        SecRecoResult() : VERTEX(), MASS(-999.), debugTag(0)
        { printf(" Hey, you need to initialize '''SecRecoResult''' \n"); }


        const reco::Vertex* vertex() const
        { return &VERTEX; }
        double mass() const
        { return MASS; }
        const GlobalVector* momentum() const
        { return MOMENTUM; }
        bool isValid() const
        { return flag; }
        int bugCode() const
        { return debugTag; }

    
    private:
        reco::Vertex  VERTEX;
        double MASS;
        GlobalVector* MOMENTUM;
        bool flag;
        int debugTag;
    }; // SecRecoResult end }}}



    // _inputName is the particle used in the reconstruction.
    // Totally process: J/Psi + inputParticle.
    // if not set, reconstruction the process: J/Psi + Tk + Tk.
    // This reconstruction uses J/Psi constraint.
    // The input lLamb need to contain a J/psi compositeCanddiate
    //template <typename T> const reco::Vertex* findPrimaryVertex( const T& list );

    SecRecoResult secondaryReconstruction(const BPHRecoCandidate*& cand, const std::vector<std::string>* _inputName, MultiTrackKinematicConstraint* constraint = 0);
    const BPHRecoCandidate secondVertexReconstruct_paolo( const edm::EventSetup& es, const BPHRecoCandidate* cand, const std::vector<std::string>* inputName, MultiTrackKinematicConstraint* kc = 0);
    const pat::CompositeCandidate* compCandRef(const BPHRecoCandidate* cand);
    //void fillTree();
    const reco::Vertex* findPrimaryVertex( const BPHRecoCandidate* list );
    std::string findCompName(const std::string& input );
    std::string findDaugName(const std::string& input );
    void fillTree( const edm::EventSetup& es );


    std::map<std::string, float>       _dmassMAP;
    std::map<std::string, float>       _sigmaMAP;


    // my code added end }}}






    void setRecoParameters( const edm::ParameterSet& ps );

    // write {{{
    template <class T>
    edm::OrphanHandle<pat::CompositeCandidateCollection> write( edm::Event& ev, const edm::EventSetup& es,
            const std::vector<T>& list, const std::string& name , 
            bool useJpsiConstrVertexFit = false) 
            //bool approveSecondaryReconstruction = false)
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

// Find out ccrIter and ccrIter in JpsoIter from input list[], the record it {{{
            // if find ccrIter, record it
            // &&
            // if fint ccrIter recorded in JpsoIter, record it, too.
            // Input a particle,  Add reference subparticle.
            //    first find out the name of subparticle.
            //    second search for the map, what is map of the subparticle stored.
            //    if find out the map stored. add the ref-subparticle to UserData.
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
// Find out ccrIter and ccrIter in JpsoIter, the record it end }}}
            const BPHPlusMinusCandidate* pmp =
                dynamic_cast<const BPHPlusMinusCandidate*>( ptr.get() );
            if ( pmp != 0 ) cc.addUserData( "cowboy", pmp->isCowboy() );
            if ( ptr->isEmpty() )
            {
                if ( writeVertex ) cc.addUserData( "vertex" , ptr->vertex() );
                continue;
            }


// My jpsi constraint fitVertex
                
            if ( useJpsiConstrVertexFit )//{{{
            {
                ParticleMass jpsi_mass = 3.096916;
                TwoTrackMassKinematicConstraint *jpsi_const = new TwoTrackMassKinematicConstraint(jpsi_mass);
                RefCountedKinematicTree compTree = ptr.get()->kinematicTree("JPsi", jpsi_const);
                if ( writeVertex ) cc.addUserData( "fitVertex", reco::Vertex( *compTree->currentDecayVertex() ) );
                if ( ptr->isValidFit() )
                {
                    const RefCountedKinematicParticle kinPart = compTree->currentParticle();
                    const           KinematicState    kinStat = kinPart->currentState();
                    cc.addUserFloat( "fitMass", kinStat.mass() );
                    cc.addUserData ( "fitMomentum",
                                     kinStat.kinematicParameters().momentum() );
                }

                // store transient track to fit at further analysis
                



                bool approveSecondaryReconstruction = true;
                if ( approveSecondaryReconstruction )
                {
                    const BPHRecoCandidate* cand = ptr.get();
                    const reco::Candidate* mPos   = cand->originalReco(
                                                    cand->getDaug( "JPsi/MuPos" ) );
                    const reco::Candidate* mNeg   = cand->originalReco(
                                                    cand->getDaug( "JPsi/MuNeg" ) );
                    const reco::Candidate* proton = cand->originalReco(
                                                    cand->getDaug( "Proton"     ) );
                    if ( !proton ) printf("no proton candidate\n");
                    else
                    {
                        BPHRecoCandidatePtr njp( new BPHPlusMinusCandidate( &es ) );
                        njp->add( "MuPos", mPos,
                                  BPHParticleMasses::muonMass,
                                  BPHParticleMasses::muonMSigma );
                        njp->add( "MuNeg", mNeg,
                                  BPHParticleMasses::muonMass,
                                  BPHParticleMasses::muonMSigma );
                        BPHRecoCandidate nSecFit( &es );
                        nSecFit.add( "JPsi", njp );
                        nSecFit.add( "Proton", proton,
                                  BPHParticleMasses::protonMass,
                                  BPHParticleMasses::protonMSigma );
                        compTree = nSecFit.kinematicTree( "JPsi", jpsi_const );
                        RefCountedKinematicTree compTree = ptr.get()->kinematicTree("JPsi", jpsi_const);
                        if ( writeVertex ) cc.addUserData( "secfitVertex", reco::Vertex( *compTree->currentDecayVertex() ) );
                        if ( ptr->isValidFit() )
                        {
                            printf("secondary vertex fit processed\n");
                            const RefCountedKinematicParticle kinPart = compTree->currentParticle();
                            const           KinematicState    kinStat = kinPart->currentState();
                            cc.addUserFloat( "secfitMass", kinStat.mass() );
                            cc.addUserData ( "secfitMomentum",
                                             kinStat.kinematicParameters().momentum() );
                        }
                    }
                }
            }// if(useJpsiConstrVertexFit) end}}}
            else
            {
                if ( writeVertex ) cc.addUserData( "fitVertex", reco::Vertex( *ptr->currentDecayVertex() ) );
                if ( ptr->isValidFit() )
                {
                    const RefCountedKinematicParticle kinPart = ptr->currentParticle();
                    const           KinematicState    kinStat = kinPart->currentState();
                    cc.addUserFloat( "fitMass", kinStat.mass() );
                    cc.addUserData ( "fitMomentum",
                                     kinStat.kinematicParameters().momentum() );

                }
            }

        }
        typedef std::unique_ptr<pat::CompositeCandidateCollection> ccc_pointer;
        edm::OrphanHandle<pat::CompositeCandidateCollection> ccHandle =
            ev.put( ccc_pointer( ccList ), name );

        //  putProducts().push_back(std::make_pair(edp, &desc));

        // product.release(); // The object has been copied into the Wrapper.
        // The old copy must be deleted, so we cannot release ownership.

        //return(OrphanHandle<PROD>(prod, makeProductID(desc)));


        for ( i = 0; i < n; ++i )
        {
            const BPHRecoCandidate* ptr = list[i].get();
            edm::Ref<pat::CompositeCandidateCollection> ccRef( ccHandle, i );
            ccRefMap[ptr] = ccRef;
        }
        return ccHandle;
    }

    // write end }}}
};

std::string testProducer::findDaugName(const std::string& input )
{
    std::string daugName;
    if ( input.empty() ) return std::string();
    size_t slash = input.find_first_of("/");
    if ( slash > input.size() ) 
        daugName = input;
    else
        for(size_t i=slash+1;i < input.size() ; ++i)
            daugName += input[i];
    return daugName;
}
std::string testProducer::findCompName(const  std::string& input )
{
    std::string compName;
    if ( input.empty() ) return std::string();
    size_t slash = input.find_first_of("/");
    if ( slash > input.size() ) 
        return std::string();
    else
        for(size_t i=0;i < slash ; ++i)
            compName += input[i];
    return compName;
}
// secondaryReconstruction functions {{{
testProducer::SecRecoResult testProducer::secondaryReconstruction(const BPHRecoCandidate*& cand, const std::vector<std::string>* _inputName, MultiTrackKinematicConstraint* constraint)
{
    // _inputName example:
    // vector<string> inputName;
    // inputName.push_back("JPsi/MuPos");
    // inputName.push_back("JPsi/MuNeg");
    // inputName.push_back("Lam0/Proton");
    // intpuName.push_back("Lam0/Pion");
    std::vector<std::string> compName;
    // find composite particle name
    for (const auto& inName : *_inputName)
    {
        std::string a;
        a = findCompName(inName);
        bool addName = true;
        if ( a.empty() ) addName = false;
        for( const auto& checkName : compName )
            if ( a == checkName ) addName = false;
        if ( addName ) compName.push_back(a);
    }
    std::vector<std::string> daugName;
    // find daughter particle name
    for (const auto& inName : *_inputName)
        if ( findCompName(inName).empty() ) daugName.push_back(inName);

    std::vector<RefCountedKinematicParticle> kComp;


    // debugCode: (written in int, replacing 8 boolin.)
    // 0. all events
    // 1. sucess to reconstruct composite particle
    // 2. sucess to add daughter particle in kComp
    // 3. sucess to fit vertex (compTree is not empty)
    // 4. sucess to get kinematicState
    // 5. vertex probability is bigger than 0
    int debugCode = 0;
    debugCode += 1 << 0;
    // --------------------------------------------------
    // perform vertex fit for composite particle
    if ( constraint ) // if there is jpsi constraint, erase jpsi composite particle.
    {
        compName.erase( compName.begin() );     
        kComp = cand->getComp( "JPsi" ).get()->kinParticles();
    }
    for (const auto& ccName : compName)
    {
        const BPHRecoCandidate* _cand = cand->getComp( ccName ).get();
        std::vector<RefCountedKinematicParticle> _kComp = _cand->kinParticles();

        KinematicParticleVertexFitter _vtxFitter;
        RefCountedKinematicTree _compTree = _vtxFitter.fit( _kComp );
        if (_compTree->isEmpty()) return SecRecoResult(false, debugCode);
        _compTree->movePointerToTheTop();
        const RefCountedKinematicParticle  _kPart  = _compTree->currentParticle();
        const RefCountedKinematicVertex    _kVtx   = _compTree->currentDecayVertex();
        const KinematicState&              _kState = _kPart   ->currentState();
        if (!_kState.isValid()) return SecRecoResult(false, debugCode);

        double _vtxProb = TMath::Prob( _kVtx->chiSquared(), _kVtx->degreesOfFreedom() );
        if (_vtxProb <=0.02) return SecRecoResult(false, debugCode); 

        VirtualKinematicParticleFactory vFactory;
        float chi = _kVtx->chiSquared();
        float ndf = _kVtx->degreesOfFreedom();
        kComp.push_back(vFactory.particle(_kState,chi,ndf,_kPart)); 
    }
    debugCode += 1 << 1;
    // --------------------------------------------
    // add daughter particles in the final particle
    for (const auto& ddName : daugName)
    {
        const reco::Candidate* cand_daug = cand->getDaug(ddName);
        if(cand_daug != 0 )
        {
            reco::TransientTrack* tTrack   = cand->getTransientTrack( cand_daug );

            KinematicParticleFactoryFromTransientTrack pFactory;
            float chi = cand_daug->vertexChi2();
            float ndf = cand_daug->vertexNdof();

            kComp.push_back( pFactory.particle(*tTrack, _dmassMAP[ddName], chi, ndf, _sigmaMAP[ddName]) );
        
        }
        else printf( "daughter %s cannot be found!\n", ddName.c_str() );
    }
    debugCode += 1 << 2;



    KinematicConstrainedVertexFitter kcvFitter;
    RefCountedKinematicTree compTree;
    if( constraint )
        compTree = kcvFitter.fit(kComp, constraint);
    else
        compTree = kcvFitter.fit(kComp);

    if (compTree->isEmpty()) return SecRecoResult(false, debugCode);
    debugCode += 1 << 3;
    compTree->movePointerToTheTop();
    const RefCountedKinematicParticle kinPart = compTree->currentParticle();
    const RefCountedKinematicVertex secVertex = compTree->currentDecayVertex();
    const           KinematicState&   kinStat = kinPart->currentState();
    const GlobalVector secMom = kinStat.kinematicParameters().momentum();
    const double secMass = kinStat.mass();
    if (!kinStat.isValid()) return SecRecoResult(false, debugCode);
    debugCode += 1 << 4;

    double vtxProb = TMath::Prob( secVertex->chiSquared(), secVertex->degreesOfFreedom() );
    if (vtxProb <=0.) return SecRecoResult(false, debugCode);
    debugCode += 1 << 5;

    return SecRecoResult( *secVertex, secMass, &secMom , debugCode);
}
const pat::CompositeCandidate* testProducer::compCandRef(const BPHRecoCandidate* cand)
{

    // --------------------------------------------------
    // refit the composite candidate in the particle.
    const pat::CompositeCandidate* cand_ref = NULL;
    std::map<const BPHRecoCandidate*,compcc_ref>::const_iterator ccrIter;
    std::map<const BPHRecoCandidate*,compcc_ref>::const_iterator ccrIend = ccRefMap.end();
    if ( ( ccrIter = ccRefMap.find( cand ) ) != ccrIend )
    {
        compcc_ref cref = ccrIter->second;
        if ( cref.isNonnull() ) cand_ref = ccrIter->second.get();
    }

    return cand_ref;
} 
const BPHRecoCandidate testProducer::secondVertexReconstruct_paolo( const edm::EventSetup& es, const BPHRecoCandidate* cand, const std::vector<std::string>* inputName, MultiTrackKinematicConstraint* kc)
{
    // input Name example:
    // vector<string> input;
    // input.push_back("JPsi/MuPos");
    // input.push_back("JPsi/MuNeg");
    // input.push_back("Proton");

    // What is this function doing:
    //
    //BPHRecoCandidate secParticleReco( &es );
    //

    //const reco::Candidate* mPos = cand->originalReco(
    //                              cand->getDaug( "JPsi/MuPos" ) );
    //const reco::Candidate* mNeg = cand->originalReco(
    //                              cand->getDaug( "JPsi/MuNeg" ) );
    //const reco::Candidate* kaon = cand->originalReco(
    //                              cand->getDaug( "kaon"       ) );
    //BPHRecoCandidatePtr njp( new BPHPlusMinusCandidate( &es ) );
    //njp->add( "MuPos", mPos,
    //          BPHParticleMasses::muonMass,
    //          BPHParticleMasses::muonMSigma );
    //njp->add( "MuNeg", mNeg,
    //          BPHParticleMasses::muonMass,
    //          BPHParticleMasses::muonMSigma );
    //BPHRecoCandidate nbu( &es );
    //nbu.add( "JPsi", njp );
    //nbu.add( "Kaon", kaon,
    //         BPHParticleMasses::kaonMass,
    //         BPHParticleMasses::kaonMSigma );
    //nbu.kinematicTree( "JPsi", 
    
    BPHRecoCandidatePtr compReco( new BPHPlusMinusCandidate( &es ) );
    
    std::string nnn = (*inputName)[0];
    //std::string _cName = findCompName( (*inputName)[0] );
    std::string _cName = findCompName( nnn );
    // find composite candidate name.
    // e.g.
    //     JPsi/MuPos, JPsi/MuNeg   -> JPsi
    //     Lam0/Proton,Lam0/Pion    -> Lam0
    std::vector<std::string> cName;
    for ( const auto& _cName : *inputName )
    { 
        std::string tmpcName = findCompName( _cName );
        bool ownSameCompName = false;
        for ( const auto& __cName : cName )
            if ( __cName == tmpcName )
                ownSameCompName = true;
        if ( !ownSameCompName )
            cName.push_back( tmpcName );
    }
    // find daughter, not form a particle.
    // e.g.
    //     Kaon, Proton
    std::vector<std::string> dName;
    for ( const auto& _dName : *inputName )
        if( findCompName(_dName).empty() )
            dName.push_back( _dName );

    
    // reconstruct a composite particle
    // e.g.
    //      JPsi, Lam0
    //BPHRecoCandidatePtr compPtr[cName.size()]( new BPHPlusMinusCandidate(&es) );
    BPHRecoCandidatePtr compPtr[cName.size()];
    for (size_t i = 0; i < cName.size(); ++i)
        compPtr[i] = BPHRecoCandidatePtr( new BPHPlusMinusCandidate( &es ) );
    for(size_t i=0;i != cName.size(); ++i)
    {
        for ( const auto& _dInput : *inputName )
        {
            std::string _ccName = findCompName( _dInput );
            // if find composite particle
            if ( !_ccName.empty() && cName[i] == _ccName )
            {
                const reco::Candidate* dComp= cand->originalReco(
                                              cand->getDaug( _dInput.c_str() ) );
                std::string dcName = findDaugName( _dInput );
                compPtr[i]->add( dcName.c_str(), dComp, _dmassMAP[dcName], _sigmaMAP[dcName] );
            }
        }
    }

    BPHRecoCandidate finalComp( &es );
    if ( kc )
    {
        // if there is kinematic constraint
        // means that there is JPsi particle.
        for(size_t i=0;i<cName.size();++i)
            finalComp.add(findDaugName( cName[i] ).c_str(), compPtr[i]);
        for(size_t i=0;i<dName.size();++i)
        {
            const reco::Candidate* daug = cand->originalReco(
                                          cand->getDaug( dName[i].c_str() ) );
            std::string ddName = findDaugName( dName[i] );
            finalComp.add( dName[i].c_str(), daug, _dmassMAP[ddName], _sigmaMAP[ddName] );
        }
        finalComp.kinematicTree( "JPsi", kc );
    }
    else
    {
        // if there is no kinematic constraint,
        // all particles are at the same level.
        //
        // e.g. 
        //     mu+, mu-, proton. three of them has the same contribution on vertex.
        for(size_t i=0;i<inputName->size();++i)
        {
            const reco::Candidate* daug = cand->originalReco(
                                          cand->getDaug( (*inputName)[i].c_str() ) );

            std::string ddName = findDaugName( (*inputName)[i] );
            finalComp.add( ddName.c_str(), daug, _dmassMAP[ddName], _sigmaMAP[ddName] );
        }
    }


        return finalComp;
}
// secondaryReconstruction functions }}}

//template <typename T>
//const reco::Vertex* testProducer::findPrimaryVertex( const T& list )
const reco::Vertex* testProducer::findPrimaryVertex( const BPHRecoCandidate* cand )
{
    std::map<const BPHRecoCandidate*,vertex_ref>::const_iterator pvrIter;
    std::map<const BPHRecoCandidate*,vertex_ref>::const_iterator pvrIend = pvRefMap.end();
    std::map<const BPHRecoCandidate*,const BPHRecoCandidate*>::const_iterator jpOiter;
    std::map<const BPHRecoCandidate*,const BPHRecoCandidate*>::const_iterator jpOiend = jPsiOMap.end  ();
    const BPHRecoCandidate* cand_onia = NULL;
    if ( ( jpOiter = jPsiOMap.find( cand ) ) != jpOiend )
        cand_onia = jpOiter->second;

    if ( ( pvrIter = pvRefMap.find( cand_onia ) ) != pvrIend )
    { return pvrIter->second.get(); }
    return NULL;
}
#endif
