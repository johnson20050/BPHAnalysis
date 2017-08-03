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

#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

class TH1F;
class BPHRecoCandidate;

class lbWriteSpecificDecay:
      public BPHAnalyzerWrapper<BPHModuleWrapper::one_producer> {

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
  BPHTokenWrapper< std::vector<reco::Vertex                > > pVertexToken;
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
  std::map< recoType, std::map<parType,double> > parMap;

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

  void setRecoParameters( const edm::ParameterSet& ps );

    // write {{{
    template <class T>
    edm::OrphanHandle<pat::CompositeCandidateCollection> write( edm::Event& ev, const edm::EventSetup& es, 
            const std::vector<T>& list, const std::string& name , 
            bool useJpsiConstrVertexFit = false,
            bool approveSecondaryReconstruction = false)
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
                if ( !ptr->isValidFit() ) continue;
                if ( compTree->isEmpty() ) continue;
                compTree->movePointerToTheTop();
                const RefCountedKinematicParticle kinPart = compTree->currentParticle();
                const           KinematicState    kinStat = kinPart->currentState();
                if ( !kinStat.isValid() ) continue;

                if ( writeVertex ) cc.addUserData( "fitVertex", reco::Vertex( *compTree->currentDecayVertex() ) );
                cc.addUserFloat( "fitMass", kinStat.mass() );
                cc.addUserData ( "fitMomentum",
                                 kinStat.kinematicParameters().momentum() );
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
                        RefCountedKinematicTree _compTree = ptr.get()->kinematicTree("JPsi", jpsi_const);
                        if ( !ptr->isValidFit() ) continue;
                        if ( _compTree->isEmpty() ) continue;
                        compTree->movePointerToTheTop();
                        const RefCountedKinematicParticle _kinPart = _compTree->currentParticle();
                        const           KinematicState    _kinStat = _kinPart->currentState();
                        if ( !_kinStat.isValid() ) continue;
                        if ( writeVertex ) cc.addUserData( "secfitVertex", reco::Vertex( *_compTree->currentDecayVertex() ) );
                        cc.addUserFloat( "secfitMass", _kinStat.mass() );
                        cc.addUserData ( "secfitMomentum",
                                         _kinStat.kinematicParameters().momentum() );
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
#endif
