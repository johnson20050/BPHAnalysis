#include "BPHAnalysis/SpecificDecay/plugins/lbWriteSpecificDecay.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "BPHAnalysis/RecoDecay/interface/BPHRecoBuilder.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHRecoCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHPlusMinusCandidate.h"
#include "BPHAnalysis/RecoDecay/interface/BPHMomentumSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHVertexSelect.h"
#include "BPHAnalysis/RecoDecay/interface/BPHTrackReference.h"

#include "BPHAnalysis/SpecificDecay/interface/BPHMuonPtSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMuonEtaSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticlePtSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleNeutralVeto.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHParticleMasses.h"
#include "BPHAnalysis/RecoDecay/interface/BPHMultiSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHChi2Select.h"

#include "BPHAnalysis/SpecificDecay/interface/BPHOniaToMuMuBuilder.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHLambda0ToPPiBuilder.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHLambda0_bToJPsiLambda0Builder.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <set>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

#define SET_PAR(TYPE,NAME,PSET) ( NAME = PSET.getParameter< TYPE >( #NAME ) )
// SET_PAR(string,xyz,ps);
// is equivalent to
// ( xyz = ps.getParameter< string >( "xyx" ) )

lbWriteSpecificDecay::lbWriteSpecificDecay( const edm::ParameterSet& ps ) {

    //  Check if there is the label used in reading data or not.
    useBS = ( SET_PAR( string, bsPointLabel, ps ) != "" );
    usePV = ( SET_PAR( string, pVertexLabel, ps ) != "" );
    usePM = ( SET_PAR( string, patMuonLabel, ps ) != "" );
    useCC = ( SET_PAR( string, ccCandsLabel, ps ) != "" );
    usePF = ( SET_PAR( string, pfCandsLabel, ps ) != "" );
    usePC = ( SET_PAR( string, pcCandsLabel, ps ) != "" );
    useGP = ( SET_PAR( string, gpCandsLabel, ps ) != "" );
    SET_PAR( string, dedxHrmLabel, ps );
    SET_PAR( string, dedxPLHLabel, ps );


    SET_PAR( bool  ,        isAOD, ps );
    // Set the label used in output product.
    SET_PAR( string,     oniaName, ps );
    SET_PAR( string,     Lam0Name, ps );
    SET_PAR( string, LbToLam0Name, ps );
    SET_PAR( bool  ,writeMomentum, ps );
    SET_PAR( bool  ,writeVertex  , ps );

    rMap["Onia"    ] = Onia;
    rMap["JPsi"    ] = Psi1;
    rMap["Psi2"    ] = Psi2;
    rMap["Lam0"    ] = Lam0;
    rMap["LbToLam0"] = LbToLam0;

    pMap["ptMin"      ] = ptMin;
    pMap["etaMax"     ] = etaMax;
    pMap["mPsiMin"    ] = mPsiMin;
    pMap["mPsiMax"    ] = mPsiMax;
    pMap["massMin"    ] = massMin;
    pMap["massMax"    ] = massMax;
    pMap["probMin"    ] = probMin;
    pMap["massFitMin" ] = mFitMin;
    pMap["massFitMax" ] = mFitMax;
    pMap["constrMass" ] = constrMass;
    pMap["constrSigma"] = constrSigma;

    fMap["constrMJPsi"   ] = constrMJPsi;
    fMap["writeCandidate"] = writeCandidate;

    recoOnia      =
    recoLam0      = writeLam0     =   true;
    recoLbToLam0  = writeLbToLam0 =  false;

    writeOnia = true;
    const vector<edm::ParameterSet> recoSelect =
        ps.getParameter< vector<edm::ParameterSet> >( "recoSelect" );
    int iSel;
    int nSel = recoSelect.size();
    for ( iSel = 0; iSel != nSel; ++iSel ) setRecoParameters( recoSelect[iSel] );
    if ( !recoOnia ) writeOnia = false;

    if (  recoLbToLam0 )  recoOnia =  recoLam0   = true;
    if ( writeLbToLam0 ) writeOnia = writeLam0   = true;

    // Get data by label/token
    if ( useBS ) consume< reco::BeamSpot                       >( bsPointToken,
                                                                  bsPointLabel );
    if ( usePV ) consume< vector<reco::Vertex                > >( pVertexToken,
                                                                  pVertexLabel );
    if ( usePM ) consume< pat::MuonCollection                  >( patMuonToken,
                                                                  patMuonLabel );
    if ( useCC ) consume< vector<pat::CompositeCandidate     > >( ccCandsToken,
                                                                  ccCandsLabel );
    if ( usePF ) consume< vector<reco::PFCandidate           > >( pfCandsToken,
                                                                  pfCandsLabel );
    if ( usePC ) consume< vector<BPHTrackReference::candidate> >( pcCandsToken,
                                                                  pcCandsLabel );
    if ( useGP ) consume< vector<pat::GenericParticle        > >( gpCandsToken,
                                                                  gpCandsLabel );

    if ( isAOD ) consume< edm::ValueMap<reco::DeDxData>        >( dedxHrmToken,
                                                                  dedxHrmLabel );
    if ( isAOD ) consume< edm::ValueMap<reco::DeDxData>        >( dedxPLHToken,
                                                                  dedxPLHLabel );

  if ( writeOnia     ) produces<pat::CompositeCandidateCollection>( oniaName     );
  if ( writeLam0     ) produces<pat::CompositeCandidateCollection>( Lam0Name     );
  if ( writeLbToLam0 ) produces<pat::CompositeCandidateCollection>( LbToLam0Name );

}


lbWriteSpecificDecay::~lbWriteSpecificDecay() {
}


void lbWriteSpecificDecay::fillDescriptions(
                            edm::ConfigurationDescriptions& descriptions ) {
    edm::ParameterSetDescription desc;
    desc.add<string>( "bsPointLabel", "" );
    desc.add<string>( "pVertexLabel", "" );
    desc.add<string>( "patMuonLabel", "" );
    desc.add<string>( "ccCandsLabel", "" );
    desc.add<string>( "pfCandsLabel", "" );
    desc.add<string>( "pcCandsLabel", "" );
    desc.add<string>( "gpCandsLabel", "" );
    desc.add<string>( "dedxHrmLabel", "" );
    desc.add<string>( "dedxPLHLabel", "" );
    desc.add<string>(     "oniaName",   "oniaCand" );
    desc.add<string>(     "Lam0Name",   "Lam0Cand" );
    desc.add<string>( "LbToLam0Name", "LbToLam0Fitted" );
    

    desc.add<bool>  ( "writeVertex"  , true );
    desc.add<bool>  ( "writeMomentum", true );
    desc.add<bool>  (         "isAOD",false );
    edm::ParameterSetDescription dpar;
    dpar.add<string>(           "name" );
    dpar.add<double>(          "ptMin", -2.0e35 );
    dpar.add<double>(         "etaMax", -2.0e35 );
    dpar.add<double>(        "mPsiMin", -2.0e35 );
    dpar.add<double>(        "mPsiMax", -2.0e35 );
    dpar.add<double>(        "massMin", -2.0e35 );
    dpar.add<double>(        "massMax", -2.0e35 );
    dpar.add<double>(        "probMin", -2.0e35 );
    dpar.add<double>(     "massFitMin", -2.0e35 );
    dpar.add<double>(     "massFitMax", -2.0e35 );
    dpar.add<double>(     "constrMass", -2.0e35 );
    dpar.add<double>(    "constrSigma", -2.0e35 );
    dpar.add<  bool>(    "constrMJPsi",    true );
    dpar.add<  bool>( "writeCandidate",    true );
    dpar.add<  bool>(          "isAOD",   false );

    vector<edm::ParameterSet> rpar;
    desc.addVPSet( "recoSelect", dpar, rpar );
    descriptions.add( "lbSpecificDecay", desc );
    return;
}


void lbWriteSpecificDecay::beginJob() {
  return;
}


void lbWriteSpecificDecay::produce( edm::Event& ev,
                                     const edm::EventSetup& es ) {
  fill( ev, es );
  //bool writeEvent = lLbToLam0.size();
  bool writeEvent = true;
  if ( writeOnia     ) write( ev, lFull     ,     oniaName , writeEvent ); 
  if ( writeLam0     ) write( ev, lLam0     ,     Lam0Name , writeEvent ); 
  if ( writeLbToLam0 ) write( ev, lLbToLam0 , LbToLam0Name , writeEvent ); 
  return;
}


void lbWriteSpecificDecay::fill( edm::Event& ev,
                                  const edm::EventSetup& es ) {
    // clean up {{{
    lFull    .clear();
    lJPsi    .clear();
    lLam0    .clear();
    lLbToLam0.clear();
    jPsiOMap .clear();
    pvRefMap .clear();
    ccRefMap .clear();
    // clean up end }}}

    // get magnetic field
    edm::ESHandle<MagneticField> magneticField;
    es.get<IdealMagneticFieldRecord>().get( magneticField );

    std::map<const BPHRecoCandidate*,const reco::Vertex*> oniaVtxMap;
    // find muon and get MuMu Onia full list, use them to decide primary vertex {{{
    // get object collections
    // collections are got through "BPHTokenWrapper" interface to allow
    // uniform access in different CMSSW versions
    
    edm::Handle< vector<reco::Vertex> > pVertices;
    int npv = 0;
    if ( usePV )
    {
        pVertexToken.get( ev, pVertices );
        npv = pVertices->size();
    }

    int nrc = 0;

    // get reco::PFCandidate collection (in full AOD )
    edm::Handle< vector<reco::PFCandidate> > pfCands;
    if ( usePF )
    {
        pfCandsToken.get( ev, pfCands );
        nrc = pfCands->size();
    }

    // get pat::PackedCandidate collection (in MiniAOD)
    // pat::PackedCandidate is not defined in CMSSW_5XY, so a
    // typedef (BPHTrackReference::candidate) is used, actually referring
    // to pat::PackedCandidate only for CMSSW versions where it's defined
    edm::Handle< vector<BPHTrackReference::candidate> > pcCands;
    if ( usePC )
    {
        pcCandsToken.get( ev, pcCands );
        nrc = pcCands->size();
    }

    // get pat::GenericParticle collection (in skimmed data)
    edm::Handle< vector<pat::GenericParticle> > gpCands;
    if ( useGP )
    {
        gpCandsToken.get( ev, gpCands );
        nrc = gpCands->size();
    }
    if( nrc > 0 ) ;

    // get pat::Muon collection (in full AOD and MiniAOD)
    edm::Handle<pat::MuonCollection> patMuon;
    if ( usePM )
    {
        patMuonToken.get( ev, patMuon );
    }
    // get muons from pat::CompositeCandidate objects describing onia;
    // muons from all composite objects are copied to an unique std::vector
    vector<const reco::Candidate*> muDaugs;
    set<const pat::Muon*> muonSet;
    typedef multimap<const reco::Candidate*,
            const pat::CompositeCandidate*> mu_cc_map;
    mu_cc_map muCCMap;
    if ( useCC )
    {
        // take out data
        edm::Handle< vector<pat::CompositeCandidate> > ccCands;
        ccCandsToken.get( ev, ccCands );
        int n = ccCands->size();
        muDaugs.clear();
        muDaugs.reserve( n );
        muonSet.clear();
        // add muon track from data{{{
        set<const pat::Muon*>::const_iterator iter;
        set<const pat::Muon*>::const_iterator iend;
        int i;
        for ( i = 0; i < n; ++i )
        {
            const pat::CompositeCandidate& cc = ccCands->at( i );
            int j;
            int m = cc.numberOfDaughters();
            for ( j = 0; j < m; ++j )
            {
                const reco::Candidate* dp = cc.daughter( j );
                const pat::Muon* mp = dynamic_cast<const pat::Muon*>( dp );
                iter = muonSet.begin();
                iend = muonSet.end();
                // if there is "daughter(j)" in "mp" not listed in "muonSet", push it into "muonSet"
                // while loop: check for the tracks are the same or not in "muonSet" and "mp"
                bool add = ( mp != 0 ) && ( muonSet.find( mp ) == iend );
                while ( add && ( iter != iend ) )
                {
                    if ( BPHRecoBuilder::sameTrack( mp, *iter++, 1.0e-5 ) ) add = false;
                }
                if ( add ) muonSet.insert( mp );
                // associate muon to the CompositeCandidate containing it
                muCCMap.insert( pair<const reco::Candidate*,
                                const pat::CompositeCandidate*>( dp, &cc ) );
            }
        }// add muon track from data end}}}

        iter = muonSet.begin();
        iend = muonSet.end();
        while ( iter != iend ) muDaugs.push_back( *iter++ );
    }

    map< recoType, map<parType,double> >::const_iterator rIter = parMap.begin();
    map< recoType, map<parType,double> >::const_iterator rIend = parMap.end  ();

    // reconstruct quarkonia

    BPHOniaToMuMuBuilder* onia = 0;
    if (recoOnia )
    {
        if ( usePM )
            onia = new BPHOniaToMuMuBuilder( es,
                    BPHRecoBuilder::createCollection( patMuon, "cfmign" ),
                    BPHRecoBuilder::createCollection( patMuon, "cfmign" ) );
        else if ( useCC ) 
            onia = new BPHOniaToMuMuBuilder( es,
                    BPHRecoBuilder::createCollection( muDaugs, "cfmig" ),
                    BPHRecoBuilder::createCollection( muDaugs, "cfmig" ) );
    }


    // Check mapList and pass particle type && cut value into onia{{{


    // Check for mapList and assign a particle type into "onia"
    if ( onia != 0 )
    {
        while ( rIter != rIend )
        {
            const map< recoType, map<parType,double> >::value_type& rEntry = *rIter++;
            recoType                   rType = rEntry.first;
            const map<parType,double>& pMap  = rEntry.second;//asdf
            BPHOniaToMuMuBuilder::oniaType type;
            switch( rType )
            {
                case Psi1: type = BPHOniaToMuMuBuilder::Psi1; break;
                case Psi2: type = BPHOniaToMuMuBuilder::Psi2; break;
                //case Ups : type = BPHOniaToMuMuBuilder::Ups ; break;
                //case Ups1: type = BPHOniaToMuMuBuilder::Ups1; break;
                //case Ups2: type = BPHOniaToMuMuBuilder::Ups2; break;
                //case Ups3: type = BPHOniaToMuMuBuilder::Ups3; break;
            default:
                continue;
            }
            map<parType,double>::const_iterator pIter = pMap.begin();
            map<parType,double>::const_iterator pIend = pMap.end();
 
            // pass particle typename and cutvalue into "onia"
            while ( pIter != pIend )
            {
                const map<parType,double>::value_type& pEntry = *pIter++;
                parType id = pEntry.first;
                double  pv = pEntry.second;
                switch( id )
                {
                    case ptMin      : onia->setPtMin  ( type, pv ); break;
                    case etaMax     : onia->setEtaMax ( type, pv ); break;
                    case massMin    : onia->setMassMin( type, pv ); break;
                    case massMax    : onia->setMassMax( type, pv ); break;
                    case probMin    : onia->setProbMin( type, pv ); break;
                    case constrMass : onia->setConstr ( type, pv, onia->getConstrSigma( type )); break;
                    case constrSigma: onia->setConstr ( type, onia->getConstrMass ( type ), pv ); break;
                    default:
                        break;
                }
            }
        }
        lFull = onia->build();
    }// Check mapList and pass particle type && cut value into onia end}}}

    // associate onia to primary vertex
    // (take them into map)

    int iFull;
    int nFull = lFull.size();

    // Vertex information{{{
    typedef mu_cc_map::const_iterator mu_cc_iter;
    for ( iFull = 0; iFull < nFull; ++iFull )
    {

        const reco::Vertex* pVtx = 0;
        int pvId = 0;
        const BPHRecoCandidate* ptr = lFull[iFull].get();
        const std::vector<const reco::Candidate*>& daugs = ptr->daughters();

        // try to recover primary vertex association in skim data:
        // get the CompositeCandidate containing both muons
        pair<mu_cc_iter,mu_cc_iter> cc0 = muCCMap.equal_range(
                                              ptr->originalReco( daugs[0] ) );
        pair<mu_cc_iter,mu_cc_iter> cc1 = muCCMap.equal_range(
                                              ptr->originalReco( daugs[1] ) );
        mu_cc_iter iter0 = cc0.first;
        mu_cc_iter iend0 = cc0.second;
        mu_cc_iter iter1 = cc1.first;
        mu_cc_iter iend1 = cc1.second;
        // if you use "useCC", that will fill the map "muCCMap",
        // if(useCC), you can find PV in this way {{{
        while ( ( iter0 != iend0 ) && ( pVtx == 0 )  )
        {
            // if useCC, you can find pV in this loop
            const pat::CompositeCandidate* ccp = iter0++->second;
            while ( iter1 != iend1 )
            {
                if ( ccp != iter1++->second ) continue;
                // If iter0->second == iter1->second

                // get the vertex information, then there is no any loop by this while
                pVtx = ccp->userData<reco::Vertex>( "PVwithmuons" );
                const reco::Vertex* sVtx = 0;
                const reco::Vertex::Point& pPos = pVtx->position();
                float dMin = 999999.;
                int ipv;
                for ( ipv = 0; ipv < npv; ++ipv )
                {
                    const reco::Vertex* tVtx = &pVertices->at( ipv );
                    const reco::Vertex::Point& tPos = tVtx->position();
                    float dist = pow( pPos.x() - tPos.x(), 2 ) +
                                 pow( pPos.y() - tPos.y(), 2 ) +
                                 pow( pPos.z() - tPos.z(), 2 );
                    if ( dist < dMin )
                    {
                        dMin = dist;
                        sVtx = tVtx;
                        pvId = ipv;
                    }
                }
                // find the pointer for smallest distance primaryV to vertex
                pVtx = sVtx;
                break;
            }
        } // if(useCC), you can find PV in this way end }}}
        // if not found, as ofr other type of inut data,
        // try to get the nearest primary vertex in z direction
        // if not useCC {{{
        if ( pVtx == 0 )
        {
            const reco::Vertex::Point& sVtp = ptr->vertex().position();
            GlobalPoint  cPos( sVtp.x(), sVtp.y(), sVtp.z() );
            const pat::CompositeCandidate& sCC = ptr->composite();
            GlobalVector cDir( sCC.px(), sCC.py(), sCC.pz() );
            GlobalPoint  bPos( 0.0, 0.0, 0.0 );
            GlobalVector bDir( 0.0, 0.0, 1.0 );
            TwoTrackMinimumDistance ttmd;
            bool state = ttmd.calculate( GlobalTrajectoryParameters( cPos, cDir,
                                         TrackCharge( 0 ), &( *magneticField ) ),
                                         GlobalTrajectoryParameters( bPos, bDir,
                                         TrackCharge( 0 ), &( *magneticField ) ) );
            float minDz = 999999.;
            float extrapZ = ( state ? ttmd.points().first.z() : -9e20 );
            int ipv;
            for ( ipv = 0; ipv < npv; ++ipv )
            {
                const reco::Vertex& tVtx = pVertices->at( ipv );
                float deltaZ = fabs( extrapZ - tVtx.position().z() ) ;
                if ( deltaZ < minDz )
                {
                    minDz = deltaZ;
                    pVtx = &tVtx;
                    pvId = ipv;
                }
            }
        }// if not useCC end}}}
        pvRefMap[ptr] = vertex_ref( pVertices, pvId );
        oniaVtxMap[ptr] = pVtx;


    }// vertex information end}}}


    // get JPsi subsample and associate JPsi candidate to original
    // generic onia candidate
    if ( nFull ) lJPsi = onia->getList( BPHOniaToMuMuBuilder::Psi1 );

    int nJPsi = lJPsi.size();
    delete onia; // chk pV end }}}

    if ( !nJPsi ) return;
    // Search for the map of lJPsi and lFull {{{
    if ( !nrc   ) return;

    int ij;
    int io;
    int nj = lJPsi.size();
    int no = lFull.size();
    for ( ij = 0; ij < nj; ++ij )
    {
        const BPHRecoCandidate* jp = lJPsi[ij].get();
        for ( io = 0; io < no; ++io )
        {
            const BPHRecoCandidate* oc = lFull[io].get();
            if ( ( jp->originalReco( jp->getDaug( "MuPos" ) ) ==
                    oc->originalReco( oc->getDaug( "MuPos" ) ) ) &&
                    ( jp->originalReco( jp->getDaug( "MuNeg" ) ) ==
                      oc->originalReco( oc->getDaug( "MuNeg" ) ) ) )
            {
                jPsiOMap[jp] = oc;
                break;
            }
        }
    } 
    // Search for the map of lJPsi and lFull end}}}

    // Build Lam0 {{{
    BPHLambda0ToPPiBuilder* tktk = 0;

    if ( recoLam0 )
    {
        if      ( usePF ) tktk = new BPHLambda0ToPPiBuilder( es,
                    BPHRecoBuilder::createCollection( pfCands ),
                    BPHRecoBuilder::createCollection( pfCands )
                    );
        else if ( usePC ) tktk = new BPHLambda0ToPPiBuilder( es,
                    BPHRecoBuilder::createCollection( pcCands ),
                    BPHRecoBuilder::createCollection( pcCands )
                    );
        else if ( useGP ) tktk = new BPHLambda0ToPPiBuilder( es,
                    BPHRecoBuilder::createCollection( gpCands ),
                    BPHRecoBuilder::createCollection( gpCands )
                    );
    }
    // Set cut value 
    // use "Name"(recorded in enum) to find particle, then there is a subMap
    // The subMap is list of the cuts used in  the particle
    if ( tktk != 0 )
    {
        rIter = parMap.find( Lam0 );

        // if find something
        if ( rIter != rIend )
        {
            const map<parType,double>& _parMap = rIter->second;
            map<parType,double>::const_iterator _parIter = _parMap.begin();
            map<parType,double>::const_iterator _parIend = _parMap.end();
            while ( _parIter != _parIend )
            {
                // set cut value by switch
                const map<parType,double>::value_type& _parEntry = *_parIter++;
                parType _parId      = _parEntry.first;
                double  _parValue   = _parEntry.second;
                switch( _parId )
                {
                    case ptMin          : tktk->setPtMin       ( _parValue ); break;
                    case etaMax         : tktk->setEtaMax      ( _parValue ); break;
                    case probMin        : tktk->setProbMin     ( _parValue ); break;
                    case writeCandidate : writeLam0 =          ( _parValue > 0 ); break;
                    default:
                        break;
                }
            }
        }
        lLam0 = tktk->build();
        delete   tktk;


    } // set cut value end  
    unsigned nLam0 = lLam0.size();
    // Build Lam0 end }}}

    // Build and dump Lb->Jpsi+Lam0 {{{
    if ( nLam0 && recoLbToLam0 )
    {
        BPHLambda0_bToJPsiLambda0Builder* _lb = new BPHLambda0_bToJPsiLambda0Builder( es, lJPsi, lLam0);
        // Set cut value 
        // use "Name"(recorded in enum) to find particle, then there is a subMap
        // The subMap is list of the cuts used in  the particle
        rIter = parMap.find( LbToLam0 );

        // if find something
        if ( rIter != rIend )
        {
            const map<parType,double>& _parMap = rIter->second;
            map<parType,double>::const_iterator _parIter = _parMap.begin();
            map<parType,double>::const_iterator _parIend = _parMap.end();
            while ( _parIter != _parIend )
            {
                // set cut value by switch
                const map<parType,double>::value_type& _parEntry = *_parIter++;
                parType _parId      = _parEntry.first;
                double  _parValue   = _parEntry.second;
                switch( _parId )
                {
                    case mPsiMin        : _lb->setJPsiMassMin ( _parValue ); break;
                    case mPsiMax        : _lb->setJPsiMassMax ( _parValue ); break;
                    case mLam0Min       : _lb->setLam0MassMin ( _parValue ); break;
                    case mLam0Max       : _lb->setLam0MassMax ( _parValue ); break;
                                          
                    case massMin        : _lb->setMassMin     ( _parValue ); break;
                    case massMax        : _lb->setMassMax     ( _parValue ); break;
                    case probMin        : _lb->setProbMin     ( _parValue ); break;
    
                    case mFitMin        : _lb->setMassFitMin  ( _parValue ); break;
                    case mFitMax        : _lb->setMassFitMax  ( _parValue ); break;
                    case constrMJPsi    : _lb->setConstr        ( _parValue ); break;
                    case writeCandidate : writeLbToLam0 =     ( _parValue > 0 ); break;
                    default:
                        break;
                }
            }
        }
        lLbToLam0 = _lb->build();
        delete   _lb;
        // set cut value end  
    }
    // Build LbToLam0 end }}}
}


void lbWriteSpecificDecay::endJob() {
  return;
}


void lbWriteSpecificDecay::setRecoParameters( const edm::ParameterSet& ps ) {
    const string& name = ps.getParameter<string>( "name" );
    bool writeCandidate = ps.getParameter<bool>( "writeCandidate" );
    switch( rMap[name] )
    {
        case Onia     : recoOnia      = true; writeOnia       = writeCandidate; break;
        case Psi1     :
        case Psi2     : recoOnia      = true;                                   break;
        case Lam0     : recoLam0      = true; writeLam0       = writeCandidate; break;
        case LbToLam0 : recoLbToLam0  = true; writeLbToLam0   = writeCandidate; break;
    }


    map<string,parType>::const_iterator pIter = pMap.begin();
    map<string,parType>::const_iterator pIend = pMap.end();
    while ( pIter != pIend )
    {
        const map<string,parType>::value_type& entry = *pIter++;
        const string& pn = entry.first;
        parType       id = entry.second;
        double        pv = ps.getParameter<double>( pn );
        if ( pv > -1.0e35 ) edm::LogVerbatim( "Configuration" )
                    << "LbSpecificDecay::setRecoParameters: set " << pn
                    << " for " << name << " : "
                    << ( parMap[rMap[name]][id] = pv );
    }

    map<string,parType>::const_iterator fIter = fMap.begin();
    map<string,parType>::const_iterator fIend = fMap.end();
    while ( fIter != fIend )
    {
        const map<string,parType>::value_type& entry = *fIter++;
        const string& fn = entry.first;
        parType       id = entry.second;
        edm::LogVerbatim( "Configuration" )
                << "LbSpecificDecay::setRecoParameters: set " << fn
                << " for " << name << " : "
                << ( parMap[rMap[name]][id] =
                         ( ps.getParameter<bool>( fn ) ? 1 : -1 ) );
    }

}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( lbWriteSpecificDecay );
