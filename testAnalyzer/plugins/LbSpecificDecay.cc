#include "BPHAnalysis/testAnalyzer/plugins/LbSpecificDecay.h"

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
#include "BPHAnalysis/RecoDecay/interface/BPHMultiSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHMassSelect.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHChi2Select.h"

#include "BPHAnalysis/SpecificDecay/interface/BPHOniaToMuMuBuilder.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHLambda0ToPPiBuilder.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHKx0ToKPiBuilder.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHLambda0_bToJPsiLambda0Builder.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHTkTkBuilder.h"
#include "BPHAnalysis/SpecificDecay/interface/BPHLambda0_bToJPsiTkTkBuilder.h"

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
// xyz = ps.getParameter< string >( "xyx" )

LbSpecificDecay::LbSpecificDecay( const edm::ParameterSet& ps )
{

    //  Check if there is the label used in reading data or not.
    usePV = ( SET_PAR( string, pVertexLabel, ps ) != "" );
    usePM = ( SET_PAR( string, patMuonLabel, ps ) != "" );
    useCC = ( SET_PAR( string, ccCandsLabel, ps ) != "" );
    usePF = ( SET_PAR( string, pfCandsLabel, ps ) != "" );
    usePC = ( SET_PAR( string, pcCandsLabel, ps ) != "" );
    useGP = ( SET_PAR( string, gpCandsLabel, ps ) != "" );

    // Set the label used in output product.
    SET_PAR( string,     oniaName, ps );
    SET_PAR( string,     Lam0Name, ps );
    SET_PAR( string,     TkTkName, ps );
    SET_PAR( string, LbToLam0Name, ps );
    SET_PAR( string, LbToTkTkName, ps );
    SET_PAR( bool  ,writeMomentum, ps );
    SET_PAR( bool  ,writeVertex  , ps );

    rMap["Onia"    ] = Onia;
    rMap["JPsi"    ] = Psi1;
    rMap["Psi2"    ] = Psi2;
    rMap["Lam0"    ] = Lam0;
    rMap["TkTk"    ] = TkTk;
    rMap["LbToLam0"] = LbToLam0;
    rMap["LbToTkTk"] = LbToTkTk;

    pMap["ptMin"      ] = ptMin;
    pMap["etaMax"     ] = etaMax;
    pMap["mPsiMin"    ] = mPsiMin;
    pMap["mPsiMax"    ] = mPsiMax;
    pMap["mLam0Min"   ] = mLam0Min;
    pMap["mLam0Max"   ] = mLam0Max;
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
    recoLam0      = writeLam0     =  
    recoTkTk      = writeTkTk     =   true;
    recoLbToLam0  = writeLbToLam0 =
    recoLbToTkTk  = writeLbToTkTk =  false;

    writeOnia = true;
    const vector<edm::ParameterSet> recoSelect =
        ps.getParameter< vector<edm::ParameterSet> >( "recoSelect" );
    int iSel;
    int nSel = recoSelect.size();
    for ( iSel = 0; iSel != nSel; ++iSel ) setRecoParameters( recoSelect[iSel] );
    if ( !recoOnia ) writeOnia = false;

    if (  recoLbToLam0 )  recoOnia =  recoLam0   = true;
    if (  recoLbToTkTk )  recoOnia =  recoTkTk   = true;
    if ( writeLbToLam0 ) writeOnia = writeLam0   = true;
    if ( writeLbToTkTk ) writeOnia = writeTkTk   = true;

    // Get data by label/token
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

    _dmassMAP[ "MuPos"  ] = BPHParticleMasses::muonMass;
    _sigmaMAP[ "MuPos"  ] = BPHParticleMasses::muonMSigma;
    _dmassMAP[ "MuNeg"  ] = BPHParticleMasses::muonMass;
    _sigmaMAP[ "MuNeg"  ] = BPHParticleMasses::muonMSigma;
    _dmassMAP[ "Pion"   ] = BPHParticleMasses::pionMass;
    _sigmaMAP[ "Pion"   ] = BPHParticleMasses::pionMSigma;
    _dmassMAP[ "Kaon"   ] = BPHParticleMasses::kaonMass;
    _sigmaMAP[ "Kaon"   ] = BPHParticleMasses::kaonMSigma;
    _dmassMAP[ "Proton" ] = BPHParticleMasses::protonMass;
    _sigmaMAP[ "Proton" ] = BPHParticleMasses::protonMSigma;
    eventNo = 0;
}


LbSpecificDecay::~LbSpecificDecay()
{
    printf("----analysis ended, total events are: %u\n", eventNo);
}

void LbSpecificDecay::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions )
{
    edm::ParameterSetDescription desc;
    desc.add<string>( "pVertexLabel", "" );
    desc.add<string>( "patMuonLabel", "" );
    desc.add<string>( "ccCandsLabel", "" );
    desc.add<string>( "pfCandsLabel", "" );
    desc.add<string>( "pcCandsLabel", "" );
    desc.add<string>( "gpCandsLabel", "" );
    desc.add<string>(     "oniaName",   "oniaCand" );
    desc.add<string>(     "Lam0Name",   "Lam0Cand" );
    desc.add<string>(     "TkTkName",   "TkTkCand" );
    desc.add<string>( "LbToLam0Name", "LbToLam0Fitted" );
    desc.add<string>( "LbToTkTkName", "LbToTkTkFitted" );

    desc.add<bool>  ( "writeVertex"  , true );
    desc.add<bool>  ( "writeMomentum", true );
    edm::ParameterSetDescription dpar;
    dpar.add<string>(           "name" );
    dpar.add<double>(          "ptMin", -2.0e35 );
    dpar.add<double>(         "etaMax", -2.0e35 );
    dpar.add<double>(        "mPsiMin", -2.0e35 );
    dpar.add<double>(        "mPsiMax", -2.0e35 );
    dpar.add<double>(       "mLam0Min", -2.0e35 );
    dpar.add<double>(       "mLam0Max", -2.0e35 );
    dpar.add<double>(        "massMin", -2.0e35 );
    dpar.add<double>(        "massMax", -2.0e35 );
    dpar.add<double>(        "probMin", -2.0e35 );
    dpar.add<double>(     "massFitMin", -2.0e35 );
    dpar.add<double>(     "massFitMax", -2.0e35 );
    dpar.add<double>(     "constrMass", -2.0e35 );
    dpar.add<double>(    "constrSigma", -2.0e35 );
    dpar.add<  bool>(    "constrMJPsi",    true );
    dpar.add<  bool>( "writeCandidate",    true );

    vector<edm::ParameterSet> rpar;
    desc.addVPSet( "recoSelect", dpar, rpar );
    descriptions.add( "LbSpecificDecay", desc );
    return;
}




void LbSpecificDecay::beginJob()
{
    jpsiTree = fs->make<TTree>("jpsiTree","jpsiTree");
    lam0Tree = fs->make<TTree>("lam0Tree","lam0Tree");
    lambTree = fs->make<TTree>("LbToLam0Tree","LbToLam0Tree");
    lbv1Tree = fs->make<TTree>("LbToTkTree","LbToTkTree");
    lbv2Tree = fs->make<TTree>("LbTotkTree","LbTotkTree");


    jpsiBr.RegTree(jpsiTree);
    lam0Br.RegTree(lam0Tree);
    lambBr.RegTree(lambTree);
    lbv1Br.RegTree(lbv1Tree);
    lbv2Br.RegTree(lbv2Tree);

    jpsiBr.RegTree(jpsiTree);

    return;
}


void LbSpecificDecay::analyze( const edm::Event& ev,
                               const edm::EventSetup& es )
{
    fill( ev, es );
    //fillTree();
    fillTree(es);
    ++eventNo;
    return;
}


void LbSpecificDecay::fill( const edm::Event& ev,
                            const edm::EventSetup& es )
{

    // clean up {{{
    lFull    .clear();
    lJPsi    .clear();
    lLam0    .clear();
    lLbToLam0.clear();
    lLbToTkTk.clear();
    jPsiOMap .clear();
    pvRefMap .clear();
    ccRefMap .clear();
    // clean up end }}}

    // get magnetic field
    edm::ESHandle<MagneticField> magneticField;
    es.get<IdealMagneticFieldRecord>().get( magneticField );

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
                //case Pmm : type = BPHOniaToMuMuBuilder::Phi ; break;
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
    map<const BPHRecoCandidate*,const reco::Vertex*> oniaVtxMap;

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
    BPHLambda0ToPPiBuilder* lam0 = 0;

    if ( recoLam0 )
    {
        if      ( usePF ) lam0 = new BPHLambda0ToPPiBuilder( es,
                    BPHRecoBuilder::createCollection( pfCands ),
                    BPHRecoBuilder::createCollection( pfCands ) );
        else if ( usePC ) lam0 = new BPHLambda0ToPPiBuilder( es,
                    BPHRecoBuilder::createCollection( pcCands ),
                    BPHRecoBuilder::createCollection( pcCands ) );
        else if ( useGP ) lam0 = new BPHLambda0ToPPiBuilder( es,
                    BPHRecoBuilder::createCollection( gpCands ),
                    BPHRecoBuilder::createCollection( gpCands ) );
    }
    // Set cut value 
    // use "Name"(recorded in enum) to find particle, then there is a subMap
    // The subMap is list of the cuts used in  the particle
    if ( lam0 != 0 )
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
                    case ptMin          : lam0->setPtMin       ( _parValue ); break;
                    case etaMax         : lam0->setEtaMax      ( _parValue ); break;
                    case massMin        : lam0->setMassMin     ( _parValue ); break;
                    case massMax        : lam0->setMassMax     ( _parValue ); break;
                    case probMin        : lam0->setProbMin     ( _parValue ); break;
                    case constrMass     : lam0->setConstr      ( _parValue, lam0->getConstrSigma() ); break;
                    case constrSigma    : lam0->setConstr      ( lam0->getConstrMass(), _parValue  ); break;
                    case writeCandidate : writeLam0 =          ( _parValue > 0 ); break;
                    default:
                        break;
                }
            }
        }
        lLam0 = lam0->build();
        delete   lam0;


    } // set cut value end  
    // Build Lam0 end }}}
    // Build TkTk {{{
    BPHTkTkBuilder* tktk = 0;
    std::string pName = "Proton";
    std::string nName = "Kaon";
    double pMass = BPHParticleMasses::protonMass;
    double nMass = BPHParticleMasses::  kaonMass;
    double pSigma= BPHParticleMasses::protonMSigma;
    double nSigma= BPHParticleMasses::  kaonMSigma;

    if ( recoTkTk )
    {
        if      ( usePF ) tktk = new BPHTkTkBuilder( es,
                    BPHRecoBuilder::createCollection( pfCands ),
                    BPHRecoBuilder::createCollection( pfCands ),
                    pName, pMass, pSigma,
                    nName, nMass, nSigma
                    );
        else if ( usePC ) tktk = new BPHTkTkBuilder( es,
                    BPHRecoBuilder::createCollection( pcCands ),
                    BPHRecoBuilder::createCollection( pcCands ),
                    pName, pMass, pSigma,
                    nName, nMass, nSigma
                    );
        else if ( useGP ) tktk = new BPHTkTkBuilder( es,
                    BPHRecoBuilder::createCollection( gpCands ),
                    BPHRecoBuilder::createCollection( gpCands ),
                    pName, pMass, pSigma,
                    nName, nMass, nSigma
                    );
    }
    // Set cut value 
    // use "Name"(recorded in enum) to find particle, then there is a subMap
    // The subMap is list of the cuts used in  the particle
    if ( tktk != 0 )
    {
        rIter = parMap.find( TkTk );

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
                    case writeCandidate : writeTkTk =          ( _parValue > 0 ); break;
                    default:
                        break;
                }
            }
        }
        lTkTk = tktk->build();
        delete   tktk;


    } // set cut value end  
    unsigned nTkTk = lTkTk.size();
    // Build TkTk end }}}

    // build and dump Lb->Jpsi+Lam0{{{

    
    int nLam0 = lLam0.size();


    // Lambda0 is built, start to build Lambda0_b->J/Psi+Lambda0
    if ( recoLam0 && nLam0 )
    {
        BPHLambda0_bToJPsiLambda0Builder* _lb = new BPHLambda0_bToJPsiLambda0Builder( es, lJPsi, lLam0 );
        
        // set cut value 
        rIter = parMap.find( LbToLam0 );
        if ( rIter != rIend )
        {
            const map<parType,double>& _parMap = rIter->second;
            map<parType,double>::const_iterator _parIter = _parMap.begin();
            map<parType,double>::const_iterator _parIend = _parMap.end  ();
            while ( _parIter != _parIend )
            {
                const map<parType,double>::value_type& _parEntry = *_parIter++;
                parType  _parId      = _parEntry.first;
                double   _parValue   = _parEntry.second;

                switch( _parId )
                {
                case mPsiMin        : _lb->setJPsiMassMin   ( _parValue ); break;
                case mPsiMax        : _lb->setJPsiMassMax   ( _parValue ); break;
                case mLam0Min       : _lb->setLam0MassMin   ( _parValue ); break;
                case mLam0Max       : _lb->setLam0MassMax   ( _parValue ); break;
                case massMin        : _lb->setMassMin       ( _parValue ); break;
                case massMax        : _lb->setMassMax       ( _parValue ); break;
                case probMin        : _lb->setProbMin       ( _parValue ); break;
                case mFitMin        : _lb->setMassFitMin    ( _parValue ); break;
                case mFitMax        : _lb->setMassFitMax    ( _parValue ); break;
                case constrMJPsi    : _lb->setConstr        ( _parValue > 0); break;
                case writeCandidate : writeLbToLam0 =       ( _parValue > 0); break;
                default:
                    break;
                }
            }
        }
        lLbToLam0 = _lb->build();
        delete _lb;
        // set cut value end 
    }// build lambda b end }}}

    // Build and dump Lb->Jpsi+TkTk {{{
    if ( nJPsi && nTkTk && recoLbToTkTk )
    {
        BPHLambda0_bToJPsiTkTkBuilder* _lb = 0;
        if ( usePF ) _lb = new BPHLambda0_bToJPsiTkTkBuilder( es, lJPsi, lTkTk );
        else if ( usePC ) _lb = new BPHLambda0_bToJPsiTkTkBuilder( es, lJPsi, lTkTk );
        else if ( useGP ) _lb = new BPHLambda0_bToJPsiTkTkBuilder( es, lJPsi, lTkTk );
        // Set cut value 
        // use "Name"(recorded in enum) to find particle, then there is a subMap
        // The subMap is list of the cuts used in  the particle
        if ( _lb != 0 )
        {
            rIter = parMap.find( LbToTkTk );

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
                        case massMin        : _lb->setMassMin     ( _parValue ); break;
                        case massMax        : _lb->setMassMax     ( _parValue ); break;
                        case probMin        : _lb->setProbMin     ( _parValue ); break;
    
                        case mFitMin        : _lb->setMassFitMin  ( _parValue ); break;
                        case mFitMax        : _lb->setMassFitMax  ( _parValue ); break;
                        case writeCandidate : writeLbToTkTk =     ( _parValue > 0 ); break;
                        default:
                            break;
                    }
                }
            }
            lLbToTkTk = _lb->build();
            delete   _lb;
        }
        // set cut value end  

        


    }
    // Build LbToTkTk end }}}
    // Build and dump ant Lb->Jpsi+TkTk {{{
    if ( nJPsi && nTkTk && recoLbToTkTk )
    {
        BPHLambda0_bToJPsiTkTkBuilder* _lb = 0;
        // regard Proton track as Kaon mass and sigma.
        // That is to say, change the particle track and refit it.
        // finally, the "anti-Proton" with the name "Kaon"
        // and "anti-kaon" with the name "Proton"
        if ( usePF ) _lb = new BPHLambda0_bToJPsiTkTkBuilder( es, lJPsi, lTkTk);
        else if ( usePC ) _lb = new BPHLambda0_bToJPsiTkTkBuilder( es, lJPsi, lTkTk);
        else if ( useGP ) _lb = new BPHLambda0_bToJPsiTkTkBuilder( es, lJPsi, lTkTk);
        // Set cut value 
        // use "Name"(recorded in enum) to find particle, then there is a subMap
        // The subMap is list of the cuts used in  the particle
        if ( _lb != 0 )
        {
            rIter = parMap.find( LbToTkTk );

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
                        case massMin        : _lb->setMassMin     ( _parValue ); break;
                        case massMax        : _lb->setMassMax     ( _parValue ); break;
                        case probMin        : _lb->setProbMin     ( _parValue ); break;
    
                        case mFitMin        : _lb->setMassFitMin  ( _parValue ); break;
                        case mFitMax        : _lb->setMassFitMax  ( _parValue ); break;
                        case writeCandidate : writeLbToTkTk =     ( _parValue > 0 ); break;
                        default:
                            break;
                    }
                }
            }
            laLbToTkTk = _lb->build();
            delete   _lb;
        }
        // set cut value end  

        


    }
    // Build LbToTkTk end }}}

    return;
}

void LbSpecificDecay::fillTree( const edm::EventSetup& es )
{
    // apply kinematic fit, and fill the tree
    
    // clean up
    memset(&jpsiBr,0x00,sizeof(jpsiBr));
    memset(&lam0Br,0x00,sizeof(lam0Br));
    memset(&lambBr,0x00,sizeof(lambBr));
    memset(&lbv1Br,0x00,sizeof(lbv1Br));
    memset(&lbv2Br,0x00,sizeof(lbv2Br));

    std::cout << "writeOnia\n";

    if (writeOnia) // {{{
    {
        unsigned int nJPsi = lJPsi.size();
        for ( unsigned int iJPsi = 0; iJPsi < nJPsi; ++iJPsi )
        {
            memset(&jpsiBr,0x00,sizeof(jpsiBr));

            const BPHRecoCandidate* cand = lJPsi[iJPsi].get();

            //const reco::Vertex& vx = cand->vertex();
            //TVector3 bs_to_vtx(vx.position().X()-BeamSpotVtx.position().X(),
            //                   vx.position().Y()-BeamSpotVtx.position().Y(),
            //                   vx.position().Z()-BeamSpotVtx.position().Z());

            jpsiBr.momentum.mass    = cand->composite().mass();
            jpsiBr.momentum.pt      = cand->composite().pt();
            jpsiBr.momentum.eta     = cand->composite().eta();
            jpsiBr.momentum.phi     = cand->composite().phi();
            jpsiBr.eventNumber      = eventNo;

            jpsiTree->Fill();
        }
    } // writeOnia end }}}
    std::cout << "writeLam0\n";
    if (writeLam0) // {{{
    {
        unsigned nLam0 = lLam0.size();
        // full the Lambda tree
        for ( unsigned int iLam = 0; iLam < nLam0; ++iLam )
        {
            memset(&lam0Br,0x00,sizeof(lam0Br));
            const BPHRecoCandidate* cand = lLam0[iLam].get();
            const std::vector<std::string> lam0Name = { "Proton", "Pion" };
            const SecRecoResult lam0Reco = secondaryReconstruction(cand, &lam0Name);

            const reco::Candidate* proton = cand->getDaug("Proton");
            const reco::Candidate* pion   = cand->getDaug("Pion"  );
    
            lam0Br.momentum.mass    = cand->composite().mass();
            lam0Br.momentum.pt      = cand->composite().pt();
            lam0Br.momentum.eta     = cand->composite().eta();
            lam0Br.momentum.phi     = cand->composite().phi();
    
            if ( lam0Reco.isValid() )
            {
                lam0Br.refitMom.mass    = lam0Reco.mass();
                lam0Br.refitMom.pt      = lam0Reco.momentum()->transverse();  
                lam0Br.refitMom.eta     = lam0Reco.momentum()->eta(); 
                lam0Br.refitMom.phi     = lam0Reco.momentum()->phi();

                lam0Br.refitPos.x       = lam0Reco.vertex()->position().X(); 
                lam0Br.refitPos.y       = lam0Reco.vertex()->position().Y(); 
                lam0Br.refitPos.z       = lam0Reco.vertex()->position().Z(); 
                lam0Br.refitPos.vtxprob = TMath::Prob( lam0Reco.vertex()->chi2(), lam0Reco.vertex()->ndof() );
            }
            lam0Br.proton.px    = proton->px();
            lam0Br.proton.py    = proton->py();
            lam0Br.proton.pz    = proton->pz();
            lam0Br.proton.charge= proton->charge();
            lam0Br.pion  .px    = pion  ->px();
            lam0Br.pion  .py    = pion  ->py();
            lam0Br.pion  .pz    = pion  ->pz();
            lam0Br.pion  .charge= pion  ->charge();
            lam0Br.eventNumber  = eventNo;

            lam0Tree->Fill();
        }
    }// writeLam0 end }}}
    std::cout << "writeLbL0\n";
    // LbToLam0 fill {{{
    unsigned int nLamb = lLbToLam0.size();
    for ( unsigned int iLamb = 0; iLamb < nLamb; ++iLamb )
    {
        memset(&lambBr,0x00,sizeof(lambBr));
        const BPHRecoCandidate* cand = lLbToLam0[iLamb].get();
        const BPHRecoCandidate* cand_jpsi = cand->getComp( "JPsi"   ).get();
        const BPHRecoCandidate* cand_lam0 = cand->getComp( "Lam0"   ).get();
        const std::vector<std::string> lambName = { "JPsi/MuPos", "JPsi/MuNeg", "Lam0/Proton", "Lam0/Pion" };
        ParticleMass jpsi_mass = 3.096916;
        TwoTrackMassKinematicConstraint *jpsi_const = new TwoTrackMassKinematicConstraint(jpsi_mass);
        const SecRecoResult lambReco = secondaryReconstruction(cand, &lambName, jpsi_const);
        const std::vector<std::string> lam0Name = { "Proton", "Pion" };
        const SecRecoResult lam0Reco = secondaryReconstruction(cand_lam0, &lam0Name);
        const reco::Vertex* PVptr = findPrimaryVertex(cand_jpsi);


        if ( PVptr )
        {
            lambBr.primaryV.x        = PVptr->position().X();
            lambBr.primaryV.y        = PVptr->position().Y();
            lambBr.primaryV.z        = PVptr->position().Z();
        }

        if ( lambReco.isValid() )             
        {
            lambBr.refitMom.mass       = lambReco.mass(); 
            lambBr.refitMom.pt         = lambReco.momentum()->transverse(); 
            lambBr.refitMom.eta        = lambReco.momentum()->eta(); 
            lambBr.refitMom.phi        = lambReco.momentum()->phi(); 
                          
            lambBr.refitPos.x          = lambReco.vertex()->position().X(); 
            lambBr.refitPos.y          = lambReco.vertex()->position().Y(); 
            lambBr.refitPos.z          = lambReco.vertex()->position().Z(); 
            lambBr.refitPos.vtxprob    = TMath::Prob( lambReco.vertex()->chi2(), lambReco.vertex()->ndof() );
        }


        lambBr.jpsiMom.mass      = cand_jpsi->composite().mass(); 
        lambBr.jpsiMom.pt        = cand_jpsi->composite().pt(); 
        lambBr.jpsiMom.eta       = cand_jpsi->composite().pt(); 
        lambBr.jpsiMom.phi       = cand_jpsi->composite().pt(); 

        lambBr.jpsiPos.x         = cand_jpsi->vertex().position().X(); 
        lambBr.jpsiPos.y         = cand_jpsi->vertex().position().Y(); 
        lambBr.jpsiPos.z         = cand_jpsi->vertex().position().Z(); 
        lambBr.jpsiPos.vtxprob   = TMath::Prob( cand_jpsi->vertex().chi2(), cand_jpsi->vertex().ndof() );


        if ( lam0Reco.isValid() )
        {
            lambBr.lam0Mom.mass      = lam0Reco.mass(); 
            lambBr.lam0Mom.pt        = lam0Reco.momentum()->transverse(); 
            lambBr.lam0Mom.eta       = lam0Reco.momentum()->eta(); 
            lambBr.lam0Mom.phi       = lam0Reco.momentum()->phi(); 
                                                                                                        
            lambBr.lam0Pos.x        = lam0Reco.vertex()->position().X(); 
            lambBr.lam0Pos.y        = lam0Reco.vertex()->position().Y(); 
            lambBr.lam0Pos.z        = lam0Reco.vertex()->position().Z(); 
            lambBr.lam0Pos.vtxprob  = TMath::Prob( lam0Reco.vertex()->chi2(), lam0Reco.vertex()->ndof() );
        }
        lambBr.eventNumber      = eventNo;
        // 0: lambda0_b reconstruction
        // 1: lambda0   reconstruction
        int tagInfo = 0;
        if ( lambReco.isValid() )   tagInfo += 1 << 0;
        if ( lam0Reco.isValid() )   tagInfo += 1 << 1;
        if ( PVptr )                tagInfo += 1 << 2;

        if ( tagInfo == (1<<3)-1 )
            lambTree->Fill();
    } // LbToLam0 fill end }}}

    std::cout << "writeLbTk\n";
    // LbToTkTk fill {{{
    if ( writeLbToTkTk )
    {
        nLamb = lLbToTkTk.size();
        for ( unsigned int iLamb = 0; iLamb < nLamb; ++iLamb )
        {
            memset(&lbv1Br,0x00,sizeof(lbv1Br));
            const BPHRecoCandidate* cand = lLbToTkTk[iLamb].get();
            const BPHRecoCandidate* cand_jpsi = cand->getComp( "JPsi" ).get();

            const std::vector<std::string> lambName = { "JPsi/MuPos", "JPsi/MuNeg", "TkTk/Proton", "TkTk/Kaon" };
            ParticleMass jpsi_mass = 3.096916;
            TwoTrackMassKinematicConstraint *jpsi_const = new TwoTrackMassKinematicConstraint(jpsi_mass);
            const SecRecoResult lambReco = secondaryReconstruction(cand, &lambName, jpsi_const);
            //const std::vector<std::string> penQName = { "JPsi/MuPos", "JPsi/MuNeg", "Proton"};
            //const SecRecoResult penQReco = secondaryReconstruction(cand, &penQName);
            const reco::Vertex* PVptr = findPrimaryVertex(cand_jpsi);

            if ( PVptr ) 
            {
                lbv1Br.primaryV.x = PVptr->position().X();
                lbv1Br.primaryV.y = PVptr->position().Y();
                lbv1Br.primaryV.z = PVptr->position().Z();
            }
            if ( cand->isValidFit() )
            {
                    const RefCountedKinematicParticle kinPart = cand->currentParticle();
                    const           KinematicState    kinStat = kinPart->currentState();
                lbv1Br.momentum.mass    = kinStat.mass();
                lbv1Br.momentum.pt      = kinStat.kinematicParameters().momentum().transverse();
                lbv1Br.momentum.eta     = kinStat.kinematicParameters().momentum().eta();
                lbv1Br.momentum.phi     = kinStat.kinematicParameters().momentum().phi();
                lbv1Br.position.x       = cand->currentDecayVertex()->position().x();
                lbv1Br.position.y       = cand->currentDecayVertex()->position().y();
                lbv1Br.position.z       = cand->currentDecayVertex()->position().z();
                //lbv1Br.refitPos.vtxprob = TMath::Prob( lambReco.vertex()->chi2(), lambReco.vertex()->ndof() );
                //}

            }

            if( lambReco.isValid() )
            {
                lbv1Br.refitMom.mass    = lambReco.mass();
                lbv1Br.refitMom.pt      = lambReco.momentum()->transverse();
                lbv1Br.refitMom.eta     = lambReco.momentum()->eta();
                lbv1Br.refitMom.phi     = lambReco.momentum()->phi();
                lbv1Br.refitPos.x       = lambReco.vertex()->position().X();
                lbv1Br.refitPos.y       = lambReco.vertex()->position().Y();
                lbv1Br.refitPos.z       = lambReco.vertex()->position().Z();
                lbv1Br.refitPos.vtxprob = TMath::Prob( lambReco.vertex()->chi2(), lambReco.vertex()->ndof() );
            }


            lbv1Br.jpsiMom.mass      = cand_jpsi->composite().mass();
            lbv1Br.jpsiMom.pt        = cand_jpsi->composite().pt();
            lbv1Br.jpsiMom.eta       = cand_jpsi->composite().eta();
            lbv1Br.jpsiMom.phi       = cand_jpsi->composite().phi();
            lbv1Br.jpsiPos.x         = cand_jpsi->vertex().position().X();
            lbv1Br.jpsiPos.y         = cand_jpsi->vertex().position().Y();
            lbv1Br.jpsiPos.z         = cand_jpsi->vertex().position().Z();
            lbv1Br.jpsiPos.vtxprob   = TMath::Prob( cand_jpsi->vertex().chi2(), cand_jpsi->vertex().ndof() );

            //if ( penQReco.isValid() )
            //{
            //    lbv1Br.penQMom.mass    = penQReco.mass();
            //    lbv1Br.penQMom.pt      = penQReco.momentum()->transverse();
            //    lbv1Br.penQMom.eta     = penQReco.momentum()->eta();
            //    lbv1Br.penQMom.phi     = penQReco.momentum()->phi();
            //    lbv1Br.penQPos.x       = penQReco.vertex()->position().X();
            //    lbv1Br.penQPos.y       = penQReco.vertex()->position().Y();
            //    lbv1Br.penQPos.z       = penQReco.vertex()->position().Z();
            //    lbv1Br.penQPos.vtxprob = TMath::Prob( penQReco.vertex()->chi2(), penQReco.vertex()->ndof() );
            //}

            const reco::Candidate* pTk = cand->getDaug("TkTk/Proton");
            const reco::Candidate* nTk = cand->getDaug("TkTk/Kaon"  );
            lbv1Br.pTk.px     = pTk->px(); 
            lbv1Br.pTk.py     = pTk->py(); 
            lbv1Br.pTk.pz     = pTk->pz(); 
            lbv1Br.pTk.charge = pTk->charge();
            lbv1Br.nTk.px     = nTk->px();
            lbv1Br.nTk.py     = nTk->py();
            lbv1Br.nTk.pz     = nTk->pz();
            lbv1Br.nTk.charge = nTk->charge();
            lbv1Br.eventNumber= eventNo;
            if ( pTk->charge() * nTk->charge() != -1 ) cout << "Lb is not neutral!\n";

            int tagInfo = 0;
            if ( PVptr )                tagInfo += 1 << 0;
            if ( lambReco.isValid() )   tagInfo += 1 << 1;
            //if ( penQReco.isValid() )   tagInfo += 1 << 1;

            if ( tagInfo == (1<<2)-1 )
                lbv1Tree->Fill();
        }
    } // LbToTkTk fill end }}}
//    // anti LbToTkTk fill {{{
//    if ( writeLbToTkTk )
//   {
//        nLamb = laLbToTkTk.size();
//        for ( unsigned int iLamb = 0; iLamb < nLamb; ++iLamb )
//        {
//            memset(&lbv2Br,0x00,sizeof(lbv2Br));
//            const BPHRecoCandidate* cand = laLbToTkTk[iLamb].get();
//            const BPHRecoCandidate* cand_jpsi = cand->getComp( "JPsi" ).get();
//
//            const std::vector<std::string> lambName = { "JPsi/MuPos", "JPsi/MuNeg", "Proton", "Kaon" };
//            ParticleMass jpsi_mass = 3.096916;
//            TwoTrackMassKinematicConstraint *jpsi_const = new TwoTrackMassKinematicConstraint(jpsi_mass);
//            const SecRecoResult lambReco = secondaryReconstruction(cand, &lambName, jpsi_const);
//            const std::vector<std::string> penQName = { "JPsi/MuPos", "JPsi/MuNeg", "Kaon" }; // "Kaon" is the name of "anti-Proton"
//            const SecRecoResult penQReco = secondaryReconstruction(cand, &penQName);
//            const reco::Vertex* PVptr = findPrimaryVertex(cand_jpsi);
//
//            if ( PVptr ) 
//            {
//                lbv2Br.primaryV.x = PVptr->position().X();
//                lbv2Br.primaryV.y = PVptr->position().Y();
//                lbv2Br.primaryV.z = PVptr->position().Z();
//            }
//
//            if( lambReco.isValid() )
//            {
//                lbv2Br.refitMom.mass    = lambReco.mass();
//                lbv2Br.refitMom.pt      = lambReco.momentum()->transverse();
//                lbv2Br.refitMom.eta     = lambReco.momentum()->eta();
//                lbv2Br.refitMom.phi     = lambReco.momentum()->phi();
//                lbv2Br.refitPos.x       = lambReco.vertex()->position().X();
//                lbv2Br.refitPos.y       = lambReco.vertex()->position().Y();
//                lbv2Br.refitPos.z       = lambReco.vertex()->position().Z();
//                lbv2Br.refitPos.vtxprob = TMath::Prob( lambReco.vertex()->chi2(), lambReco.vertex()->ndof() );
//            }
//
//
//            lbv2Br.jpsiMom.mass      = cand_jpsi->composite().mass();
//            lbv2Br.jpsiMom.pt        = cand_jpsi->composite().pt();
//            lbv2Br.jpsiMom.eta       = cand_jpsi->composite().eta();
//            lbv2Br.jpsiMom.phi       = cand_jpsi->composite().phi();
//            lbv2Br.jpsiPos.x         = cand_jpsi->vertex().position().X();
//            lbv2Br.jpsiPos.y         = cand_jpsi->vertex().position().Y();
//            lbv2Br.jpsiPos.z         = cand_jpsi->vertex().position().Z();
//            lbv2Br.jpsiPos.vtxprob   = TMath::Prob( cand_jpsi->vertex().chi2(), cand_jpsi->vertex().ndof() );
//
//            if ( penQReco.isValid() )
//            {
//                lbv2Br.penQMom.mass    = penQReco.mass();
//                lbv2Br.penQMom.pt      = penQReco.momentum()->transverse();
//                lbv2Br.penQMom.eta     = penQReco.momentum()->eta();
//                lbv2Br.penQMom.phi     = penQReco.momentum()->phi();
//                lbv2Br.penQPos.x       = penQReco.vertex()->position().X();
//                lbv2Br.penQPos.y       = penQReco.vertex()->position().Y();
//                lbv2Br.penQPos.z       = penQReco.vertex()->position().Z();
//                lbv2Br.penQPos.vtxprob = TMath::Prob( penQReco.vertex()->chi2(), penQReco.vertex()->ndof() );
//            }
//
//            const reco::Candidate* pTk = cand->getDaug("Kaon"  );
//            const reco::Candidate* nTk = cand->getDaug("Proton");
//            lbv2Br.pTk.px     = pTk->px(); 
//            lbv2Br.pTk.py     = pTk->py(); 
//            lbv2Br.pTk.pz     = pTk->pz(); 
//            lbv2Br.pTk.charge = pTk->charge();
//            lbv2Br.nTk.px     = nTk->px();
//            lbv2Br.nTk.py     = nTk->py();
//            lbv2Br.nTk.pz     = nTk->pz();
//            lbv2Br.nTk.charge = nTk->charge();
//            lbv2Br.eventNumber= eventNo;
//
//            int tagInfo = 0;
//            if ( lambReco.isValid() )   tagInfo += 1 << 0;
//            if ( penQReco.isValid() )   tagInfo += 1 << 1;
//            if ( PVptr )                tagInfo += 1 << 2;
//
//            if ( tagInfo == (1<<3)-1 )
//                lbv2Tree->Fill();
//        }
//    } // LbToTkTk fill end }}}
}
void LbSpecificDecay::endJob()
{
    return;
}




void LbSpecificDecay::setRecoParameters( const edm::ParameterSet& ps )
{



    const string& name = ps.getParameter<string>( "name" );
    bool writeCandidate = ps.getParameter<bool>( "writeCandidate" );
    switch( rMap[name] )
    {
        case Onia     : recoOnia      = true; writeOnia       = writeCandidate; break;
        case Psi1     :
        case Psi2     : recoOnia      = true;                                   break;
        case Lam0     : recoLam0      = true; writeLam0       = writeCandidate; break;
        case TkTk     : recoTkTk      = true; writeTkTk       = writeCandidate; break;
        case LbToLam0 : recoLbToLam0  = true; writeLbToLam0   = writeCandidate; break;
        case LbToTkTk : recoLbToTkTk  = true; writeLbToTkTk   = writeCandidate; break;
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

DEFINE_FWK_MODULE( LbSpecificDecay );
