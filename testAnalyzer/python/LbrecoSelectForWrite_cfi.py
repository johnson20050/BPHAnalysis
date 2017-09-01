import FWCore.ParameterSet.Config as cms

recoSelect = cms.VPSet(
    cms.PSet( name = cms.string( 'JPsi' ),
        # ptMin is muon pt cut
          ptMin = cms.double( 4.0 ),
         etaMax = cms.double( 10.0 ),
        massMin = cms.double( 1.50   ),
        massMax = cms.double( 5.50   ),
        probMin = cms.double( 0.002 ),
        constrMass  = cms.double( 3.096900 ),
        constrSigma = cms.double( 0.000006 )
    ),
    cms.PSet( name = cms.string( 'Lam0'  ),
        # ptMin is proton pt cut
          ptMin = cms.double( 1.0 ),
         etaMax = cms.double( 10.0 ),
        massMin = cms.double(  0.80  ),
        massMax = cms.double(  1.80  ),
        probMin = cms.double( -1.0  ),
        constrMass  = cms.double(  1.115683 ),
        constrSigma = cms.double( -1.0 )
    ),
    cms.PSet( name = cms.string( 'TkTk'  ),
        # ptMin is proton pt cut
          ptMin = cms.double( 0.8 ),
         etaMax = cms.double(  2.5 ),
        massMin = cms.double(  1.40  ),
        massMax = cms.double(  3.00  ),
        probMin = cms.double( 0.005 ),
    ),
    cms.PSet( name = cms.string( 'LbToLam0' ),
        mPsiMin = cms.double( 3.0019 ),
        mPsiMax = cms.double( 3.2019 ),
        mLam0Min= cms.double( 1.100 ),
        mLam0Max= cms.double( 1.140 ),
        massMin = cms.double( 5.00 ),
        massMax = cms.double( 6.00 ),
        probMin = cms.double( 0.02 ),
        massFitMin = cms.double( 5.00 ),
        massFitMax = cms.double( 6.00 ),
        constrMJPsi = cms.bool(  True )
    ),
    cms.PSet( name = cms.string( 'LbToTkTk' ),
        mPsiMin = cms.double( 3.0019 ),
        mPsiMax = cms.double( 3.2019 ),
        # ptMin is proton pt cut
        massMin = cms.double( 4.50 ),
        massMax = cms.double( 6.50 ),
        probMin = cms.double( 0.02 ),
     massFitMin = cms.double( 5.00 ),
     massFitMax = cms.double( 6.00 ),
     #compCharge = cms.double( 0 )
    ),
)
