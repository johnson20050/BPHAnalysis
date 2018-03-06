import FWCore.ParameterSet.Config as cms

recoSelect = cms.VPSet(
    cms.PSet( name = cms.string( 'JPsi' ),
        # ptMin is muon pt cut
          ptMin = cms.double( 4.0 ),
         etaMax = cms.double( 10.0 ),
        massMin = cms.double( 1.50  ),
        massMax = cms.double( 5.50  ),
        probMin = cms.double( 0.1   ),
        constrMass  = cms.double( 3.096900 ),
        constrSigma = cms.double( 0.000006 )
    ),
    cms.PSet( name = cms.string( 'TkTk'  ),
          ptMin = cms.double( -1.0 ),
         etaMax = cms.double(  2.5 ),
        massMin = cms.double(  1.00  ),
        massMax = cms.double(  2.25  ),
        probMin = cms.double( 0.10  ),
        #constrMass  = cms.double(  1.115683 ),
        #constrSigma = cms.double( -1.0 )
    ),
    cms.PSet( name = cms.string( 'LbToTkTk' ),
        mPsiMin = cms.double( 2.80   ),
        mPsiMax = cms.double( 3.40   ),
        #mTkTkMin= cms.double( 1.00 ),
        #mTkTkMax= cms.double( 2.25 ),
        # ptMin is proton pt cut
        massMin = cms.double( 5.00 ),
        massMax = cms.double( 6.00 ),
        probMin = cms.double( 0.10 ),
     massFitMin = cms.double( 5.00 ),
     massFitMax = cms.double( 6.00 ),
        constrMJPsi = cms.bool(  True )
    ),
)
