import FWCore.ParameterSet.Config as cms

recoSelect = cms.VPSet(
    cms.PSet( name = cms.string( 'JPsi' ),
        # ptMin is muon pt cut
          ptMin = cms.double( 4.0 ),
         etaMax = cms.double( 10.0 ),
        massMin = cms.double( 1.50  ),
        massMax = cms.double( 5.50  ),
        probMin = cms.double( 0.15  ),
        constrMass  = cms.double( 3.096900 ),
        constrSigma = cms.double( 0.000006 ),
        writeCandidate = cms.bool( True )
    ),
    cms.PSet( name = cms.string( 'pTks'  ),
          ptMin = cms.double(  0.8 ),
         etaMax = cms.double(  2.5 ),
        massMin = cms.double(  1.00  ),
        massMax = cms.double(  2.50  ),
        probMin = cms.double( 0.01  ),
        writeCandidate = cms.bool( True )
    ),
    cms.PSet( name = cms.string( 'pL0B' ),
        mPsiMin = cms.double( 3.096 - 0.150 ),
        mPsiMax = cms.double( 3.096 + 0.150 ),
        # ptMin is proton pt cut
        massMin = cms.double( 5.00 ),
        massMax = cms.double( 6.00 ),
        probMin = cms.double( 0.02 ),
     massFitMin = cms.double( 5.40 ),
     massFitMax = cms.double( 5.60 ),
        constrMJPsi = cms.bool(  True ),
        writeCandidate = cms.bool( True )
    ),
    cms.PSet( name = cms.string( 'nTks'  ),
          ptMin = cms.double(  0.8 ),
         etaMax = cms.double(  2.5 ),
        massMin = cms.double(  1.00  ),
        massMax = cms.double(  2.50  ),
        probMin = cms.double( 0.01  ),
        writeCandidate = cms.bool( True )
    ),
    cms.PSet( name = cms.string( 'nL0B' ),
        mPsiMin = cms.double( 3.096 - 0.150 ),
        mPsiMax = cms.double( 3.096 + 0.150 ),
        # ptMin is proton pt cut
        massMin = cms.double( 5.00 ),
        massMax = cms.double( 6.00 ),
        probMin = cms.double( 0.02 ),
     massFitMin = cms.double( 5.40 ),
     massFitMax = cms.double( 5.90 ),
        constrMJPsi = cms.bool(  True ),
        writeCandidate = cms.bool( True )
    ),
    cms.PSet( name = cms.string( 'Lam0'  ),
          ptMin = cms.double(  0.8 ),
         etaMax = cms.double(  2.5 ),
        massMin = cms.double(  0.95  ),
        massMax = cms.double(  1.30  ),
        probMin = cms.double( 0.02  ),
        writeCandidate = cms.bool( True )
    ),
    cms.PSet( name = cms.string( 'LbL0' ),
        mPsiMin = cms.double( 3.096 - 0.150 ),
        mPsiMax = cms.double( 3.096 + 0.150 ),
        mLam0Min= cms.double( 1.115683 - 0.1 ),
        mLam0Max= cms.double( 1.115683 + 0.1 ),
        # ptMin is proton pt cut
        massMin = cms.double( 5.00 ),
        massMax = cms.double( 6.00 ),
        probMin = cms.double( 0.02 ),
     massFitMin = cms.double( 5.40 ),
     massFitMax = cms.double( 5.90 ),
        constrMJPsi = cms.bool(  True ),
        writeCandidate = cms.bool( True )
    ),
    cms.PSet( name = cms.string( 'Lamo'  ),
          ptMin = cms.double(  0.8 ),
         etaMax = cms.double(  2.5 ),
        massMin = cms.double(  0.95  ),
        massMax = cms.double(  1.30  ),
        probMin = cms.double( 0.02  ),
        writeCandidate = cms.bool( True )
    ),
    cms.PSet( name = cms.string( 'LbLo' ),
        mPsiMin = cms.double( 3.096 - 0.150 ),
        mPsiMax = cms.double( 3.096 + 0.150 ),
        mLam0Min= cms.double( 1.115683 - 0.1 ),
        mLam0Max= cms.double( 1.115683 + 0.1 ),
        # ptMin is proton pt cut
        massMin = cms.double( 5.00 ),
        massMax = cms.double( 6.00 ),
        probMin = cms.double( 0.02 ),
     massFitMin = cms.double( 5.40 ),
     massFitMax = cms.double( 5.90 ),
        constrMJPsi = cms.bool(  True ),
        writeCandidate = cms.bool( True )
    ),
)
