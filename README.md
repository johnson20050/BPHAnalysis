# BPHAnalysis
This analysis code is modified from Paolo.<br>
This is the original site:<br>
https://github.com/ronchese/BPHAnalysis<br>
And the Twiki page is:<br>
https://twiki.cern.ch/twiki/bin/viewauth/CMS/BPHReco_V2<br>
<br>
<br>
The analysis code is for these channels:<br>
  Lambda0_b->J/Psi + Lambda0<br>
<br>
  Lambda0_b->PentaQuark + K-<br>
             PentaQuark->J/Psi + p<br>
                         J/Psi-> Mu + Mu<br>
<br>
Usage:
1. RecoDecay is the basic class.
2. SpecificDecay is the derived class.
    it defined EDAnalyzer and EDProducer
3. testAnalyzer is the code I used. It is the EDanalyzer.
    * in plugin/LbSpecificDecay:
        1. Use BPHRecoCandidate and BPHPlusMinusCandidate to reconstruct Lambda0_b.
        2. Use LbSpecificDecay::secondaryReconstruction() to recalculate momentum of each track in the particle reconstructed by step1.
        3. record the result in Tree and output the TFile.
    * in python/LbrecoSelectForWrite_cfi:<br />
        1. record the preselection used in LbSpecificDecay. <br />
    * in workspace/python/LbrecoWrite.py:<br />
        1. The configuration file to analyze 'BPHSkim'.<br />
    * in workspace/python/LbrecoWritevHLT.py:<br />
        1. The configuration file to analyze 'BPHSkim'.<br />
        2. HLT is set.<br />
    * in workspace/python/LbrecoWrite_AOD.py:<br />
        1. The configuration file to analyze 'AOD'.<br />
    * in workspace/python/LbrecoWrite_AODvHLT.py:<br />
        1. The configuration file to analyze 'AOD'.<br />
        2. HLT is set.<br />
    * in workspace/python/LbrecoWrite_miniAOD.py:<br />
        1. The configuration file to analyze 'miniAOD'.<br />

4. The user added variables:
    * refit:
        1. GlobalVector             fitMomentum
        2. reco::Vertex             fitVertex
        3. float                    fitMass
    * get information by Ref<std::vector<T> >: ( load data from other label )
        1. reco::Vertex             primaryVertex
        2. pat::CompositeCandidate  refToJPsi
        3. pat::CompositeCandidate  refTo????
    * other information:
        1. bool                     cowboy
        2. GlobalVector             JPsi/MuPos.fitMom
        3. GlobalVector             JPsi/MuNeg.fitMom
        2. GlobalVector             ????/??tk1.fitMom
        3. GlobalVector             ????/??tk2.fitMom
        2. GlobalVector             JPsi/MuPos.IPt
        3. GlobalVector             JPsi/MuNeg.IPt
        2. GlobalVector             ????/??tk1.IPt
        3. GlobalVector             ????/??tk2.IPt
        2. GlobalVector             JPsi/MuPos.IPt_Error
        3. GlobalVector             JPsi/MuNeg.IPt_Error
        2. GlobalVector             ????/??tk1.IPt_Error
        3. GlobalVector             ????/??tk2.IPt_Error
        2. GlobalVector             JPsi/MuPos.IPz
        3. GlobalVector             JPsi/MuNeg.IPz
        2. GlobalVector             ????/??tk1.IPz
        3. GlobalVector             ????/??tk2.IPz
        2. GlobalVector             JPsi/MuPos.IPz_Error
        3. GlobalVector             JPsi/MuNeg.IPz_Error
        2. GlobalVector             ????/??tk1.IPz_Error
        3. GlobalVector             ????/??tk2.IPz_Error
<br>
<br>
<br>
## Modified by 29/07/2017<br>
## Author: Lian-Sheng, Tsai. in NTUHEP<br>
