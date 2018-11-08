#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/Vector3D.h"

        //message logger
#include "FWCore/MessageLogger/interface/MessageLogger.h"
        //magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/Common/interface/Ref.h"
        //ROOT
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
        //others
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include <iostream>
        //Bfinder
#include "BPHAnalysis/Bfinder/interface/format.h"

#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TBranch.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TObject.h>
#include <TFile.h>
#include <math.h>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <list>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677
#define KSHORT_MASS 0.497614
#define KSTAR_MASS  0.89594
#define PHI_MASS    1.019455
#define JPSI_MASS   3.096916
#define PSI2S_MASS  3.686109
#define PROTON_MASS 0.938272046
#define LAMBDA_MASS 1.115683

//
// class declaration
//

class DumpGenInfo : public edm::EDAnalyzer
{//{{{
    public:
        explicit DumpGenInfo(const edm::ParameterSet&);
        ~DumpGenInfo();
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 
    private:
        virtual void beginJob() ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
 
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        
        virtual bool GetAncestor(const reco::Candidate* p);
 
        // ----------member data ---------------------------
        edm::ParameterSet theConfig;
        edm::InputTag genLabel_;

        edm::Service<TFileService> fs;
        TTree *root;
        GenInfoBranches     GenInfo;
    
};//}}}

void DumpGenInfo::beginJob()
{//{{{
    root = fs->make<TTree>("root","root");
    GenInfo.regTree(root);
}//}}}


DumpGenInfo::DumpGenInfo(const edm::ParameterSet& iConfig):theConfig(iConfig)
{//{{{
    //now do what ever initialization is needed
//  TriggersForMatching_= iConfig.getUntrackedParameter<std::vector<std::string> >("TriggersForMatching");
    genLabel_           = iConfig.getParameter<edm::InputTag>("GenLabel");

}//}}}

DumpGenInfo::~DumpGenInfo()
{//{{{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}//}}}

//
// member functions
//

// ------------ method called for each event  ------------
void DumpGenInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //std::cout << "*************************\nReconstructing event number: " << iEvent.id() << "\n";
    using namespace edm;

    //CLEAN all memory
    memset(&GenInfo     ,0x00,sizeof(GenInfo)   );

    GenInfo.size    = 0;
    

    edm::Handle<reco::GenParticleCollection> gens;
    iEvent.getByLabel(genLabel_, gens);

    std::vector<const reco::Candidate *> sel_cands;
    
    for(std::vector<reco::GenParticle>::const_iterator it_gen=gens->begin();
        it_gen != gens->end(); it_gen++){
        if (it_gen->status() > 2 && it_gen->status() != 8) continue;//only status 1, 2, 8(simulated)
        if(GenInfo.size >= MAX_GEN){
            fprintf(stderr,"ERROR: number of gens exceeds the size of array.\n");
            break;;
        }

                bool isGenSignal = false;
                //save target intermediat state particle
                if (
                    abs(int(it_gen->pdgId()/100) % 100) == 4  ||//c menson
                    abs(int(it_gen->pdgId()/100) % 100) == 5  ||//b menson
                    abs(it_gen->pdgId()) == 511 ||//B_0
                    abs(it_gen->pdgId()) == 521 ||//B_+-
                    abs(it_gen->pdgId()) == 531 ||//B_s
                    abs(it_gen->pdgId()) == 130 ||//KL
                    it_gen->pdgId() == 443      ||//Jpsi
                    it_gen->pdgId() == 100443   ||//Psi(2S)
                    it_gen->pdgId() == 553      ||//Upsilon
                    it_gen->pdgId() == 100553     //Upsilon(2S)
                   ) isGenSignal = true;//b, c, s mesons

                if (abs(it_gen->pdgId()) == 13) isGenSignal = true;//all mu

                if (
                    abs(int(it_gen->pdgId()/100) % 100) == 3  ||//s menson
                    abs(it_gen->pdgId()) == 111 || //pi0
                    abs(it_gen->pdgId()) == 211 //pi+
                    ){
                    reco::GenParticle _deRef = (*it_gen);
                    reco::Candidate* Myself = dynamic_cast<reco::Candidate*>(&_deRef);
                    //std::cout<<Myself->pdgId()<<"-----------"<<std::endl;
                    isGenSignal = GetAncestor(Myself);
                }//all pi from a b meson

                if (!isGenSignal) continue;


                GenInfo.index[GenInfo.size]         = GenInfo.size;
                GenInfo.handle_index[GenInfo.size]  = it_gen-gens->begin();
                GenInfo.pt[GenInfo.size]            = it_gen->pt();
                GenInfo.eta[GenInfo.size]           = it_gen->eta();
                GenInfo.phi[GenInfo.size]           = it_gen->phi();
                GenInfo.mass[GenInfo.size]          = it_gen->mass();
                GenInfo.vx[GenInfo.size]            = it_gen->vx();
                GenInfo.vy[GenInfo.size]            = it_gen->vy();
                GenInfo.vz[GenInfo.size]            = it_gen->vz();
                GenInfo.pdgId[GenInfo.size]         = it_gen->pdgId();
                GenInfo.status[GenInfo.size]        = it_gen->status();
                GenInfo.nMo[GenInfo.size]           = it_gen->numberOfMothers();
                GenInfo.nDa[GenInfo.size]           = it_gen->numberOfDaughters();
                GenInfo.mo1[GenInfo.size]           = -1;//To be matched later.
                GenInfo.mo2[GenInfo.size]           = -1;
                GenInfo.da1[GenInfo.size]           = -1;
                GenInfo.da2[GenInfo.size]           = -1;
                GenInfo.size++;
                sel_cands.push_back(&*it_gen);
            }
            //printf("-----*****DEBUG:End of gens loop.\n");

            int geninfo_idx = 0;
            for(std::vector<const reco::Candidate *>::iterator sCands = sel_cands.begin();
                sCands != sel_cands.end(); sCands++){
                geninfo_idx = int(sCands-sel_cands.begin());
                for(int nGenMo = 0; nGenMo < std::min(2,int((*sCands)->numberOfMothers())); nGenMo++){
                //if((*sCands)->numberOfMothers()==1){
                    for(std::vector<const reco::Candidate *>::iterator mCands = sel_cands.begin();
                    mCands != sel_cands.end(); mCands++){
                        if((*sCands)->mother(nGenMo) == *mCands){
                        //if((*sCands)->mother(0) == *mCands){
                            if(nGenMo == 0) GenInfo.mo1[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenMo == 1) GenInfo.mo2[geninfo_idx] = int(mCands-sel_cands.begin());
                        }
                    }
                }
                for(int nGenDa = 0; nGenDa < std::min(2,int((*sCands)->numberOfDaughters())); nGenDa++){
                    for(std::vector<const reco::Candidate *>::iterator mCands = sel_cands.begin();
                    mCands != sel_cands.end(); mCands++){
                        if((*sCands)->daughter(nGenDa) == *mCands){
                            if(nGenDa == 0) GenInfo.da1[geninfo_idx] = int(mCands-sel_cands.begin());
                            if(nGenDa == 1) GenInfo.da2[geninfo_idx] = int(mCands-sel_cands.begin());
                        }
                    }
                }
            }
    
    root->Fill();
    //std::cout<<"filled!\n";
}

// ------------ method called once each job just after ending the event loop  ------------
void DumpGenInfo::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void DumpGenInfo::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void DumpGenInfo::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void DumpGenInfo::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void DumpGenInfo::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DumpGenInfo::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

bool DumpGenInfo::GetAncestor(const reco::Candidate* p)
{
    if(p->numberOfMothers()==0) return false;
    else{
        const reco::Candidate* MyMom = p->mother(0);
        int mpid = abs(MyMom->pdgId());
        if(abs(int(mpid/100) % 100) == 5) return true;
        else return GetAncestor(MyMom);
    }
}

//}}}
//define this as a plug-in
DEFINE_FWK_MODULE(DumpGenInfo);
