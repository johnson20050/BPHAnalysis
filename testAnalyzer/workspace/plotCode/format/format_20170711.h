#ifndef __FORMAT_H__
#define __FORMAT_H__
#include "TTree.h"
#include <string>
// Modified :
// 21/March/2017
// 27June6/2017 Modified for new  kind of the Lb->Jpsi Tk Tk code

// dataStore{{{
class momentumStore
{
public:
    double mass, pt, eta, phi;
    void RegTree(TTree* root, const std::string& branchName)
    {
        std::string name[4] = { ".mass", ".pt", ".eta", ".phi" };
        root->Branch( (branchName + name[0]).c_str(), &mass   , (branchName + name[0]+"/D").c_str() );
        root->Branch( (branchName + name[1]).c_str(), &pt     , (branchName + name[1]+"/D").c_str() );
        root->Branch( (branchName + name[2]).c_str(), &eta    , (branchName + name[2]+"/D").c_str() );
        root->Branch( (branchName + name[3]).c_str(), &phi    , (branchName + name[3]+"/D").c_str() );
    }
    void SetBranchAddress(TTree* root, const std::string& branchName)
    {
        std::string name[4] = { ".mass", ".pt", ".eta", ".phi" };
        root->SetBranchAddress( (branchName + name[0]).c_str(), &mass   );
        root->SetBranchAddress( (branchName + name[1]).c_str(), &pt     );
        root->SetBranchAddress( (branchName + name[2]).c_str(), &eta    );
        root->SetBranchAddress( (branchName + name[3]).c_str(), &phi    );
    }
    momentumStore* getAdd() { return this; }
};
class positionStore
{
public:
    double x, y, z, vtxprob;
    void RegTree(TTree* root, const std::string& branchName)
    {
        std::string name[4] = { ".x", ".y", ".z", ".vtxprob" };
        root->Branch( (branchName + name[0]).c_str(), &x      , (branchName + name[0]+"/D").c_str() );
        root->Branch( (branchName + name[1]).c_str(), &y      , (branchName + name[1]+"/D").c_str() );
        root->Branch( (branchName + name[2]).c_str(), &z      , (branchName + name[2]+"/D").c_str() );
        root->Branch( (branchName + name[3]).c_str(), &vtxprob, (branchName + name[3]+"/D").c_str() );
    }
    void SetBranchAddress(TTree* root, const std::string& branchName)
    {
        std::string name[4] = { ".x", ".y", ".z", ".vtxprob" };
        root->SetBranchAddress( (branchName + name[0]).c_str(), &x      );
        root->SetBranchAddress( (branchName + name[1]).c_str(), &y      );
        root->SetBranchAddress( (branchName + name[2]).c_str(), &z      );
        root->SetBranchAddress( (branchName + name[3]).c_str(), &vtxprob);
    }
    positionStore* getAdd() { return this; }
};
class daughterStore
{
public:
    double px, py, pz; 
    int    charge;
    void RegTree(TTree* root, const std::string& branchName)
    {
        std::string name[4] = { ".px", ".py", ".pz" , "charge"};
        root->Branch( (branchName + name[0]).c_str(), &px     , (branchName + name[0]+"/D").c_str() );
        root->Branch( (branchName + name[1]).c_str(), &py     , (branchName + name[1]+"/D").c_str() );
        root->Branch( (branchName + name[2]).c_str(), &pz     , (branchName + name[2]+"/D").c_str() );
        root->Branch( (branchName + name[3]).c_str(), &charge , (branchName + name[3]+"/I").c_str() );
    }
    void SetBranchAddress(TTree* root, const std::string& branchName)
    {
        std::string name[4] = { ".px", ".py", ".pz" , "charge"};
        root->SetBranchAddress( (branchName + name[0]).c_str(), &px     );
        root->SetBranchAddress( (branchName + name[1]).c_str(), &py     );
        root->SetBranchAddress( (branchName + name[2]).c_str(), &pz     );
        root->SetBranchAddress( (branchName + name[3]).c_str(), &charge );
    }
    daughterStore* getAdd() { return this; }
};

// }}}



class JpsiBranches // {{{
{
public:
    momentumStore momentum;
    positionStore position;
    UInt_t eventNumber;

    void RegTree(TTree *root)
    {
        momentum.RegTree(root, "momentum");
        position.RegTree(root, "position");
        root->Branch( "eventNumber", &eventNumber, "eventNumber/i" );
    }
    void SetBranchAddress(TTree *root)
    {
        momentum.SetBranchAddress(root, "momentum");
        position.SetBranchAddress(root, "position");
        root->SetBranchAddress("eventNumber", &eventNumber );
    }
};// }}}

class Lam0Branches// {{{
{
public:
    momentumStore momentum;
    positionStore position;
    momentumStore refitMom;
    positionStore refitPos;
    daughterStore proton, pion;
    UInt_t eventNumber;

    void RegTree(TTree *root)
    {
        momentum.RegTree(root, "momentum");
        position.RegTree(root, "position");
        refitMom.RegTree(root, "refitMom");
        refitPos.RegTree(root, "refitPos");
        proton  .RegTree(root, "proton"  );
        pion    .RegTree(root, "pion"    );  
        root->Branch( "eventNumber", &eventNumber, "eventNumber/i" );
    }
    void SetBranchAddress(TTree *root)
    {
        momentum.SetBranchAddress(root, "momentum");
        position.SetBranchAddress(root, "position");
        refitMom.SetBranchAddress(root, "refitMom");
        refitPos.SetBranchAddress(root, "refitPos");
        proton  .SetBranchAddress(root, "proton"  );
        pion    .SetBranchAddress(root, "pion"    );  
        root->SetBranchAddress("eventNumber", &eventNumber );
    }
};// }}}




class LambToLam0Branches // {{{
{
public:
    positionStore primaryV;
    momentumStore refitMom;
    positionStore refitPos;
    momentumStore lam0Mom;
    positionStore lam0Pos;
    momentumStore jpsiMom;
    positionStore jpsiPos;
    UInt_t eventNumber;
    void RegTree(TTree *root)
    {
        primaryV.RegTree(root, "primaryV");
        refitMom.RegTree(root, "refitMom");
        refitPos.RegTree(root, "refitPos");
        lam0Mom .RegTree(root, "lam0Mom" );
        lam0Pos .RegTree(root, "lam0Pos" );
        jpsiMom .RegTree(root, "jpsiMom" );
        jpsiPos .RegTree(root, "jpsiPos" );
        root->Branch( "eventNumber", &eventNumber, "eventNumber/i" );

    }                          
    void SetBranchAddress(TTree *root)
    {
        primaryV.SetBranchAddress(root, "primaryV");
        refitMom.SetBranchAddress(root, "refitMom");
        refitPos.SetBranchAddress(root, "refitPos");
        lam0Mom .SetBranchAddress(root, "lam0Mom" );
        lam0Pos .SetBranchAddress(root, "lam0Pos" );
        jpsiMom .SetBranchAddress(root, "jpsiMom" );
        jpsiPos .SetBranchAddress(root, "jpsiPos" );
        root->SetBranchAddress("eventNumber", &eventNumber );
    }

}; // }}}
class LambToTkTkBranches // {{{
{
public:
    positionStore primaryV;
    momentumStore refitMom;
    positionStore refitPos;
    momentumStore jpsiMom;
    positionStore jpsiPos;
    momentumStore penQMom;
    positionStore penQPos;
    daughterStore pTk;
    daughterStore nTk;  
    UInt_t eventNumber;
    
    void RegTree(TTree *root)
    {
        primaryV.RegTree(root, "primaryV");
        refitMom.RegTree(root, "refitMom");
        refitPos.RegTree(root, "refitPos");
        jpsiMom .RegTree(root, "jpsiMom" );
        jpsiPos .RegTree(root, "jpsiPos" );
        penQMom .RegTree(root, "penQMom" );
        penQPos .RegTree(root, "penQPos" );
        pTk     .RegTree(root, "pTk"     );
        nTk     .RegTree(root, "nTk"     );
        root->Branch( "eventNumber", &eventNumber, "eventNumber/i" );

    }                          
    void SetBranchAddress(TTree *root)
    {
        primaryV.SetBranchAddress(root, "primaryV");
        refitMom.SetBranchAddress(root, "refitMom");
        refitPos.SetBranchAddress(root, "refitPos");
        jpsiMom .SetBranchAddress(root, "jpsiMom" );
        jpsiPos .SetBranchAddress(root, "jpsiPos" );
        penQMom .SetBranchAddress(root, "penQMom" );
        penQPos .SetBranchAddress(root, "penQPos" );
        pTk     .SetBranchAddress(root, "pTk"     );
        nTk     .SetBranchAddress(root, "nTk"     );
        root->SetBranchAddress("eventNumber", &eventNumber );

    }
}; // }}}

#endif
