#ifndef __FORMAT_H__
#define __FORMAT_H__
#include "TTree.h"
#include <string>
// Modified :
// 21/March/2017

// dataStore{{{
class momentumStore
{
public:
    double mass, pt, eta, phi;
    int bugCode;
    void RegTree(TTree* root, const std::string& branchName)
    {
        std::string name[5] = { ".mass", ".pt", ".eta", ".phi" , ".bugCode"};
        root->Branch( (branchName + name[0]).c_str(), &mass   , (branchName + name[0]+"/D").c_str() );
        root->Branch( (branchName + name[1]).c_str(), &pt     , (branchName + name[1]+"/D").c_str() );
        root->Branch( (branchName + name[2]).c_str(), &eta    , (branchName + name[2]+"/D").c_str() );
        root->Branch( (branchName + name[3]).c_str(), &phi    , (branchName + name[3]+"/D").c_str() );
        root->Branch( (branchName + name[4]).c_str(), &bugCode, (branchName + name[4]+"/I").c_str() );
    }
    void SetBranchAddress(TTree* root, const std::string& branchName)
    {
        std::string name[5] = { ".mass", ".pt", ".eta", ".phi" , ".bugCode"};
        root->SetBranchAddress( (branchName + name[0]).c_str(), &mass   );
        root->SetBranchAddress( (branchName + name[1]).c_str(), &pt     );
        root->SetBranchAddress( (branchName + name[2]).c_str(), &eta    );
        root->SetBranchAddress( (branchName + name[3]).c_str(), &phi    );
        root->SetBranchAddress( (branchName + name[4]).c_str(), &bugCode);
    }
};
class positionStore
{
public:
    double x, y, z, vtxprob;
    int bugCode;
    void RegTree(TTree* root, const std::string& branchName)
    {
        std::string name[5] = { ".x", ".y", ".z", ".vtxprob", ".bugCode" };
        root->Branch( (branchName + name[0]).c_str(), &x      , (branchName + name[0]+"/D").c_str() );
        root->Branch( (branchName + name[1]).c_str(), &y      , (branchName + name[1]+"/D").c_str() );
        root->Branch( (branchName + name[2]).c_str(), &z      , (branchName + name[2]+"/D").c_str() );
        root->Branch( (branchName + name[3]).c_str(), &vtxprob, (branchName + name[3]+"/D").c_str() );
        root->Branch( (branchName + name[4]).c_str(), &bugCode, (branchName + name[4]+"/I").c_str() );
    }
    void SetBranchAddress(TTree* root, const std::string& branchName)
    {
        std::string name[5] = { ".x", ".y", ".z", ".vtxprob", ".bugCode" };
        root->SetBranchAddress( (branchName + name[0]).c_str(), &x      );
        root->SetBranchAddress( (branchName + name[1]).c_str(), &y      );
        root->SetBranchAddress( (branchName + name[2]).c_str(), &z      );
        root->SetBranchAddress( (branchName + name[3]).c_str(), &vtxprob);
        root->SetBranchAddress( (branchName + name[4]).c_str(), &bugCode);
    }
};
class daughterStore
{
public:
    double pt;
    int    charge;
    void RegTree(TTree* root, const std::string& branchName)
    {
        std::string name[2] = { ".pt", ".charge" };
        root->Branch( (branchName + name[0]).c_str(), &pt     , (branchName + name[0]+"/D").c_str() );
        root->Branch( (branchName + name[1]).c_str(), &charge , (branchName + name[1]+"/I").c_str() );
    }
    void SetBranchAddress(TTree* root, const std::string& branchName)
    {
        std::string name[2] = { ".pt", ".charge" };
        root->SetBranchAddress( (branchName + name[0]).c_str(), &pt     );
        root->SetBranchAddress( (branchName + name[1]).c_str(), &charge );
    }
};

// }}}



class JpsiBranches // {{{
{
public:
    momentumStore momentum;
    positionStore position;

    void RegTree(TTree *root)
    {
        momentum.RegTree(root, "momentum");
        position.RegTree(root, "position");
    }
    void SetBranchAddress(TTree *root)
    {
        momentum.SetBranchAddress(root, "momentum");
        position.SetBranchAddress(root, "position");
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

    void RegTree(TTree *root)
    {
        momentum.RegTree(root, "momentum");
        position.RegTree(root, "position");
        refitMom.RegTree(root, "refitMom");
        refitPos.RegTree(root, "refitPos");
        proton  .RegTree(root, "proton"  );
        pion    .RegTree(root, "pion"    );  
    }
    void SetBranchAddress(TTree *root)
    {
        momentum.SetBranchAddress(root, "momentum");
        position.SetBranchAddress(root, "position");
        refitMom.SetBranchAddress(root, "refitMom");
        refitPos.SetBranchAddress(root, "refitPos");
        proton  .SetBranchAddress(root, "proton"  );
        pion    .SetBranchAddress(root, "pion"    );  
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
    int validInfo;
    void RegTree(TTree *root)
    {
        primaryV.RegTree(root, "primaryV");
        refitMom.RegTree(root, "refitMom");
        refitPos.RegTree(root, "refitPos");
        lam0Mom .RegTree(root, "lam0Mom" );
        lam0Pos .RegTree(root, "lam0Pos" );
        jpsiMom .RegTree(root, "jpsiMom" );
        jpsiPos .RegTree(root, "jpsiPos" );

        root->Branch("validInfo", &validInfo, "validInfo/I");
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

        root->SetBranchAddress("validInfo", &validInfo);
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
    daughterStore proton;
    daughterStore kaon;  
    int validInfo;
    
    void RegTree(TTree *root)
    {
        primaryV.RegTree(root, "primaryV");
        refitMom.RegTree(root, "refitMom");
        refitPos.RegTree(root, "refitPos");
        jpsiMom .RegTree(root, "jpsiMom" );
        jpsiPos .RegTree(root, "jpsiPos" );
        penQMom .RegTree(root, "penQMom" );
        penQPos .RegTree(root, "penQPos" );
        proton  .RegTree(root, "proton"  );
        kaon    .RegTree(root, "kaon"    );

        root->Branch("validInfo", &validInfo, "validInfo/I");
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
        proton  .SetBranchAddress(root, "proton"  );
        kaon    .SetBranchAddress(root, "kaon"    );

        root->SetBranchAddress("validInfo", &validInfo);
    }
}; // }}}

#endif
