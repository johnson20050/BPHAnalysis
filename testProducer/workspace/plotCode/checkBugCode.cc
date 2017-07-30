#include <iostream>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1I.h"
#include "TLegend.h"
using namespace std;


TH1I* getMass(const string& fileName, const string& treeName, const string& branchName, int min = 0., int max = 10.)
{
    TFile* f1 = TFile::Open( ("file:" + fileName).c_str() );
    TTree* tree = (TTree*)f1->Get( ("lbSpecificDecay/" + treeName ).c_str() );
    //TTree* tree = (TTree*)f1->Get("lbSpecificDecay/LbToLam0Tree");

    int value;
    tree->SetBranchAddress( branchName.c_str(), &value );
    TH1I* histo = new TH1I("", "",  20, min, max  );
    unsigned iter = 0;
    unsigned iend = tree->GetEntries();
    while( iter != iend )
    {
        tree->GetEntry(iter++);
        // -1  = total events
        histo->Fill(-1);
        for(int i=0;i<10;++i)
            if((value >> i)%2)
                histo->Fill(i);
    }
    return histo;
}

int main()
{
    string treeLabel    = "LbToLam0Tree";
    string branchLabel  = "refitMom.bugCode";
    vector<string> fileLabel = { "HLT_DoubleMu4_3_Jpsi_Displaced_v3" };
    TH1I* h = getMass("result"+fileLabel[0]+".root" , treeLabel, branchLabel);
    TCanvas c1("", "", 800, 600);

    h->Draw();
    c1.SaveAs( (treeLabel + "." + branchLabel + ".png").c_str() );
}
