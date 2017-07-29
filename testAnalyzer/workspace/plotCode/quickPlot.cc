#include <iostream>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
using namespace std;


TH1D* getMass(const string& fileName, const string& treeName, const string& branchName, double min = 0., double max = 10.)
{
    TFile* f1 = TFile::Open( ("file:" + fileName).c_str() );
    TTree* tree = (TTree*)f1->Get( ("lbSpecificDecay/" + treeName ).c_str() );
    //TTree* tree = (TTree*)f1->Get("lbSpecificDecay/LbToLam0Tree");

    double value;
    tree->SetBranchAddress( branchName.c_str(), &value );
    TH1D* histo = new TH1D("", "",  50, min, max  );
    unsigned iter = 0;
    unsigned iend = tree->GetEntries();
    while( iter != iend )
    {
        tree->GetEntry(iter++);
        histo->Fill(value);
    }
    return histo;
}
void setHistoPlot(TH1* histo, short lineColor, short lineWidth, short fillColor, short fillStyle)
{
    if(lineColor != 999) histo->SetLineColor(lineColor);
    if(lineWidth != 999) histo->SetLineWidth(lineWidth);
    if(fillColor != 999) histo->SetFillColor(fillColor);
    if(fillStyle != 999) histo->SetFillStyle(fillStyle);
}
void setHistoYRange(TH1* histo, Double_t min, Double_t max = -1111.)
{
    histo->SetMinimum(min);
    histo->SetMaximum(max);
}
int main()
{
    string treeLabel    = "lam0Tree";
    string branchLabel  = "refitMom.mass";
    vector<string> fileLabel = {
        "HLT_DoubleMu4_3_Jpsi_Displaced_v3" , 
    };
  //TH1D* h =  getMass("result"+fileLabel[0]+".root" , treeLabel, branchLabel, 5. , 6.  );
    TH1D* h =  getMass("result"+fileLabel[0]+".root" , treeLabel, branchLabel, 1.1, 1.15);


    TCanvas c1("", "", 800, 600);

    h->Draw();
    c1.SaveAs( ("qplot_"+treeLabel + "." + branchLabel + ".png").c_str() );

}





