#include <iostream>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/format.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/cutFuncs.h"

using namespace std;

template<typename T>
TTree* initialize(T& inputData, const string& fileName, const string& treeName)
{
    TFile* f1 = TFile::Open( ("file:" + fileName).c_str() );
    TTree* tree = (TTree*)f1->Get( ("lbSpecificDecay/" + treeName ).c_str() );
    //TTree* tree = (TTree*)f1->Get("lbSpecificDecay/LbToLam0Tree");

    inputData.SetBranchAddress(tree);
    return tree;
}


TH1D* getMass(const string& fileName, const string& treeName, const string& branchName, double min = 0., double max = 10.)
{
    TFile* f1 = TFile::Open( ("file:" + fileName).c_str() );
    TTree* tree = (TTree*)f1->Get( ("lbSpecificDecay/" + treeName ).c_str() );
    //TTree* tree = (TTree*)f1->Get("lbSpecificDecay/LbToLam0Tree");

    double value;
    tree->SetBranchAddress( branchName.c_str(), &value );
    TH1D* histo = new TH1D("", "",  25, min, max  );
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
    string branchLabel  = "momentum.mass";
    vector<string> fileLabel = {
        ""                                  , 
        "HLT_DoubleMu4_3_Jpsi_Displaced_v3" , 
        "HLT_Dimuon16_Jpsi_v3"              , 
        "HLT_DoubleMu4_JpsiTrk_Displaced_v3" 
    };
    Lam0Branches lam0;
    //TTree* tree = initialize<Lam0Branches>(lam0, "result"+fileLabel[0]+".root", "lam0Tree");

    TFile* f1 = TFile::Open( "file:result.root" );
    //TTree* tree = (TTree*)f1->Get( ("lbSpecificDecay/" + "lam0Tree" ).c_str() );
    TTree* tree = (TTree*)f1->Get( "lbSpecificDecay/lam0Tree" );

    lam0.SetBranchAddress(tree);
    TH1D* h             = new TH1D("", "", 25, 1.1, 1.15);
    TH1D* hproperTime   = new TH1D("", "", 25, 0. , 0.001);
            
    unsigned i = 0;
    unsigned n = tree->GetEntries();
    while ( i != n )
    {
        tree->GetEntry(i++);
        
        //if( !properTimeCut(lam0.refitPos, lam0.refitMom, cutVal) ) continue;
        hproperTime->Fill( properTime(lam0.refitPos, lam0.refitMom) );
        cout << lam0.refitPos.x << endl;
    }


    //setHistoPlot(h[0], kBlack  , 2,     999,  999);
    //setHistoPlot(h[1], kYellow , 2, kYellow,  999);
    //setHistoPlot(h[2], kRed    , 2, kRed   , 3244);
    //setHistoPlot(h[3], kBlue+1 , 2, kBlue+1, 3351);

    //setHistoYRange(h[0], 0);

    //TLegend* legend = new TLegend(0.15, 0.65, 0.5, 0.85, "HLT applied: ", "NDC");
    //for(int i=0;i<4;++i)
    //    legend->AddEntry(h[i], fileLabel[i].c_str(), "lepf");
    //legend->SetBorderSize(0);



    TCanvas c1("", "", 800, 600);

    hproperTime->Draw();
    //h[1]->Draw("same");
    //h[2]->Draw("same");
    //h[3]->Draw("same");
    //legend->Draw("same");
    c1.SaveAs( ("properTime_"+treeLabel + "." + branchLabel + ".png").c_str() );

}




