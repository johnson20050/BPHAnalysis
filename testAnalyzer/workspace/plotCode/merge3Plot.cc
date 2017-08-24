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
    string treeLabel    = "jpsiTree";
    string branchLabel  = "momentum.mass";
    vector<string> fileLabel = {
        ""                                  , 
        "HLT_DoubleMu4_3_Jpsi_Displaced_v3" , 
        "HLT_Dimuon16_Jpsi_v3"              , 
        "HLT_DoubleMu4_JpsiTrk_Displaced_v3" 
    };
    TH1D* h[4];
    //for(int i=0;i<4;++i)
    //    h[i] = getMass("result"+fileLabel[i]+".root" , treeLabel, branchLabel, 5. , 6.  );
    //  h[i] = getMass("result"+fileLabel[i]+".root" , treeLabel, branchLabel, 2.5, 3.5 );
    //  h[i] = getMass("result"+fileLabel[i]+".root" , treeLabel, branchLabel, 1.1, 1.15);
    //TH1D* histo = new TH1D("", "",  50, 5.,   6. );
    //TH1D* histo = new TH1D("", "",  50, 2.5, 3.5  );

    h[0] = getMass("result"+fileLabel[3]+".root", "lam0Tree", "momentum.mass", 1.1, 1.15);
    h[1] = getMass("result"+fileLabel[3]+".root", "lam0Tree", "refitMom.mass", 1.1, 1.15);
    h[2] = getMass("result"+fileLabel[3]+".root", "lam0Tree_paolo", "refitMom.mass", 1.1, 1.15);
    setHistoPlot(h[0], kBlack  , 2,     999,  999);
    setHistoPlot(h[1], kYellow , 2, kYellow,  999);
    setHistoPlot(h[2],    0    , 2, kRed   , 3244);
    setHistoPlot(h[3], kBlue+1 , 2, kBlue+1, 3351);

    setHistoYRange(h[1], 0, 180);

    TLegend* legend = new TLegend(0.15, 0.65, 0.5, 0.85, "HLT applied: ", "NDC");
//    for(int i=0;i<4;++i)
//        legend->AddEntry(h[i], fileLabel[i].c_str(), "lepf");
    legend->AddEntry(h[0], "lam0Tree.momentum.mass", "lepf");
    legend->AddEntry(h[1], "lam0Tree.refitMom.mass", "lepf");
    legend->AddEntry(h[2], "lam0Tree_paolo.refitMom.mass", "lepf");

    legend->SetBorderSize(0);



    TCanvas c1("", "", 800, 600);

    h[1]->Draw();
    h[2]->Draw("same");
    h[0]->Draw("same");
//    h[3]->Draw("same");
    legend->Draw("same");
//    c1.SaveAs( (treeLabel + "." + branchLabel + ".png").c_str() );
    c1.SaveAs( "comparePlot_lam0Tree.differentReconstruction.png" );

    //delete[] h;
    delete legend;
}




