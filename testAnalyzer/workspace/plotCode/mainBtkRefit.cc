#include <iostream>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/format.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/fourMom.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/cutFuncs.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/filePath.h"
//#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/filePathmc.h"
#include <map>
using namespace std;

string defaultFilePath = "/home/ltsai/Data/CRABdata/CRABdata_6_Jun_2017/result/";
const double   jPsiMass =  3.096916;
const double protonMass = 0.938272046;
const double   pionMass =  0.13957018;
const double   kaonMass =  0.493667;

// input value need to change with tree Entry.
// setHistoPlot {{{
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
} // setHistoPlot end }}}
// plotMassWithCut {{{
template<typename myDataBranch>
TH1* plotMassWithCut(TTree* tree, const myDataBranch& root)
{
    //TH1D* histo = new TH1D("", "",  80, 1.1, 1.15 );
    //TH1D* histo = new TH1D("", "",  150/2, 4.50, 6.00 );
    TH1D* histo = new TH1D("", "",  250/6, 3.50, 6.00 );
    //TH1D* histo = new TH1D("", "",  90, 0.6, 1.50 );
    //TH1D* histo = new TH1D("pt", "",  50, 0, 25.00 );
    unsigned i = 0;
    unsigned n = tree->GetEntries();
    while( i != n )
    {
        tree->GetEntry(i++);
        bool fillTag = true;
        if (  emptyMom(root.refitMom) ) fillTag = false;
        if (  emptyPos(root.refitPos) ) fillTag = false;
        if (  emptyPos(root.primaryV) ) fillTag = false;

        if ( !mmassCut(root.refitMom, 5.0  ) ) fillTag = false;
        if ( !MmassCut(root.refitMom, 6.0  ) ) fillTag = false;
        if ( !mvtxprobCut(root.refitPos, 0.15) ) fillTag = false;
        if ( !mcosa2d(root.refitPos, root.refitMom, root.primaryV, 0.99) ) fillTag = false;
        if ( !mflightDist2DCut(root.refitPos, root.primaryV, 0.2) ) fillTag = false;
        if ( !MflightDist2DCut(root.refitPos, root.primaryV, 0.6) ) fillTag = false;

        fourMom jpsi( root.jpsiMom, jPsiMass );
        fourMom  ptk( root.    pTk,   kaonMass );
        fourMom  ntk( root.    nTk,   kaonMass );
        fourMom  tot, tot2;
        tot = jpsi + ptk + ntk;
        tot2 = ptk+ntk;
        double val = tot2.Mag();
        //if ( val < 0.75 || val > 0.95 ) fillTag = false; //K Pi resonence cut
        //if ( val < 0.3  || val > 0.4  ) fillTag = false; //PiPi resonence cut
        //if ( val < 1.00 || val > 1.20 ) fillTag = false; //K K  resonence cut
        //if ( val < 1.11 || val > 1.12 ) fillTag = false; //p Pi resonence cut
        
        if ( fillTag )
        {
            //histo->Fill(root.refitMom.mass);
            histo->Fill( tot.Mag() );
        }
     
    }

    return histo;
} // plotMassWithCut end }}}
// plotFile {{{
template<typename myDataBranch>
TH1* plotFromFile(const string& branchName, const string& fileName)
{
    TFile* f1 = TFile::Open( ("file:"+fileName).c_str() );
    TTree* tree = (TTree*)f1->Get( ("lbSpecificDecay/"+branchName).c_str() );
    myDataBranch root;

    root.SetBranchAddress(tree);
    TH1* h = plotMassWithCut(tree, root);
    h->SetDirectory(0);  // prevent histogram is deleted by TFile.
    f1->Close();
    return h;
} 
// plotFile}}}
// main() {{{
int main() // mass with cuts plot
{
    string branchName  = "LbToTkTree";
    string branchName2 = "LbTotkTree";

    vector<string>::const_iterator iterFileName = totalFileName.begin();
    vector<string>::const_iterator iend         = totalFileName.end  ();

    TH1D* hTot = new TH1D("", "",  250/6, 3.50, 6.00 );
    //TH1D* hTot = new TH1D("", "",  150/2, 4.50, 6.00 );
    //TH1D* hTot = new TH1D("", "",  90, 0.6, 1.50 );
    while ( iterFileName != iend )
    {
        TH1* hTmp = plotFromFile<LambToTkTkBranches>(branchName, *iterFileName  );
        TH1* hTmp2= plotFromFile<LambToTkTkBranches>(branchName2, *iterFileName++);
        hTot->Add(hTmp);
        hTot->Add(hTmp2);
        delete hTmp;
    }
    string name = "BToJPsiKK_withCut";
    hTot->SetTitle( name.c_str() );
    TCanvas c("", "", 1600, 1000);
    //c.SetFillColor(39);
    hTot->Draw();
    //c.SaveAs( ("cut."+branchName+".refitTkTktoKK_withCuts.png").c_str() );
    c.SaveAs( ("refitplot_"+name+"_twoTree.png").c_str() );
    delete hTot;
}
// main() end }}}
