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
    //TH1D* histo = new TH1D("", "",  50, 5.0, 6.00 );
    //TH1D* histo = new TH1D("pt", "",  50, 0, 25.00 );
    TH1D* histo = new TH1D("", "", 180/2, 0.2, 2.00 ); 
    //TH1D* histo = new TH1D("", "",  30*4  , 0.2, 0.50 ); 
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
        fourMom  tot;
        //tot = jpsi + ptk + ntk;
        tot = ptk+ntk;
       
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

    vector<string>::const_iterator iterFileName = totalFileName.begin();
    vector<string>::const_iterator iend         = totalFileName.end  ();

    //TH1D* hTot = new TH1D("", "",  50, 5.0, 6.00 );
    TH1D* hTot = new TH1D("", "", 180/2, 0.2, 2.00 ); 
    //TH1D* hTot = new TH1D("", "",  30*4  , 0.2, 0.50 ); 
    while ( iterFileName != iend )
    {
        TH1* hTmp = plotFromFile<LambToTkTkBranches>(branchName, *iterFileName++);
        hTot->Add(hTmp);
        delete hTmp;
    }
    string name = "daughter_KK";
    hTot->SetTitle( name.c_str() );
    TCanvas c("", "", 1600, 1000);
    hTot->Draw();
    c.SaveAs( ("refitplot_"+name+".png").c_str() );

    delete hTot;
}
// main() end }}}
