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
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/filePath.h"
//#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/filePath_20170630.h"
#include <map>
using namespace std;

string defaultFilePath = "/home/ltsai/Data/CRABdata/CRABdata_6_Jun_2017/result/";

// input value need to change with tree Entry.
// getMass {{{
template<typename myDataValues>
TH1D* getMass(TTree* tree, const myDataValues& val, int bin, double xmin = 0., double xmax = 10.)
{
    TH1D* histo = new TH1D("", "", bin, xmin, xmax  );
    unsigned i = 0;
    unsigned n = tree->GetEntries();
    while( i != n )
    {
        tree->GetEntry(i++);
        histo->Fill(val);
        //histo->Fill(root.refitMom.mass);
    }
    return histo;
} // getMass end }}}
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
    TH1D* histo = new TH1D("", "",  50, 5.0, 6.00 );
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

        //if ( !mmassCut(root.lam0Mom, 1.100 ) ) fillTag = false;
        //if ( !MmassCut(root.lam0Mom, 1.130 ) ) fillTag = false;
        if ( !mmassCut(root.refitMom, 5.0  ) ) fillTag = false;
        if ( !MmassCut(root.refitMom, 6.0  ) ) fillTag = false;
        //if ( !mptCut(root.refitMom, 15) ) fillTag = false;
        if ( !mvtxprobCut(root.refitPos, 0.1) ) fillTag = false;
        if ( !mcosa2d(root.refitPos, root.refitMom, root.primaryV, 0.99) ) fillTag = false;
        if ( !mflightDist2DCut(root.refitPos, root.primaryV, 0.1) ) fillTag = false;
        //if ( !MflightDist2DCut(root.refitPos, root.primaryV, 0.012 ) ) fillTag = false;
        //if ( !mproperTimeCut(root.refitPos, root.refitMom, root.primaryV, 0.04) ) fillTag = false;
        //if ( !MproperTimeCut(root.refitPos, root.refitMom, root.primaryV, 0.13) ) fillTag = false;
        //if ( !mproperTimeCut(root.refitPos, root.refitMom, root.primaryV, 0.03) ) fillTag = false;
        //if ( !mflightDist2DCut(root.lam0Pos, root.refitPos, 20) ) fillTag = false;

        if ( fillTag )
        {
            histo->Fill(root.refitMom.mass);
            //histo->Fill(root.lam0Mom.mass);
        }
     
    }

    return histo;
} // plotMassWithCut end }}}
// plotParameterFigure {{{
template<typename myDataBranch>
map<string, TH1*> plotParameterFigure(TTree* tree, const myDataBranch& root, const string& tag="")
{

    map<string, TH1*> h;
    //h["vtxprob"]        = new TH1D("", "", 100, 0., 0.0000003);
    h["vtxprob"]        = new TH1D("", "", 100, 0.0003, 0.3);
    h["properTime"]     = new TH1D("", "", 150, 0.05, 0.2);
    h["flightDist2D"]   = new TH1D("", "", 100, 0.0, 2.0);
    //h["flightDist2D"]   = new TH1D("", "",  40, 0.11, 0.13);
    h["cosa2d"]         = new TH1D("", "",  50,  0.996, 1.001);
    //h["lam0TolambDist"] = new TH1D("", "",  80, 20 ,  60);

    unsigned i = 0;
    unsigned n = tree->GetEntries();
    while( i != n )
    {
        tree->GetEntry(i++);
        h["vtxprob"]        ->Fill( root.refitPos.vtxprob );
        h["properTime"]     ->Fill( properTime(root.refitPos, root.refitMom, root.primaryV) );  
        h["flightDist2D"]   ->Fill( flightDist2D(root.refitPos, root.primaryV) );
        h["cosa2d"]         ->Fill( cosa2d(root.refitPos, root.refitMom, root.primaryV) );
        //h["lam0TolambDist"] ->Fill( flightDist2D(root.refitPos, root.lam0Pos) );
    }
    return h;

    //TCanvas c1("", "", 1600, 1000);
    //map<string, TH1D*>::iterator iter = h.begin();
    //map<string, TH1D*>::iterator iend = h.end  ();
    //while( iter!=iend )
    //{
    //    //setHistoPlot(iter->second, ....);
    //    iter->second->Draw();
    //    c1.SaveAs( ("parFigure_"+tag+"."+iter++->first+".png").c_str() );
    //}
}// plotParameterFigure end }}}
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

template<typename myDataBranch>
map<string, TH1*> plotParFromFile(const string& branchName, const string& fileName)
{
    TFile* f1 = TFile::Open( ("file:"+fileName).c_str() );
    TTree* tree = (TTree*)f1->Get( ("lbSpecificDecay/"+branchName).c_str() );
    myDataBranch root;
    root.SetBranchAddress(tree);
    map<string, TH1*> hmap;
    hmap["vtxprob"]        = new TH1D("", "", 100, 0.0003, 0.3);
    hmap["properTime"]     = new TH1D("", "",  35, 0.05, 0.4);
    hmap["flightDist2D"]   = new TH1D("", "", 100, 0.1, 0.2);
    hmap["flightDist2D"]   = new TH1D("", "",  40, 0.11, 0.13);
    hmap["cosa2d"]         = new TH1D("", "",  50,  0.996, 1.001);

    hmap = plotParameterFigure( tree, root );
    //size_t n = branchName.size();
    //branchName.erase( n-4, 4); // remove "Tree"
    //hmap = plotParameterFigure( tree, root, branchName );

    map<string, TH1*>::const_iterator iter = hmap.begin();
    map<string, TH1*>::const_iterator iend = hmap.end  ();
    while ( iter != iend )
    {
        iter->second->SetDirectory(0);  // prevent histogram to be deleted by TFile::Close().
        iter++;
    }
    f1->Close();
    return hmap;
} // plotFile end }}}
// main() {{{
int main() // mass with cuts plot
{
    string branchName  = "LbToTkTree";
    //string branchName2 = "LbTotkTree";

    vector<string>::const_iterator iterFileName = totalFileName.begin();
    vector<string>::const_iterator iend         = totalFileName.end  ();

    TH1D* hTot = new TH1D("", "",  50, 5.0, 6.00 );
    //TH1D* hTot = new TH1D("pt", "",  50, 0, 25.00 );
    while ( iterFileName != iend )
    {
        TH1* hTmp = plotFromFile<LambToTkTkBranches>(branchName, *iterFileName++);
        hTot->Add(hTmp);
        delete hTmp;
    }
    TLegend legend(0.15, 0.65, 0.5, 0.85, "HLT applied: ", "NDC");
    legend.AddEntry(hTot, "test legend", "lepf");

    legend.SetBorderSize(0);
    legend.SetFillColor(39);
    TCanvas c("", "", 1600, 1000);
    c.SetFillColor(39);
    hTot->Draw();
    c.SaveAs( ("cut."+branchName+".refitMom.mass_withCuts.png").c_str() );

    delete hTot;
}
//int main() // parFigure plot
//{
//    string branchName = "LbToTkTree";
//
//    vector<string>::const_iterator iterFileName = totalFileName.begin();
//    vector<string>::const_iterator iendFileName = totalFileName.end  ();
//
//    map<string, TH1*> hmap;
//    hmap["vtxprob"]        = new TH1D("", "", 100, 0.0003, 0.3);
//    hmap["properTime"]     = new TH1D("", "", 150, 0.05, 0.2);
//    hmap["flightDist2D"]   = new TH1D("", "", 100, 0.0, 2.0);
//    hmap["cosa2d"]         = new TH1D("", "",  50,  0.996, 1.001);
//    map<string, TH1*>::iterator iter;
//    map<string, TH1*>::iterator iend = hmap.end();
//    while ( iterFileName != iendFileName )
//    {
//        // get histograms from function.
//        map<string, TH1*> tmpMap = plotParFromFile<LambToTkTkBranches>(branchName, *iterFileName++);
//        map<string, TH1*>::const_iterator hIter = tmpMap.begin();
//        map<string, TH1*>::const_iterator hIend = tmpMap.end  ();
//        // sum up the histograms to final histograms
//        while ( hIter != hIend )
//            if ( (iter = hmap.find( hIter->first )) != iend )
//            {
//                iter++->second->Add(hIter->second);
//                delete hIter++->second;
//            }
//            else delete hIter++->second;
//    }
//
//
//    TCanvas c("", "", 1600, 1000);
//    iter = hmap.begin();
//    // draw all histogram in the map
//    while ( iter != iend )
//    {
//        iter->second->SetTitle( iter->first.c_str() );
//        iter->second->Draw();
//        c.SetLogy();
//        c.SaveAs( ("parFigure_"+branchName+"."+iter++->first+".png").c_str() );
//    }
//} // main() end }}}
