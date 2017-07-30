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
//#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/filePath_20170630.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/filePath.h"
#include <map>
using namespace std;
void clearHistoMap( map<string, TH1*>& inHistos );
template<typename myDataBranch>
void plotParameterFigure(map<string, TH1*>& inHistos, TTree* tree, const myDataBranch& root, const string& tag="");
template<typename myDataBranch>
void plotParFromFile(map<string, TH1*>& inHistos, const string& branchName, const string& fileName);
int main() 
{
    string branchName  = "LbToTkTree";
    //string branchName2 = "LbTotkTree";


    map<string, TH1*> totMap;
    //totMap["vtxprob"]        = new TH1D("", "", 100, 0.0003, 0.3);
    totMap["vtxprob"]        = new TH1D("", "", 100, 0.0003, 1.0003);
    //totMap["properTime"]     = new TH1D("", "", 150, 0.05, 0.2);
    totMap["properTime"]     = new TH1D("", "", 95, 0.005, 0.1);
    //totMap["flightDist2D"]   = new TH1D("", "", 100, 0.0, 2.0);
    totMap["flightDist2D"]   = new TH1D("", "", 100, 0.01, 2.00);
    totMap["cosa2d"]         = new TH1D("", "", 500, 1.001-500*0.00001, 1.001);
    totMap["eta"]            = new TH1D("", "", 100, -2.5, 2.5);
    map<string, TH1*> tmpMap;
    tmpMap["vtxprob"]        = (TH1D*)totMap["vtxprob"]     ->Clone(); 
    tmpMap["properTime"]     = (TH1D*)totMap["properTime"]  ->Clone(); 
    tmpMap["flightDist2D"]   = (TH1D*)totMap["flightDist2D"]->Clone(); 
    tmpMap["cosa2d"]         = (TH1D*)totMap["cosa2d"]      ->Clone(); 
    tmpMap["eta"]            = (TH1D*)totMap["eta"]         ->Clone(); 
    

    map<string, TH1*>::iterator iter;
    map<string, TH1*>::iterator iend = totMap.end();
    vector<string>::const_iterator iterFileName = totalFileName.begin();
    vector<string>::const_iterator iendFileName = totalFileName.end  ();

    // fill the histograms
    while ( iterFileName != iendFileName )
    {
        // get histograms from function.
        plotParFromFile<LambToTkTkBranches>(tmpMap, branchName, *iterFileName++);
        map<string, TH1*>::const_iterator _hIter = tmpMap.begin();
        map<string, TH1*>::const_iterator _hIend = tmpMap.end  ();
        // sum up the histograms to final histograms
        while ( _hIter != _hIend )
            if ( (iter = totMap.find( _hIter->first )) != iend )
                iter->second->Add( _hIter++->second );
            else 
                _hIter++;
        //plotParFromFile<LambToTkTkBranches>(tmpMap, branchName2, *iterFileName++);
        //_hIter = tmpMap.begin();
        //// sum up the histograms to final histograms
        //while ( _hIter != _hIend )
        //    if ( (iter = totMap.find( _hIter->first )) != iend )
        //        iter->second->Add( _hIter++->second );
        //    else 
        //        _hIter++;
    }


    string plotCommand = "";
    plotCommand += ".LogPlot";
    TCanvas c("", "", 1600, 1000);
    iter = totMap.begin();
    c.SetFillColor(39);
    // draw all histograms in the map
    while ( iter != iend )
    {
        iter->second->SetTitle( (iter->first+plotCommand).c_str() );
        iter->second->Draw();
        if ( !plotCommand.empty() )
            c.SetLogy();
        c.SaveAs( ("parFigure_"+branchName+"."+iter->first+plotCommand+".png").c_str() );
        ++iter;
    }
} // main() 
// plotParameterFigure {{{
void clearHistoMap( map<string, TH1*>& inHistos )
{
    map<string, TH1*>::iterator iter = inHistos.begin();
    map<string, TH1*>::iterator iend = inHistos.end  ();
    while ( iter != iend )
        iter++->second->Reset();
}
template<typename myDataBranch>
void plotParameterFigure(map<string, TH1*>& inHistos, TTree* tree, const myDataBranch& root, const string& tag)
{
    clearHistoMap( inHistos );
    unsigned i = 0;
    unsigned n = tree->GetEntries();
    while( i != n )
    {
        tree->GetEntry(i++);
        if ( emptyMom( root.refitMom ) ) continue;
        if ( emptyPos( root.refitPos ) ) continue;
        if ( emptyPos( root.primaryV ) ) continue;
        inHistos["vtxprob"]        ->Fill( root.refitPos.vtxprob );
        inHistos["properTime"]     ->Fill( properTime(root.refitPos, root.refitMom, root.primaryV) );  
        inHistos["flightDist2D"]   ->Fill( flightDist2D(root.refitPos, root.primaryV) );
        inHistos["cosa2d"]         ->Fill( cosa2d(root.refitPos, root.refitMom, root.primaryV) );
        inHistos["eta"]            ->Fill( root.refitMom.eta );
        //h["lam0TolambDist"] ->Fill( flightDist2D(root.refitPos, root.lam0Pos) );
    }
}


// testing, is it easy to be read?
//template<typename myDataBranch>
//void plotParameterFigure(map<string, TH1*>& inHistos, TTree* tree, const myDataBranch& root, const string& tag="")
//{
//    clearHistoMap( inHistos );
//    unsigned i = 0;
//    unsigned n = tree->GetEntries();
//    map<string, TH1*>::iterator iend = inHistos.end();
//    map<string, TH1*>::iterator iter;
//    while( i != n )
//    {
//        tree->GetEntry(i++);
//        if ( iter = inHistos.find("vtxprob"      ) != iend )
//            iter->second->Fill( root.refitPos.vtxprob );
//        if ( iter = inHistos.find("properTime"   ) != iend )
//            iter->second->Fill( properTime(root.refitPos, root.refitMom) );  
//        if ( iter = inHistos.find("flightDist2D" ) != iend )
//            iter->second->Fill( flightDist2D(root.refitPos, root.primaryV) );
//        if ( iter = inHistos.find("cosa2d"       ) != iend )
//            iter->second->Fill( cosa2d(root.refitPos, root.refitMom, root.primaryV) );
//        if ( iter = inHistos.find("eta"          ) != iend )
//            iter->second->Fill( root.refitMom.eta );
//        //h["lam0TolambDist"] ->Fill( flightDist2D(root.refitPos, root.lam0Pos) );
//    }
//}// plotParameterFigure end }}}

// plotFile {{{
template<typename myDataBranch>
void plotParFromFile(map<string, TH1*>& inHistos, const string& branchName, const string& fileName)
{
    TFile* f1 = TFile::Open( ("file:"+fileName).c_str() );
    TTree* tree = (TTree*)f1->Get( ("lbSpecificDecay/"+branchName).c_str() );
    myDataBranch root;
    root.SetBranchAddress(tree);
    
    plotParameterFigure(inHistos, tree, root );

    map<string, TH1*>::const_iterator iter = inHistos.begin();
    map<string, TH1*>::const_iterator iend = inHistos.end  ();
    while ( iter != iend )
    {
        iter->second->SetDirectory(0);  // prevent histogram to be deleted by TFile::Close().
        iter++;
    }
    f1->Close();
} // plotFile end }}}

