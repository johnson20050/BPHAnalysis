#include <iostream>
#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"

#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/format.h"
//#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/filePath.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/filePathmc.h"

#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/cutFuncs.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/histoMAP.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/cutLists.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/fourMom.h"
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/vtxprobCollection.h"
#include <map>
using namespace std;
// funcs{{{
// decide the content filled in histograms.
template<typename myDataBranch>
void sthToFill(histoMAP& inHistos, const myDataBranch* root);
// load tree and cut lists, and execute sthToFill
template<typename myDataBranch>
void plotWithCuts(histoMAP& inHistos, const myDataBranch* root, TTree* tree, const vector<CutList*>& cuts, unsigned numberToFill=0);
// load file by file, and execute plotWithCuts
template<typename myDataBranch>
void plotFromFile(histoMAP& inHistos, myDataBranch* root, vector<CutList*>& cuts, const string& branchName, const string& fileName, unsigned numberToFill=0);
// output a png figure.
void outputFigure( map<string, TH1D*>& hMap, const string& branchName, const string& plotCommand="", const string& anyCommand="");
// output a root file contains lots of histograms.
void outputFigure_rootFile( map<string, TH1D*>& hMap, TFile* fStore,  const string& dirName );

const double   jPsiMass =  3.096916;
const double protonMass = 0.938272046;
const double   pionMass =  0.13957018;
const double   kaonMass =  0.493667;
// funcs end }}}
int main() 
{
    string branchName[2] = { "LbToTkTree", "LbTotkTree" };
    string dirName[2] = { "particle", "antiparticle" };
    string anyCommand = "withCut";
    string fileName = "reassignMass_"+anyCommand+".root";

    LambToTkTkBranches* root = new LambToTkTkBranches();
    // cuts applied
    vector<CutList*> cutLists;
    cutLists.push_back( new      vtxprobCut(0.15,-99. , root->refitPos.getAdd()) );
    cutLists.push_back( new         massCut(5.0 ,  6.0, root->refitMom.getAdd()) );
    cutLists.push_back( new       cosa2dCut(0.99,       root->refitPos.getAdd(), root->refitMom.getAdd(), root->primaryV.getAdd()) );
    cutLists.push_back( new           ptCut(15  ,-99. , root->refitMom.getAdd()) );
    cutLists.push_back( new flightDist2DCut(0.2 ,  0.6, root->refitPos.getAdd(), root->primaryV.getAdd()) );


    TFile* fStore = new TFile( fileName.c_str(), "RECREATE" );

    for ( int i=0; i<2; ++i )
    {
        map<string, TH1D*> totMap;
        createHisto( totMap, "daughter_JPsi",  80  , 0.00, 10.0);
        createHisto( totMap, "daughter_KK"  , 110*2, 0.90, 2.00);
        createHisto( totMap, "daughter_Kp"  , 110*2, 1.40, 2.50);
        createHisto( totMap, "daughter_KPi" , 140*2, 0.60, 2.00);
        createHisto( totMap, "daughter_pPi" , 120*2, 1.00, 2.20);
        createHisto( totMap, "daughter_PiPi", 180*2, 0.20, 2.00);
        createHisto( totMap, "BToJPsiKK"    , 250  , 4.00, 6.50);
        createHisto( totMap, "BToJPsiKp"    , 250  , 4.50, 7.00);
        createHisto( totMap, "BToJPsipK"    , 250  , 4.50, 7.00);
        createHisto( totMap, "BToJPsiKPi"   , 330  , 3.70, 7.00);
        createHisto( totMap, "BToJPsipPi"   , 270  , 4.30, 7.00);
        createHisto( totMap, "BToJPsiPiPi"  , 350  , 3.50, 7.00);

        histoMAP tmpHistos;
        tmpHistos.copyHisto( totMap );
        unsigned numToFill = 2;
        // run over all files, totalFileName is declared in filePath.h
        for ( const auto& fileName : totalFileName )
        {
            plotFromFile( tmpHistos, root, cutLists, branchName[i], fileName , numToFill );
            tmpHistos.storeHistoInto( totMap );
        }

        outputFigure_rootFile( totMap, fStore, dirName[i] );
        deleteHistoMap( totMap );
    }

    printf( "created root file: %s \n", fileName.c_str() );
    for ( auto& myCut : cutLists )
        delete myCut;
    delete root;

    fStore->Close();
} 

template<typename myDataBranch>
void sthToFill(histoMAP& inHistos, const myDataBranch* root)
{
    fourMom jpsi( root->jpsiMom,   jPsiMass );
    fourMom   ptkK( root->    pTk,   kaonMass );
    fourMom   ptkp( root->    pTk, protonMass );
    fourMom   ptkP( root->    pTk,   pionMass );
    fourMom   ntkK( root->    nTk,   kaonMass );
    fourMom   ntkP( root->    nTk,   pionMass );
    fourMom   ntkp( root->    nTk, protonMass );
    fourMom  KK  ,  Kp  ,  pK  ,  KPi ,  pPi ,  PiPi;
    fourMom jKK  , jKp  , jpK  , jKPi , jpPi , jPiPi;
    KK   = ptkK+ntkK;
    Kp   = ptkK+ntkp;
    pK   = ptkp+ntkK;
    KPi  = ptkK+ntkP;
    pPi  = ptkp+ntkP;
    PiPi = ptkP+ntkP;
    jKK  = jpsi+KK;
    jKp  = jpsi+Kp;
    jpK  = jpsi+pK;
    jKPi = jpsi+KPi;
    jpPi = jpsi+pPi;
    jPiPi= jpsi+PiPi;


    inHistos.fillHisto( "daughter_JPsi", jpsi .Mag() );
    inHistos.fillHisto( "daughter_KK"  , KK   .Mag() );
    inHistos.fillHisto( "daughter_Kp"  , Kp   .Mag() );
    inHistos.fillHisto( "daughter_pK"  , pK   .Mag() );
    inHistos.fillHisto( "daughter_KPi" , KPi  .Mag() );
    inHistos.fillHisto( "daughter_pPi" , pPi  .Mag() );
    inHistos.fillHisto( "daughter_PiPi", PiPi .Mag() );
    inHistos.fillHisto( "BToJPsiKK"    , jKK  .Mag() );
    inHistos.fillHisto( "BToJPsiKp"    , jKp  .Mag() );
    inHistos.fillHisto( "BToJPsipK"    , jpK  .Mag() );
    inHistos.fillHisto( "BToJPsiKPi"   , jKPi .Mag() );
    inHistos.fillHisto( "BToJPsipPi"   , jpPi .Mag() );
    inHistos.fillHisto( "BToJPsiPiPi"  , jPiPi.Mag() );
}
template<typename myDataBranch>
void plotWithCuts(histoMAP& inHistos, const myDataBranch* root, TTree* tree, const vector<CutList*>& cuts, unsigned numberToFill)
{
    inHistos.clearHisto();
    tmpVtxprobCollection vtxprobCollection( numberToFill );
    unsigned i = 0;
    unsigned n = tree->GetEntries();
    while( i != n )
    {
        tree->GetEntry(i++);
        bool fillTag = true;

        if ( emptyMom( root->refitMom ) ) continue;
        if ( emptyPos( root->refitPos ) ) continue;
        if ( emptyPos( root->primaryV ) ) continue;

        // cuts are listed in main()
        for ( const auto& cut : cuts )
            if ( !cut->accept() ) fillTag = false;

        if ( fillTag || i==n ) // i==n is for the last event. prevent to lost the last information
        {
            std::pair<double,unsigned> tmpPair( root->refitPos.vtxprob, i-1 );
            vtxprobCollection.fillEntry( tmpPair, root->eventNumber );
            if ( vtxprobCollection.startToFill() || i==n ) 
            {
                const map<double, unsigned,greater<double>>* entry = vtxprobCollection.getEntry();
                map<double, unsigned>::const_iterator iter = entry->begin();
                map<double, unsigned>::const_iterator iend = entry->end  ();
                while ( iter != iend )
                {
                    unsigned _entry = iter->second;
                    unsigned _numberToFill = 0;
                    tree->GetEntry( _entry );

                    sthToFill( inHistos, root );
                    
                    ++iter;
                    if ( ++_numberToFill == vtxprobCollection.getFillNumber() ) break;
                }
                vtxprobCollection.clear();
            }
        }
    }
}

template<typename myDataBranch>
void plotFromFile(histoMAP& inHistos, myDataBranch* root, vector<CutList*>& cuts, const string& branchName, const string& fileName, unsigned numberToFill)
{
    TFile* f1 = TFile::Open( ("file:"+fileName).c_str() );
    TTree* tree = (TTree*)f1->Get( ("lbSpecificDecay/"+branchName).c_str() );
    root->SetBranchAddress( tree );
    plotWithCuts(inHistos, root, tree, cuts, numberToFill);
    
    inHistos.preventDeletedByTFile();
    f1->Close();
    return;
}

void outputFigure( map<string, TH1D*>& hMap, const string& branchName, const string& plotCommand, const string& anyCommand)
{
    TCanvas c("", "", 1600, 1000);
    map<string, TH1D*>::iterator iter = hMap.begin();
    map<string, TH1D*>::iterator iend = hMap.end  ();
    //c.SetFillColor(39);
    //c.GetSelectedPad()->SetFillColor(39);
    
    // draw all histograms in the map
    while ( iter != iend )
    {
        const string& hName = iter->first;
        TH1D*   histo = iter->second;
        histo->SetTitle( (hName+plotCommand).c_str() );
        histo->Draw();
        if ( !plotCommand.empty() )
            c.SetLogy();
        c.SaveAs( ("parFigure_"+branchName+"."+hName+plotCommand+"_"+anyCommand+".png").c_str() );
        ++iter;
    }
}

void outputFigure_rootFile( map<string, TH1D*>& hMap, TFile* fStore,  const string& dirName )
{
    TDirectory* myDir = fStore->mkdir( dirName.c_str() );
    myDir->cd();
    map<string, TH1D*>::iterator iter = hMap.begin();
    map<string, TH1D*>::iterator iend = hMap.end  ();
    
    // draw all histograms in the map
    while ( iter != iend )
    {
        const string& hName = iter->first;
        TH1D*   histo = iter->second;
        histo->SetTitle( hName.c_str() );
        histo->Write();
        ++iter;
    }
    fStore->cd();
    // TDirectory will be deleted by TFile.
}
