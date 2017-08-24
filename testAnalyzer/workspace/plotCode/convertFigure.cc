#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include <vector>
#include <map>
#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/plotOptions.h"
using namespace std;

int main()
{
    //string fileName = "cut.refitMom.mass_withCut.root";
    //vector<string> histoName = { "Lambda0_b Mass", "jpsi Mass","pentaQuark Mass" };
    
    //string fileName = "parFigure_withCut.root";
    //vector<string> histoName = { "vtxprob", "properTime", "flightDist2D", "cosa2d", "eta" }; 

    string fileName = "reassignMass_withCut.root";
    vector<string> histoName = { "daughter_KK" , "daughter_Kp"  , "daughter_KPi", 
                                 "daughter_pPi", "daughter_PiPi",
                                 "BToJPsiKK"   , "BToJPsiKp"    , "BToJPsipK"   , 
                                 "BToJPsiKPi"  , "BToJPsipPi"   , "BToJPsiPiPi"  };
    string fileType[] = { ".png", ".pdf", ".eps" };

    
    
    vector<string> dirName = { "particle", "antiparticle" };
    TFile* fStore = TFile::Open( fileName.c_str() );
    
    for ( const auto& dName : dirName )
    {
        TDirectory* dir = fStore->GetDirectory( dName.c_str() );
        map<string, TH1D*> hMap;
        for ( const auto& hName : histoName )
        {
            hMap[hName] = NULL;
            dir->GetObject( hName.c_str(), hMap[hName] );

            drawOption dOpt( hMap[hName], 1, 0 );
            dOpt.SetYaxisRange( 0. );
            dOpt.BinDevideBy( 2 );
            dOpt.Update();
            TCanvas c( "", "", 1600, 1000 );
            //c.Divide(2,1);
            //c.cd(1);
            //hMap[hName]->Draw();
            //c.cd(2);
            dOpt.GetHisto()->GetXaxis()->SetTitle("GeV");
            dOpt.GetHisto()->Draw();
            c.SaveAs( ("convertFigure_"+dName+"_"+hName+fileType[2]).c_str() );
        }

    }
    fStore->Close();
}




