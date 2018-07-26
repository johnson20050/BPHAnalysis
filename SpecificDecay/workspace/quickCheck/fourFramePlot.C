#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "format.h"


void fourFramePlot()
{
    TFile* f1 = TFile::Open("histo.root");
    TTree* tree = (TTree*)f1->Get("bphHistoSpecificDecay/lambTree");

    LambToPPiBranches aa;
    aa.SetBranchAddress( tree );

    TH1F* h[4];
    //for(int i=0;i<4;++i)
    //    h[i] = new TH1F("", "",  30,   5.,  8.);
    h[0] = new TH1F("", "", 30, 5., 6.);
    h[1] = new TH1F("", "", 60,  20.,80.);
    h[2] = new TH1F("", "", 20, - 4., 4.);
    h[3] = new TH1F("", "", 20, - 5., 5.);

    for(int i=0;i!=tree->GetEntries();++i)
    {
        tree->GetEntry(i);
        h[0]->Fill(aa.mass);
        h[1]->Fill(aa.pt  );
        h[2]->Fill(aa.eta );
        h[3]->Fill(aa.phi );
    }

    TCanvas c1;
    c1.Divide(2,2);

    h[0]->SetNameTitle("nameTest", "Lb mass");
    h[1]->SetNameTitle("nameTest", "Lb pt  ");
    h[2]->SetNameTitle("nameTest", "Lb eta ");
    h[3]->SetNameTitle("nameTest", "Lb phi ");
    for(int i=0;i<4;++i)
    {
        c1.cd(i+1);
        h[i]->Draw();
    }

    c1.SaveAs("figurevfourFramePlot.png");
    //c1.SaveAs("Figure/figure_smallplot.eps");

    delete tree;
    delete f1;
    
}
