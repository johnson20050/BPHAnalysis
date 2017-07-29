struct Pos
{
    Pos(double a, double b): x(a), y(b){}
    Pos(): x(0), y(0){}
    double radius(){ return sqrt(x*x+y*y);}
    double x, y;
};


double distCalc(Pos p1, Pos p2)
//{ return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) ); }
{ return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) ); }
double distCalc_Normalization(Pos target, Pos origin)
{ return ( distCalc(target, origin) / origin.radius() ); }



// distPlot {{{
void distPlot()
{
    TFile* f1 = TFile::Open("resultTMP.root");
    TTree* tree = (TTree*)f1->Get("Lamb/LambTree");

    Pos beforePos, afterPos;
    tree->SetBranchAddress("PosX"   ,&beforePos.x);
    tree->SetBranchAddress("PosY"   ,&beforePos.y);
    tree->SetBranchAddress("vtxPosX",&afterPos .x);
    tree->SetBranchAddress("vtxPosY",&afterPos .y);
    TH1D* histoPos = new TH1D("Pos", "", 50, 0., 0.006);


    Pos beforeMom, afterMom;
    tree->SetBranchAddress("Px"     ,&beforeMom.x);
    tree->SetBranchAddress("Py"     ,&beforeMom.y);
    tree->SetBranchAddress("vtxPx"  ,&afterMom .x);
    tree->SetBranchAddress("vtxPy"  ,&afterMom .y);
    TH1D* histoMom = new TH1D("Mom", "", 50, 0., 2.);

    for(int i=0;i!=tree->GetEntries();++i)
    {
        tree->GetEntry(i);
        histoPos->Fill( distCalc(beforePos, afterPos) );
        histoMom->Fill( distCalc(beforeMom, afterMom) );
    }

    TCanvas c1;
    c1.Divide(1,2);
    c1.cd(1);
    histoPos->Draw();
    c1.cd(2);
    histoMom->Draw();

    c1.SaveAs("figure_smallplot.png");
    //c1.SaveAs("Figure/figure_smallplot.eps");

    
}
// destPlot end }}}

// distPlot_compareToJPsi {{{
void distPlot_compareToJPsi()
{
    TFile* f1 = TFile::Open("result.root");
    TTree* tree = (TTree*)f1->Get("Lamb/LambTree");

    Pos beforePos, afterPos, jpsi_Pos;
    double vtxmass;
    tree->SetBranchAddress("PosX"    ,&beforePos.x);
    tree->SetBranchAddress("PosY"    ,&beforePos.y);
    tree->SetBranchAddress("vtxPosX" ,&afterPos .x);
    tree->SetBranchAddress("vtxPosY" ,&afterPos .y);
    tree->SetBranchAddress("jpsiPosX",&jpsi_Pos .x);
    tree->SetBranchAddress("jpsiPosY",&jpsi_Pos .y);

    tree->SetBranchAddress("vtxmass" ,&vtxmass    );

    TH1D* histoPos1 = new TH1D("Pos_origTojpsi (%)", "", 100, 0., 100.);
    TH1D* histoPos2 = new TH1D("Pos_fitTojpsi (%)" , "", 100, 0., 100.);
    TH1D* histoMASS = new TH1D("MASS", "", 100, 4., 6. );



    for(int i=0;i!=tree->GetEntries();++i)
    {
        tree->GetEntry(i);
        double distRatio_vtxTojpsi_beforeRefit = 100 * distCalc_Normalization(beforePos, jpsi_Pos);
        double distRatio_vtxTojpsi_afterRefit  = 100 * distCalc_Normalization( afterPos, jpsi_Pos);

        histoPos1->Fill( distRatio_vtxTojpsi_beforeRefit );
        histoPos2->Fill( distRatio_vtxTojpsi_afterRefit  );
    }

    TCanvas c1;
    c1.Divide(1,2);
    c1.cd(1);
    histoPos1->Draw();
    c1.cd(2);
    histoPos2->Draw();

    c1.SaveAs("figure_smallplot.png");
    //c1.SaveAs("Figure/figure_smallplot.eps");

    
}
// destPlot end }}}

void quicksmallPlot()
{
    distPlot_compareToJPsi();
}

