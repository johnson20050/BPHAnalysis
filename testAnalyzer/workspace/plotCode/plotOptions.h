#ifndef __plotOptions_H__
#define __plotOptions_H__
#include "TH1D.h"
#include <string>
#include <vector>
#include <map>
#include <stdio.h>
#include <string>
//#include "TColor.h"
#include "Rtypes.h"

/*
enum Rtypes::EColor
{
    kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
    kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
    kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
};
*/
class drawOption
{
public:
    TH1D* histo;
    short lineColor;
    short lineWidth;
    short fillColor;
    short fillStyle;
    double ymin;
    double ymax;
    double xmin;
    double xmax;
    
    
    
    drawOption( TH1D* h );
    drawOption( TH1D* h, unsigned lineMode, unsigned fillMode );
    void Update();
    void SetLineColor( short lColor ) { lineColor = lColor; return; }
    void SetLineWidth( short lWidth ) { lineWidth = lWidth; return; }
    void SetFillColor( short fColor ) { fillColor = fColor; return; }
    void SetFillStyle( short fStyle ) { fillStyle = fStyle; return; }
    void SetYaxisRange( double m, double M=-999. ) { ymin = m; ymax = M; return; }
    void changeLineStyleSets ( unsigned mode=0 );
    void changeFillStyleSets ( unsigned mode=0 );
    void BinDevideBy( long num );
    TH1D* GetHisto() { return histo; }
};

drawOption::drawOption( TH1D* h ) :
    histo(h),
    ymin(-999.),
    ymax(-999.),
    xmin(-999.),
    xmax(-999.)
{
    changeLineStyleSets(0);
    changeFillStyleSets(0);
}
drawOption::drawOption( TH1D* h, unsigned lineMode, unsigned fillMode ) :
    histo(h),
    ymin(-999.),
    ymax(-999.),
    xmin(-999.),
    xmax(-999.)
{
    changeLineStyleSets(lineMode);
    changeFillStyleSets(fillMode);
}
void drawOption::Update()
{
    histo->SetLineColor(lineColor);
    histo->SetLineWidth(lineWidth);
    histo->SetFillColor(fillColor);
    histo->SetFillStyle(fillStyle);
    if ( ymin != -999. )
        histo->SetMinimum(ymin);
    if ( ymax != -999. )
        histo->SetMaximum(ymax);
}
void drawOption::changeLineStyleSets ( unsigned mode )
{
    switch ( mode ) 
    {
        case 0 : // default value;
            lineColor = 602;
            lineWidth = 1;
            break;
        case 1 :
            lineColor = kRed;
            lineWidth = 2;
            break;
        case 2 :
            lineColor = kBlue+2;
            lineWidth = 2;
            break;
        case 3 :
            lineColor = kYellow;
            lineWidth = 2;
            break;
        case 4 :
            lineColor = kSpring;
            lineWidth = 2;
            break;
        case 5 :
            lineColor = kTeal;
            lineWidth = 2;
            break;
        default :
            printf ( "-- drawOption::changeLineStyleSets says you assigned sixth mode, which has no definition.\n");
            lineColor = kBlack;
            lineWidth = 3;
    }
}
void drawOption::changeFillStyleSets ( unsigned mode )
{
    switch ( mode ) 
    {
       case 0 : // default value;
            fillColor = kWhite;
            fillStyle = 1001;
            break;
        case 1 :
            fillColor = kRed;
            fillStyle = 3244;
            break;
        case 2 :
            fillColor = kBlue+1;
            fillStyle = 3351;
            break;
        case 3 :
            fillColor = kYellow;
            fillStyle = 3395;
            break;
        case 4 :
            fillColor = kSpring;
            fillStyle = 3636;
            break;
        case 5 :
            fillColor = kTeal;
            fillStyle = 1001;
            break;
        default :
            printf ( "-- drawOption::changeFillStyleSets says you assigned sixth mode, which has no definition.\n");
            fillColor = kBlack;
            fillStyle = 1001;
    }
}

void drawOption::BinDevideBy( long binFactor )
{
    const std::string name = histo->GetTitle();
    TH1D* tmpHisto = (TH1D*)histo->Rebin( binFactor, (name+"_reBin").c_str() );
    //delete histo;
    histo = tmpHisto;
    return;
}



#endif
