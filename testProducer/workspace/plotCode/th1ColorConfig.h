#ifndef __PLOT_CONFIG_H__
#define __PLOT_CONFIG_H__
#include "TColor.h"

class plotConfig
{
public:
    plotConfig()
    { 
        lineColor = lineWidth = fillColor = fillStyle = 999; 
        minimum = maximum = -1111.; 
        stats = true; 
    }
    plotConfit(short )


    void setPlotStyle(short lColor,short lWidth,short fColor,short fStyle)
    {
        if ( lColor != 999 ) lineColor = lColor;
        if ( lWidth != 999 ) lineWidth = lWidth;
        if ( fColor != 999 ) fillColor = fColor;
        if ( fStyle != 999 ) fillStyle = fStyle;
    }
    void SetYRange(Double_t min, Double_t max = -1111.)
    { minimum = min; maximum = max; }
    void SetLineColor(short num){ lineColor = num; } 
    void SetLineWidth(short num){ lineWidth = num; } 
    void SetFillColor(short num){ fillColor = num; } 
    void SetFillStyle(short num){ fillStyle = num; } 
    void update(TH1* histo)
    {
        if ( lineColor != 999 ) histo->SetLineColor(lineColor);
        if ( lineWidth != 999 ) histo->SetLineWidth(lineWidth);
        if ( fillColor != 999 ) histo->SetFillColor(fillColor);
        if ( fillStyle != 999 ) histo->SetFillStyle(fillStyle);
                                histo->SetStats(stats);
        if ( minimum != -1111.) histo->SetMinimum(minimum);
        if ( maximum != -1111.) histo->SetMaximum(maximum);
    }


private:
    short lineColor, lineWidth, fillColor, fillStyle;
    Double_t minimum, maximum; // init = -1111
    bool stats;

};


#endif
