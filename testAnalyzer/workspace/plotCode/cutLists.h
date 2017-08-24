#ifndef __cutLists__
#define __cutLists__

#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/cutFuncs.h"
class momentumStore;
class positionStore;
class CutList
{
public:
    //CutList(const momentumStore momData, const double min, const double max) : 
    CutList(const double min, const double max, const momentumStore* momData) : 
        pos()       ,mom(momData),pv()      , dmin(min), dmax(max) {}
    CutList(const double min, const double max, const positionStore* posData) :
        pos(posData),mom()       ,pv()      , dmin(min), dmax(max) {}
    CutList(const double min, const double max, const positionStore* posData, const positionStore* pvData) :
        pos(posData),mom()       ,pv(pvData), dmin(min), dmax(max) {}
    CutList(const double min, const double max, const positionStore* posData, const momentumStore* momData, const positionStore* pvData) : 
        pos(posData),mom(momData),pv(pvData), dmin(min), dmax(max) {}
    virtual bool accept()=0;

protected:
    const double  dmin, dmax;
    const momentumStore* mom;
    const positionStore* pos;
    const positionStore*  pv;
};

class vtxprobCut : public CutList
{
public:
    vtxprobCut( const double min, const double max, const positionStore* posData ) : 
        CutList( min, max, posData ) {}
    vtxprobCut( const double min, const positionStore* posData ) : 
        CutList( min, -99., posData ) {}
    bool accept()
    { 
        if ( dmax < 0. ) return ( mvtxprobCut(*pos, dmin) ); 
        return ( mvtxprobCut(*pos, dmin) && MvtxprobCut(*pos, dmax) ); 
    }
};
class massCut : public CutList
{
public:
    massCut( const double min, const double max, const momentumStore* momData ) : 
        CutList( min, max, momData ) {}
    bool accept()
    { return ( mmassCut(*mom, dmin) && MmassCut(*mom, dmax) ); }
};
class cosa2dCut : public CutList
{
public:
    cosa2dCut( const double min, const positionStore* pos, const momentumStore* mom, const positionStore* pv ) :
        CutList( min, -99., pos, mom, pv ) {}
    bool accept()
    { return (mcosa2d( *pos, *mom, *pv, dmin) ); }
};
class ptCut : public CutList
{
public:
    ptCut( const double min, const double max, const momentumStore* momData ) :
        CutList( min, max, momData ) {}
    bool accept()
    {
        if ( dmin < 0 ) return ( MptCut(*mom, dmax) );
        if ( dmax < 0 ) return ( mptCut(*mom, dmin) );
        return (mptCut(*mom, dmin) && MptCut(*mom, dmax) ); 
    }
};
class flightDist2DCut : public CutList
{
public:
    flightDist2DCut( const double min, const double max, const positionStore* posData, const positionStore* pvData ) :
        CutList( min, max, posData, pvData ) {}
    flightDist2DCut( const double min, const positionStore* posData, const positionStore* pvData ) :
        CutList( min, -99., posData, pvData ) {}
    bool accept()
    {
        if ( dmin < 0 ) return ( MflightDist2DCut(*pos, *pv, dmax) );
        if ( dmax < 0 ) return ( mflightDist2DCut(*pos, *pv, dmin) );
        return ( mflightDist2DCut(*pos, *pv, dmin) && MflightDist2DCut(*pos, *pv, dmax) ); 
    }
};

#endif
