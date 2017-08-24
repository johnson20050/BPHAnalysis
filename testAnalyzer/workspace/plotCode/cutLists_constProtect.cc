#ifndef __cutLists__
#define __cutLists__

#include "/wk_cms/ltsai/LbFrame/TEST/CMSSW_8_0_21/src/BPHAnalysis/testAnalyzer/workspace/plotCode/cutFuncs2.h"
class momentumStore;
class positionStore;
class CutList
{
public:
    CutList(const momentumStore* momData, double min, double max) : 
        pos(NULL)   ,mom(momData),pv(NULL)  , dmin(min), dmax(max) {}
    CutList(const positionStore* posData, double min, double max) :
        pos(posData),mom(NULL)   ,pv(NULL)  , dmin(min), dmax(max) {}
    CutList(const positionStore* posData, const positionStore* pvData, double min, double max) :
        pos(posData),mom(NULL)   ,pv(pvData), dmin(min), dmax(max) {}
    CutList(const positionStore* posData, const momentumStore* mom, const positionStore* pvData, double min, double max) : 
        pos(posData),mom(NULL)   ,pv(pvData), dmin(min), dmax(max) {}
    virtual bool accept()=0;

//private:
protected:
    double  dmin, dmax;
    const momentumStore* mom;
    const positionStore* pos;
    const positionStore*  pv;
};

class vtxprobCut : public CutList
{
public:
    vtxprobCut( const positionStore* posData, double min, double max= -99. ) : 
        CutList( posData, min, max ) {}
    bool accept()
    { return ( mvtxprobCut(pos, dmin) ); }
};
class massCut : public CutList
{
    massCut( const momentumStore* momData, double min, double max ) : 
        CutList( momData, min, max ) {}
    bool accept()
    { return ( mmassCut(mom, dmin) && MmassCut(mom, dmax) ); }
};
class cosa2dCut : public CutList
{
public:
    cosa2dCut( const positionStore* pos, const momentumStore* mom, const positionStore* pv, double min ) :
        CutList( pos, mom, pv, min, -99. ) {}
    bool accept()
    { return (mcosa2d( pos, mom, pv, dmin) ); }
};
class ptCut : public CutList
{
public:
    ptCut( const momentumStore* momData, double min, double max ) :
        CutList( momData, min, max ) {}
    bool accept()
    {
        if ( dmin < 0 ) return ( MptCut(mom, dmax) );
        if ( dmax < 0 ) return ( mptCut(mom, dmin) );
        return (mptCut(mom, dmin) && MptCut(mom, dmax) ); 
    }
};
class flightDist2DCut : public CutList
{
public:
    flightDist2DCut( const positionStore* posData, const positionStore* pvData, double min, double max=-99.) :
        CutList( posData, pvData, min, max ) {}
    bool accept()
    {
        if ( dmin < 0 ) return ( MflightDist2DCut(pos, pv, dmax) );
        if ( dmax < 0 ) return ( mflightDist2DCut(pos, pv, dmin) );
        return ( mflightDist2DCut(pos, pv, dmin) && MflightDist2DCut(pos, pv, dmax) ); 
    }
};

#endif
