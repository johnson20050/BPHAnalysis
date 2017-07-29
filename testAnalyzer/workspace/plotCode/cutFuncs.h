#ifndef _CUT_FUNCS_H_
#define _CUT_FUNCS_H_
#include "format.h"
#include <map>
#include <math.h>
#include <iostream>
// apply the cut
// m : set minimum value cut. ( data should be bigger  than value. )
// M : set maximum value cut. ( data should be smaller than value. )

bool emptyMom( const momentumStore& mom)
{ 
    if ( mom.pt == 0 && mom.eta == 0 && mom.phi == 0 ) return true;  
    return false;
}
bool emptyPos( const positionStore& pos)
{ 
    if ( pos.x  == 0 && pos.y   == 0 && pos.z   == 0 ) return true;  
    return false;
}
bool mptCut(const momentumStore& mom, const double cutVal) 
{ return (mom.pt > cutVal); }
bool MptCut(const momentumStore& mom, const double cutVal) 
{ return (mom.pt < cutVal); }
bool MetaCut(const momentumStore& mom, const double cutVal) 
{ return (mom.eta < cutVal); }
bool mmassCut(const momentumStore& mom, const double cutVal)
{ return (mom.mass > cutVal); }
bool MmassCut(const momentumStore& mom, const double cutVal)
{ return (mom.mass < cutVal); }
bool mvtxprobCut(const positionStore& pos, const double cutVal)
{ return (pos.vtxprob > cutVal); }
bool MvtxprobCut(const positionStore& pos, const double cutVal)
{ return (pos.vtxprob < cutVal); }
bool mflightDist2DCut(const positionStore& pos, const positionStore& pv, const double cutVal) 
{ 
    double x = pos.x - pv.x; double y = pos.y - pv.y;
    return ( (x*x+y*y) > cutVal*cutVal ); 
}
bool MflightDist2DCut(const positionStore& pos, const positionStore& pv, const double cutVal) 
{ 
    double x = pos.x - pv.x; double y = pos.y - pv.y;
    return ( (x*x+y*y) < cutVal*cutVal ); 
}
bool mproperTimeCut(const positionStore& pos, const momentumStore& mom, const positionStore& pv, const double cutVal) 
{
    double x = pos.x - pv.x;
    double y = pos.y - pv.y;
    double pos2D = x*x+y*y;
    double pt2 = mom.pt*mom.pt;
    double propertime = pos2D/pt2;
    return (propertime > cutVal*cutVal);
}
bool MproperTimeCut(const positionStore& pos, const momentumStore& mom, const positionStore& pv, const double cutVal) 
{
    double x = pos.x - pv.x;
    double y = pos.y - pv.y;
    double pos2D = x*x+y*y;
    double pt2 = mom.pt*mom.pt;
    double propertime = pos2D/pt2;
    return (propertime < cutVal*cutVal);
}
bool mcosa2d(const positionStore& pos, const momentumStore& mom, const positionStore& pv, const double cutVal)
{
    double x = pos.x - pv.x; double y = pos.y - pv.y; double r = sqrt(x*x+y*y);
    double mx = mom.pt*cos(mom.phi); double my = mom.pt*sin(mom.phi);
    return (mx*x+my*y)/(mom.pt*r) > cutVal;
}
double properTime(const positionStore& pos, const momentumStore& mom, const positionStore& pv) 
{
    double x = pos.x - pv.x;
    double y = pos.y - pv.y;
    double pos2D = sqrt( x*x+y*y );
    double propertime = pos2D/mom.pt;
    return propertime;
}
double flightDist2D(const positionStore& pos1, const positionStore& pos2)
{
    if ( ( pos1.x == 0 && pos1.y == 0 ) || ( pos2.x == 0 && pos2.y == 0 ) ) return -999.;
    double x = pos1.x - pos2.x; double y = pos1.y - pos2.y;
    return ( sqrt(x*x+y*y) ); 
}
double cosa2d(const positionStore& pos, const momentumStore& mom, const positionStore& pv)
{
    double x = pos.x - pv.x; double y = pos.y - pv.y; double r = sqrt(x*x+y*y);
    double mx = mom.pt*cos(mom.phi); double my = mom.pt*sin(mom.phi);
    return (mx*x+my*y)/(mom.pt*r);
}


#endif
