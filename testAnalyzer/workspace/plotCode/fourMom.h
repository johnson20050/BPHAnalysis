#ifndef __FOURMOM__
#define __FOURMOM__
#include "math.h"
#include "stdio.h"


class momentumStore;
class daughterStore;
struct fourMom
{
    fourMom(const momentumStore& in, double mass)
    {
        x1 = in.pt*cos(in.phi); x2 = in.pt*cos(in.phi); x3 = in.pt*sinh(in.eta); 
        x0 = sqrt( x1*x1+x2*x2+x3*x3+mass*mass);
    }
    fourMom(const daughterStore& in, double mass )
    { 
        x1 = in.px; x2 = in.py; x3 = in.pz; 
        x0 = sqrt( x1*x1+x2*x2+x3*x3+mass*mass);
    }
    fourMom(double v0, double v1, double v2, double v3) : x0(v0), x1(v1), x2(v2), x3(v3) {}
    fourMom() { x0 = x1 = x2 = x3 = 0; }

    double Mag()
    { return sqrt(x0*x0-x1*x1-x2*x2-x3*x3); }

    fourMom operator+(const fourMom& in)
    { 
        fourMom tmp;
        tmp.setVar(x0+in.x0, x1+in.x1, x2+in.x2, x3+in.x3);
        return tmp;
    }
    fourMom operator-(const fourMom& in)
    { 
        fourMom tmp;
        tmp.setVar(x0-in.x0, x1-in.x1, x2-in.x2, x3-in.x3);
        return tmp;
    }
    fourMom& operator=(const fourMom& in)
    { x0 = in.x0; x1 = in.x1; x2 = in.x2; x3 = in.x3; return *this; }

    void setVar(const double& v0, const double& v1, const double& v2, const double& v3)
    { x0 = v0; x1 = v1; x2 = v2; x3 = v3; }
    void ShowVal()
    { printf("(E,Px,Py,Pz) = (%4.3f, %4.3f, %4.3f, %4.3f)\n", x0, x1, x2, x3); }
        


    private:
    double x0,x1,x2,x3;
};

#endif
