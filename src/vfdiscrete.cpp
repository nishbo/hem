#include "vfdiscrete.h"

int VFDiscrete::buf0 = 0;
double VFDiscrete::buf1 = 0;
double VFDiscrete::buf2 = 0;

int VFDiscrete::inBetween(double t, double spc, double dt){
    //returns true if t~=spc * n, where n is some int.
    //equal with precision of dt/2
    if(t - ((double) floor((t + dt/5)  / spc)) * spc < dt / 2)
        return 1;
    else
        return 0;
}

int VFDiscrete::discreteDistanceOnCircle(int x, int y, int N){
    //x, y = two elements. N - amount of elements
    if(abs(x-y)<abs(MIN(x,y) + N - MAX(x,y)))
        return (abs(x-y));
    else
        return (abs(MIN(x,y) + N - MAX(x,y)));
}

double VFDiscrete::diracDelta(double t, double dt){
    buf1 = abs(t);
    if(buf1 < dt/2){
        return 1;
    } else
        return 0;
}

double VFDiscrete::heavisideTheta(double t){
    if(t<0)
        return 0.;
    else if(t>0)
        return 1.;
    else
        return 0.5;
}
