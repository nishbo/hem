#ifndef VFDISCRETE_H
#define VFDISCRETE_H

/// Contains various helpful functions to work with discrete metrics and stuff.

#include <cmath>
#include <cstdlib>
#include "cad.h"

class VFDiscrete{
private:
    static int buf0;
    static double buf1, buf2;

    VFDiscrete();
public:
    static int inBetween(double t, double spc, double dt);
    static int discreteDistanceOnCircle(int x, int y, int N);
    static double diracDelta(double t, double dt);
    static double heavisideTheta(double t);
};

#endif // VFDISCRETE_H
