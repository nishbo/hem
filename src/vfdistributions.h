#ifndef VFDISTRIBUTIONS_H
#define VFDISTRIBUTIONS_H

/// Contains various helpful functions to work with distributions and
/// random numbers

#include <cstdlib>
#include <cmath>

class VFDistributions{
private:
    static int buf0;
    static double buf1, buf2;

    VFDistributions();
public:
    static double drand();
    static double uniform(double min, double max);
    static double normal(double mean, double sd);
    static double normal(double mean, double sd, double from, double to);
    static int poisson(int lambda);
};

namespace vf_distributions{
    double uniform(double min, double max);
    double normal(double mean, double sd);
    double normal(double mean, double sd, double from, double to);
    int poisson(int lambda);

}

#endif // VFDISTRIBUTIONS_H
