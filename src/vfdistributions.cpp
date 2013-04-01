#include "vfdistributions.h"

int VFDistributions::buf0 = 0;
double VFDistributions::buf1 = 0;
double VFDistributions::buf2 = 0;

double VFDistributions::drand(){
    return ((double) rand() / (RAND_MAX));
}

double VFDistributions::uniform(double min, double max){
    return min + (max - min) * drand();
}

double VFDistributions::normal(double mean, double sd){
    buf0 = 12; //accuracy
    buf2 = 0.;
    for (int i = 0; i < buf0; ++i)
        buf2 += drand();
    buf2 =  buf2 - 6.;

    buf2 = mean + sd * buf2;
    return buf2;
}

double VFDistributions::normal(double mean, double sd, double from, double to){
    if(to < from)
        return 0;
    if(to == from)
        return from;

    while(1){
        buf2 = normal(mean, sd);
        if(buf2 < to && buf2 > from)
            break;
    }
    return buf2;
}

int VFDistributions::poisson(int lambda){
    buf0 = 0;
    buf1 = 0;
    while(1){
        buf1 -= log(drand());
        if(buf1 >= lambda)
            break;
        buf0++;
    }
    return buf0;
}
