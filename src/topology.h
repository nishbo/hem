#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include <ctime>
#include <conio.h>
#include <cstdio>

#include "cad.h"
#include "vfdistributions.h"
#include "vfdiscrete.h"

class Topology{
public:
    static int randomConnections(int N, int* com, double p);
    static int smallWorld(int N, int *com, int clu, double beta);
    static double setDelay(double prex, double prey, double posx, double posy,\
                           int type_of_delay, double sv, double dm, double dt);
private:
    Topology();
    Topology(Topology&p){}
    Topology& operator =(Topology&){}
};

#endif // TOPOLOGY_H
