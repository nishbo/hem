#include "topology.h"

int Topology::randomConnections(int N, int *com, double p){
    //if connection exists return 1, 0 otherwise
    for(int i = 0; i<N*N; i++)
      if(VFDistributions::drand() < p)
        com[i] = 1;
      else
        com[i] = 0;
}
int Topology::smallWorld(int N, int *com, int clu, double beta){
    int num, sm;

    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        num = i*N+j;
        if(VFDiscrete::discreteDistanceOnCircle(i,j,N)<clu+1 && i!=j)
          com[num] = 1;
        else
          com[num] = 0;
      }
    }

    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
        num = i*N+j;
        if(com[num])//connection exists
          if(VFDistributions::drand()<beta){    //and we rewrite it
            com[num] = 0;
            while(1){
              sm = rand() % N;
              if(com[i*N +sm]==0){
                com[i*N +sm] = 1;
                break;
              }
            }
          }
      }
    }
}

double Topology::setDelay(double prex, double prey, double posx, double posy, \
                          int type_of_delay, double dt){
    // returns time of signal travel between two neurons depending on
    // type_of_delay

    double d;
    switch (type_of_delay){
    case 0:
        d = dt + VFDistributions::drand() * 10;
        break;
    case 1:
        d = sqrt( (prex - posx) * (prex - posx) + \
                  (prey - posy) * (prey - posy) \
                ) / SPIKE_VELOCITY;
        break;
    case 2:
        d = dt;
        break;
    default:
        d = dt;
    }
    return d;
}
