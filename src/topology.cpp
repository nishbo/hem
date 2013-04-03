#include "topology.h"

using namespace std;

string Topology::fle;
int Topology::type;
int Topology::delay_type;
double Topology::r_p;
double Topology::smw_beta;
int Topology::smw_local;
double Topology::border;
double Topology::velocity;
double Topology::max_delay;
double Topology::min_delay;

int Topology::setTopology(const int N, double *x, double *y, int *con, double *delays){
    using namespace vf_file;
    fle = loadFileToString(FILE_TOPOLOGY, DATAFILES);

    setCoordinates(N, x, y);

    type = getParameterIni("Type_of_topology", fle);
    switch(type){
    case 0:
        randomTopology(N, con);
        break;
    case 1:
        smallWorldTopology(N, con);
        break;
    default:
        randomTopology(N, con);
    }

    delay_type = getParameterIni("Type_of_delay", fle);
    switch(delay_type){
    case 0:
        setDelaysRandom(N, con, delays);
        break;
    case 1:
        setDelaysCoord(N, con, x, y, delays);
        break;
    case 2:
        setDelaysdt(N, con, delays);
        break;
    default:
        setDelaysdt(N, con, delays);
    }

    return 0;
}

int Topology::setCoordinates(const int N, double *x, double *y){
    border = vf_file::getParameterIni("Border_length_of_box", fle);
    for(int i=0; i<N; i++){
        x[i] = vf_distributions::uniform(0, border);
        y[i] = vf_distributions::uniform(0, border);
    }
    return 0;
}

int Topology::randomTopology(const int N, int *con){
    r_p = vf_file::getParameterIni("Probability_of_connection", fle);
    for(int i = 0; i<N*N; i++)
        if(vf_distributions::uniform(0, 1) < r_p)
            con[i] = 1;
        else
            con[i] = 0;
    return 0;
}

int Topology::smallWorldTopology(const int N, int *con){
    int num, sm;
    smw_beta = vf_file::getParameterIni("smw_beta", fle);
    smw_local = vf_file::getParameterIni("smw_local", fle);

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            num = i*N + j;
            if(vf_discrete::discreteDistanceOnCircle(i, j, N) < smw_local + 1 && i!=j)
                con[num] = 1;
            else
                con[num] = 0;
        }
    }

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            num = i*N + j;
            if(con[num]){  //connection exists
                if(vf_distributions::uniform(0, 1)<smw_beta){    //and we rewrite it
                    con[num] = 0;
                    while(1){
                        sm = vf_distributions::uniform(0, N-1);
                        if(con[i*N + sm] == 0){
                            con[i*N + sm] = 1;
                            break;
                        }
                    }
                }
            }
        }
    }
    return 0;
}

int Topology::setDelaysdt(const int N, const int *con, double *delays){
    min_delay = vf_file::getParameterIni("Delay_min", fle);

    for(int i=0; i<N*N; i++)
        if(con[i]>0)
            delays[i] = min_delay;
        else
            delays[i] = 0;
    return 0;
}

int Topology::setDelaysRandom(const int N, const int *con, double *delays){
    min_delay = vf_file::getParameterIni("Delay_min", fle);
    max_delay = vf_file::getParameterIni("Delay_max", fle);

    for(int i=0; i<N*N; i++)
        if(con[i]>0)
            delays[i] = vf_distributions::uniform(min_delay, max_delay);
        else
            delays[i] = 0;
    return 0;
}

int Topology::setDelaysCoord(const int N, const int *con, const double *x, const double *y, double *delays){
    velocity = vf_file::getParameterIni("Spike_velocity", fle);
    min_delay = vf_file::getParameterIni("Delay_min", fle);

    for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            if(con[i*N+j]>0)
                delays[i*N+j] = min_delay + sqrt( (x[i] - x[j]) * (x[i] - x[j]) + \
                                                  (y[i] - y[j]) * (y[i] - y[j]) ) / velocity;
            else
                delays[i*N+j] = 0;

    return 0;
}
