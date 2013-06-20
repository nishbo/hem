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
int Topology::m;

int Topology::setCoordinates(const int N, double *x, double *y){
    fle = vf_file::loadFileToString(FILE_TOPOLOGY, DATAFILES);
    border = vf_file::getParameterIni("Border_length_of_box", fle);
    for(int i=0; i<N; i++){
        x[i] = vf_distributions::uniform(0, border);
        y[i] = vf_distributions::uniform(0, border);
    }
    return 0;
}

int Topology::setTopology(const int N, int *Mfull, int **sout, double **delays, \
                            const double *x, const double *y){
    using namespace vf_file;
    fle = loadFileToString(FILE_TOPOLOGY, DATAFILES);

    type = getParameterIni("Type_of_topology", fle);
    switch(type){
    case 0:
        randomTopology(N, Mfull, sout);
        break;
    case 1:
        smallWorldTopology(N, Mfull, sout);
        break;
    case 2:
        fromOneTopology(N, Mfull, sout);
        break;
    default:
        randomTopology(N, Mfull, sout);
    }

    delay_type = getParameterIni("Type_of_delay", fle);
    *delays = new double[*Mfull];
    switch(delay_type){
    case 0:
        setDelaysRandom(N, *Mfull, *delays);
        break;
    case 1:
        setDelaysCoord(N, sout, x, y, *delays);
        break;
    case 2:
        setDelaysdt(N, *Mfull, *delays);
        break;
    default:
        setDelaysdt(N, *Mfull, *delays);
        break;
    }

    return 0;
}

int Topology::randomTopology(const int N, int *Mfull, int **sout){
    r_p = vf_file::getParameterIni("Probability_of_connection", fle);

    int* arr = new int[N];
    int M, k;
    *Mfull = 0;

    for(int i = 0; i<N; i++){
        M = 0;
        for(int j=0; j < N; j++){
            if(vf_distributions::uniform(0, 1) < r_p){
                arr[j] = 1;
                M++;
            } else
                arr[j] = 0;
        }
        sout[i] = new int[M+1];
        sout[i][0] = M;
        *Mfull += M;
        k = 1;
        for(int j=0; j < N; j++){
            if(arr[j]){
                sout[i][k] = j;
                k++;
            }
        }
    }
    free(arr);
    return 0;
}

int Topology::smallWorldTopology(const int N, int *Mfull, int **sout){
    smw_beta = vf_file::getParameterIni("smw_beta", fle);
    smw_local = vf_file::getParameterIni("smw_local", fle);

    int* arr = new int[N];
    int M, k;
    *Mfull = 0;

    for(int i=0; i<N; i++){
        M = 0;
        for(int j=0; j<N; j++){
            if(vf_discrete::discreteDistanceOnCircle(i, j, N) < smw_local + 1 && i!=j){
                arr[j] = 1;
                M++;
            } else
                arr[j] = 0;
        }
        sout[i] = new int[M+1];
        sout[i][0] = M;
        *Mfull += M;
        k = 1;
        for(int j=0; j < N; j++){
            if(arr[j]){
                sout[i][k] = j;
                k++;
            }
        }
    }

    for(int i=0; i < N; i++){
        for(int j=0; j < sout[i][0]; j++){
            if(vf_distributions::uniform(0, 1) < smw_beta){
                sout[i][j+1] = vf_distributions::uniform(0, N-1);
            }
        }
    }

    free(arr);
    return 0;
}

int Topology::fromOneTopology(const int N, int *Mfull, int **sout){
    m = vf_file::getParameterIni("FROM_ONE_AMOUNT", fle);

    sout[0] = new int [(N-1)*m+1];
    sout[0][0] = (N-1)*m;
    for(int i=1; i < N; i++){
        for(int j=1; j < m+1; j++){
            sout[0][(i-1)*m+j] = i;
        }
    }
    for(int i=1; i < N; i++){
        sout[i] = new int[1];
        sout[i][0] = 0;
    }
    *Mfull = (N-1)*m;
    return 0;
}

int Topology::setDelaysdt(const int N, const int Mfull, double *delays){
    min_delay = vf_file::getParameterIni("Delay_min", fle);

    for(int i=0; i < Mfull; i++)
        delays[i] = min_delay;

    return 0;
}

int Topology::setDelaysRandom(const int N, const int Mfull, double *delays){
    min_delay = vf_file::getParameterIni("Delay_min", fle);
    max_delay = vf_file::getParameterIni("Delay_max", fle);

    for(int i=0; i < Mfull; i++)
        delays[i] = vf_distributions::uniform(min_delay, max_delay);

    return 0;
}

int Topology::setDelaysCoord(const int N, int const* const* sout, const double *x, \
                                const double *y, double *delays){
    velocity = vf_file::getParameterIni("Spike_velocity", fle);
    min_delay = vf_file::getParameterIni("Delay_min", fle);

    int k = 0;
    for(int i=0; i < N; i++){
        for(int j=0; j < sout[i][0]; j++){
            delays[k] = min_delay + sqrt( (x[i] - x[sout[i][j+1]]) * \
                                            (x[i] - x[sout[i][j+1]]) + \
                                          (y[i] - y[sout[i][j+1]]) * \
                                            (y[i] - y[sout[i][j+1]]) \
                                                ) / velocity;
        }
    }

    return 0;
}
