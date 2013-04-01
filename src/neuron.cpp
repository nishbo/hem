#include "neuron.h"

/// NEURON

void Neuron::addCurrent(double a){
    I += a;
}

std::string Neuron::getName(){
    return "Name not set";
}

double* Neuron::exportData(){
    return NULL;
}

int Neuron::importData(double *arr){
    exit(11);
    return 0;
}

int Neuron::initNeurons(){
    return 1;
}

int Neuron::numEssentialVariables(){
    return 0;
}

int Neuron::setCoordinates(double xo, double yo){
    x = xo;
    y = yo;
    return 0;
}

/// Leaky integrate-and-fire neuron:
std::string NeuronLeakyIAF::neurotype = "leaky_integrate-and-fire_(euler)";
double NeuronLeakyIAF::init_Rin = 1.0;
double NeuronLeakyIAF::init_tau_m = 29.0;
double NeuronLeakyIAF::init_Vth = 15.0;
double NeuronLeakyIAF::init_Vreset = 13.5;
double NeuronLeakyIAF::init_Vrest = 0.0;
double NeuronLeakyIAF::init_V = 0;
double NeuronLeakyIAF::init_tau_ref_abs_exc = 3.0;
double NeuronLeakyIAF::init_tau_ref_abs_inh = 2.0;

std::string NeuronLeakyIAF::getName(){
    return neurotype;
}

int NeuronLeakyIAF::initNeuronsLocal(){
    FILE* fid = fopen("./init/NeuronLeakyIAF.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "Rin = %f\n", &buf01);
    init_Rin = buf01;
    fscanf(fid, "tau_m = %f\n", &buf01);
    init_tau_m = buf01;
    fscanf(fid, "Vth = %f\n", &buf01);
    init_Vth = buf01;
    fscanf(fid, "Vreset = %f\n", &buf01);
    init_Vreset = buf01;
    fscanf(fid, "Vrest = %f\n", &buf01);
    init_Vrest = buf01;
    fscanf(fid, "V = %f\n", &buf01);
    init_V = buf01;
    fscanf(fid, "tau_ref_excitatory = %f\n", &buf01);
    init_tau_ref_abs_exc = buf01;
    fscanf(fid, "tau_ref_abs_inhibitory = %f\n", &buf01);
    init_tau_ref_abs_inh = buf01;

    fclose(fid);
    return 0;
}

int NeuronLeakyIAF::initNeurons(){
    return initNeuronsLocal();
}

void NeuronLeakyIAF::setExcitatory(int f){
    Rin = init_Rin;
    tau_m = init_tau_m;
    Vth = init_Vth;
    Vreset = init_Vreset;
    Vrest = init_Vrest;
    V = init_V;
    if(f){
        exc = 1;
        tau_ref_abs = init_tau_ref_abs_exc;
    } else {
        exc = 0;
        tau_ref_abs = init_tau_ref_abs_inh;
    }
    last_spiked = -(tau_ref_abs + 1);
    I = 0.0;
}

int NeuronLeakyIAF::evolve(double dt, double time){
    V += dt * ( - (V - Vrest) + Rin * I) / tau_m;

    // 'I' sets to 0 to let synapses to add current on the next time-step
    I = 0.0;

    if( V > Vth-1e-8 && time > last_spiked + tau_ref_abs){
        // spike
        V = Vreset;
        last_spiked = time;
        return 1;
    } else if ( time < last_spiked + tau_ref_abs )
        V = Vreset;

    return 0;
}

double* NeuronLeakyIAF::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = V;
    arr[2] = Rin;
    arr[3] = tau_m;
    arr[4] = Vth;
    arr[5] = Vreset;
    arr[6] = tau_ref_abs;
    arr[7] = exc;
    return arr;
}

int NeuronLeakyIAF::importData(double *arr){
    V = arr[0];
    Rin = arr[1];
    tau_m = arr[2];
    Vth = arr[3];
    Vreset = arr[4];
    tau_ref_abs = arr[5];
    exc = arr[6];
    return 0;
}

int NeuronLeakyIAF::numEssentialVariables(){
    return 7;
}

/// Leaky integrate-and-fire neuron Runge-Kutta:
std::string NeuronLeakyIAFRK::neurotype = "leaky_integrate-and-fire_(rk4)";

std::string NeuronLeakyIAFRK::getName(){
    return neurotype;
}

void NeuronLeakyIAFRK::setExcitatory(int f){
    /** Data from "The neural code between neocortical pyramidal neurons depends
    on neurotransmitter release probability" By Tsodyks and Markram, 1997
    and "Synchrony Generation in Recurrent Networks with Frequency-Dependent
    Synapses", 2000
    CHANGE: from http://nest-initiative.org/index.php/Network_burst_\
    generation_by_short-term_plasticity
    */
    Rin = 1.0;    //1 GOhm
    Cm = 30.0;      //30 pF
    tau_m = 30.0;     //30 ms
    // Rin * Cm == tau_m true
    Vth = 15.0;     //15 mV
    Vreset = 13.5;  //13.5 mV
    Vrest = 0.0;    //0 mV
    V = VFDistributions::drand() * Vreset;
    if(f){
        exc = 1;
        Vei = 0;
        tau_ref_abs = 3.;
        tau_ref_tot = 0;
    } else {
        exc = 0;
        Vei = 1;
        tau_ref_abs = 2.;
        tau_ref_tot = 0;
    }
    last_spiked = -(MAX(tau_ref_abs,tau_ref_tot)+1);
}

double NeuronLeakyIAFRK::Vr(double V1, double t, double dt){
    V1 = ( - (V1 - Vrest) + Rin * I) / tau_m;
    if(V1 >= Vth && t > last_spiked + tau_ref_abs)
        return Vreset;
    else
        return V1;
}

int NeuronLeakyIAFRK::evolve(double dt, double time){

    k_1_V = dt * Vr(V, time, dt);
    k_2_V = dt * Vr(V + k_1_V/2.0, time + dt/2.0, dt);
    k_3_V = dt * Vr(V + k_2_V/2.0, time + dt/2.0, dt);
    k_4_V = dt * Vr(V + k_3_V, time + dt, dt);

    V += (k_1_V + 2.0 * k_2_V + 2.0 * k_3_V + k_4_V) * 6.0;

    // 'I' sets to 0 to let synapses to add current on the next time-step
    I = 0.0;

    if( V >= Vth && time > last_spiked + tau_ref_abs &&\
                    time > last_spiked + tau_ref_tot){
        // spike
        V = Vreset;
        last_spiked = time;
        return 1;
    } else if ( time < last_spiked + tau_ref_abs )
        V = Vreset;

    return 0;
}

/// Hodgkin-Huxley neuron:
std::string NeuronHodgkinHuxley::neurotype = "Hodgkin-Huxley_model_(euler)";

std::string NeuronHodgkinHuxley::getName(){
    return neurotype;
}

void NeuronHodgkinHuxley::setExcitatory(int f){
    if(f){
        exc = 1;
    } else {
        exc = 0;
    }
    last_spiked = -100;
    Vth = 80.0;
    tau_spike = 3;

    Vrest = 0.0;
    V = Vrest;

    g_Na = 120;  //�S
    g_K = 36;
    g_L = 0.3;
    E_Na = -115; // mV
    E_K = 12;
    E_L = -10.5989;
    C_mem = 1;  // �F/cm^2

    n = 0.3177;
    m = 0.0529;
    h = 0.5961;
    I = 0.0;
}

void NeuronHodgkinHuxley::addCurrent(double a){
    I+= a / 1000; //transforming pA to microA
}

double NeuronHodgkinHuxley::mF(double V, double m, double dt){
    if (V == -25){
        b = 4.*exp(V/18);
        return m + dt*(1.0 * (1 - m) - b * m);
    } else {
        a = 0.1*(V + 25)/(exp(0.1*(V + 25))-1);
        b = 4*exp(V/18);
        return m + dt*(a * (1 - m) - b * m);
    }
}

double NeuronHodgkinHuxley::hF(double V, double h, double dt){
    a = 0.07 * exp(V/20);
    b = 1/(exp(0.1*(V + 30)) + 1);
    return h + dt*(a * (1 - h) - b * h);
}

double NeuronHodgkinHuxley::nF(double V, double n, double dt){
    if (V == -10){
        b = 0.125*exp(V/80);
        return n + dt*(0.1*(1 - n) - b*n);
    } else {
        a = 0.01*(V+10)/(exp(0.1*(V+10))-1);
        b = 0.125*exp(V/80);
        return n + dt*(a*(1 - n) - b*n);
    }
}

int NeuronHodgkinHuxley::evolve(double dt, double time){
    m = mF(V, m, dt);
    h = hF(V, h, dt);
    n = nF(V, n, dt);

    I_Na = g_Na * m*m*m * h * (V - E_Na);
    I_K  = g_K * n*n*n*n * (V - E_K);
    I_L  = g_L * (V - E_L);
    I += I_Na + I_K + I_L;

    V -= I * dt / C_mem;

    I = 0;
    if(-V>=Vth && last_spiked + tau_spike <time){
        last_spiked = time;
        return 1;
    }
    return 0;
}

/// Hodgkin-Huxley neuron RK4:
std::string NeuronHodgkinHuxleyRK::neurotype = "Hodgkin-Huxley_model_(rk4)";
double NeuronHodgkinHuxleyRK::init_Vth = 80.0;
double NeuronHodgkinHuxleyRK::init_Vrest = 0.0;
double NeuronHodgkinHuxleyRK::init_V = 0;
double NeuronHodgkinHuxleyRK::init_tau_spike = 3;
double NeuronHodgkinHuxleyRK::init_g_Na = 120;
double NeuronHodgkinHuxleyRK::init_g_K = 36;
double NeuronHodgkinHuxleyRK::init_g_L = 0.3;
double NeuronHodgkinHuxleyRK::init_E_Na = -115;
double NeuronHodgkinHuxleyRK::init_E_K = 12;
double NeuronHodgkinHuxleyRK::init_E_L = -10.5989;
double NeuronHodgkinHuxleyRK::init_C_mem = 1;
double NeuronHodgkinHuxleyRK::init_n = 0.3177;
double NeuronHodgkinHuxleyRK::init_m = 0.0529;
double NeuronHodgkinHuxleyRK::init_h = 0.5961;
double NeuronHodgkinHuxleyRK::init_divider = 2.2;

std::string NeuronHodgkinHuxleyRK::getName(){
    return neurotype;
}

int NeuronHodgkinHuxleyRK::initNeuronsLocal(){
    FILE* fid = fopen("./init/NeuronHodgkinHuxley.ini", "r");
    if(!fid){
        exit(29);
    }
    float buf01;
    fscanf(fid, "Vth = %f\n", &buf01);
    init_Vth = buf01;
    fscanf(fid, "Vrest = %f\n", &buf01);
    init_Vrest = buf01;
    fscanf(fid, "V = %f\n", &buf01);
    init_V = buf01;
    fscanf(fid, "tau_spike = %f\n", &buf01);
    init_tau_spike = buf01;
    fscanf(fid, "g_Na = %f\n", &buf01);
    init_g_Na = buf01;
    fscanf(fid, "g_K = %f\n", &buf01);
    init_g_K = buf01;
    fscanf(fid, "g_L = %f\n", &buf01);
    init_g_L = buf01;
    fscanf(fid, "E_Na = %f\n", &buf01);
    init_E_Na = buf01;
    fscanf(fid, "E_K = %f\n", &buf01);
    init_E_K = buf01;
    fscanf(fid, "E_L = %f\n", &buf01);
    init_E_L = buf01;
    fscanf(fid, "C_mem = %f\n", &buf01);
    init_C_mem = buf01;
    fscanf(fid, "n = %f\n", &buf01);
    init_n = buf01;
    fscanf(fid, "m = %f\n", &buf01);
    init_m = buf01;
    fscanf(fid, "h = %f\n", &buf01);
    init_h = buf01;
    fscanf(fid, "divider = %f\n", &buf01);
    init_divider = buf01;

    fclose(fid);
    return 0;
}

int NeuronHodgkinHuxleyRK::initNeurons(){
    return initNeuronsLocal();
}

void NeuronHodgkinHuxleyRK::setExcitatory(int f){
    if(f){
        exc = 1;
    } else {
        exc = 0;
    }
    last_spiked = -100;

    Vth = init_Vth;
    tau_spike = init_tau_spike;
    Vrest = init_Vrest;
    V = init_V;
    g_Na = init_g_Na;
    g_K = init_g_K;
    g_L = init_g_L;
    E_Na = init_E_Na;
    E_K = init_E_K;
    E_L = init_E_L;
    C_mem = init_C_mem;
    n = init_n;
    m = init_m;
    h = init_h;
    divider = init_divider;
    I = 0.0;
}

void NeuronHodgkinHuxleyRK::addCurrent(double a){
    I += a ;
}

double NeuronHodgkinHuxleyRK::mFRK(double V, double m, double h, \
                                   double n, double dt){
    if (V == 25){
        b = 4.*exp(-V/18);
        return 1.0 * (1 - m) - b * m;
    } else {
        a = 0.1*(25-V)/(exp(0.1*(25-V))-1);
        b = 4*exp(-V/18);
        return a * (1 - m) - b * m;
    }
}

double NeuronHodgkinHuxleyRK::hFRK(double V, double m, double h, \
                                   double n, double dt){
    a = 0.07 * exp(-V/20);
    b = 1./(exp(0.1*(30-V))+1);
    return a * (1 - h) - b * h;
}

double NeuronHodgkinHuxleyRK::nFRK(double V, double m, double h, \
                                   double n, double dt){
    if (V == 10){
        b = 0.125*exp(-V/80);
        return 0.1*(1 - n) - b*n;
    } else {
        a = 0.01*(10-V)/(exp(0.1*(10-V))-1);
        b = 0.125*exp(-V/80);
        return a*(1 - n) - b*n;
    }
}

double NeuronHodgkinHuxleyRK::VFRK(double V, double m, double h, \
                                   double n, double dt){
    return I/C_mem - (g_Na * m*m*m*h * (V - E_Na) + \
                      g_K  * n*n*n*n * (V - E_K) + \
                      g_L  * (V - E_L)) / C_mem;
}

int NeuronHodgkinHuxleyRK::subEvolve(double t, double dt){
    if(dt > 0.01){
        Ndt += subEvolve(t, dt/2) + subEvolve(t, dt/2);
        return Ndt;
    } else {
        k_1_h = dt * hFRK(V, m, h, n, dt);
        k_1_n = dt * nFRK(V, m, h, n, dt);
        k_1_m = dt * mFRK(V, m, h, n, dt);
        k_1_V = dt * VFRK(V, m, h, n, dt);
        k_2_h = dt * hFRK(V + k_1_V/2.0, m + k_1_m/2.0, h + k_1_h/2.0, n + k_1_n/2.0, dt);
        k_2_n = dt * nFRK(V + k_1_V/2.0, m + k_1_m/2.0, h + k_1_h/2.0, n + k_1_n/2.0, dt);
        k_2_m = dt * mFRK(V + k_1_V/2.0, m + k_1_m/2.0, h + k_1_h/2.0, n + k_1_n/2.0, dt);
        k_2_V = dt * VFRK(V + k_1_V/2.0, m + k_1_m/2.0, h + k_1_h/2.0, n + k_1_n/2.0, dt);
        k_3_h = dt * hFRK(V + k_2_V/2.0, m + k_2_m/2.0, h + k_2_h/2.0, n + k_2_n/2.0, dt);
        k_3_n = dt * nFRK(V + k_2_V/2.0, m + k_2_m/2.0, h + k_2_h/2.0, n + k_2_n/2.0, dt);
        k_3_m = dt * mFRK(V + k_2_V/2.0, m + k_2_m/2.0, h + k_2_h/2.0, n + k_2_n/2.0, dt);
        k_3_V = dt * VFRK(V + k_2_V/2.0, m + k_2_m/2.0, h + k_2_h/2.0, n + k_2_n/2.0, dt);
        k_4_h = dt * hFRK(V + k_3_V, m + k_3_m, h + k_3_h, n + k_3_n, dt);
        k_4_n = dt * nFRK(V + k_3_V, m + k_3_m, h + k_3_h, n + k_3_n, dt);
        k_4_m = dt * mFRK(V + k_3_V, m + k_3_m, h + k_3_h, n + k_3_n, dt);
        k_4_V = dt * VFRK(V + k_3_V, m + k_3_m, h + k_3_h, n + k_3_n, dt);

        h += (k_1_h + 2.0 * k_2_h + 2.0 * k_3_h + k_4_h) / 6.0;
        n += (k_1_n + 2.0 * k_2_n + 2.0 * k_3_n + k_4_n) / 6.0;
        m += (k_1_m + 2.0 * k_2_m + 2.0 * k_3_m + k_4_m) / 6.0;
        V += (k_1_V + 2.0 * k_2_V + 2.0 * k_3_V + k_4_V) / 6.0;

        if(V>=Vth && last_spiked+tau_spike < t){
            last_spiked = t;
            return 1;
        }
        return 0;
    }
}

int NeuronHodgkinHuxleyRK::evolve(double dt, double time){
    Ndt = 0;
    subEvolve(time, dt);
    I = 0;

    if(Ndt>0)
        return 1;
    return 0;
}

double* NeuronHodgkinHuxleyRK::exportData(){
    double* arr = new double[numEssentialVariables()+1];
    arr[0] = numEssentialVariables()+1;
    arr[1] = V;
    arr[2] = Vth;
    arr[3] = Vrest;
    arr[4] = tau_spike;
    arr[5] = g_Na;
    arr[6] = g_K;
    arr[7] = g_L;
    arr[8] = E_Na;
    arr[9] = E_K;
    arr[10] = E_L;
    arr[11] = C_mem;
    arr[12] = n;
    arr[13] = m;
    arr[14] = h;
    arr[15] = divider;
    arr[16] = exc;
    return arr;
}

int NeuronHodgkinHuxleyRK::importData(double *arr){
    V = arr[0];
    Vth = arr[1];
    Vrest = arr[2];
    tau_spike = arr[3];
    g_Na = arr[4];
    g_K = arr[5];
    g_L = arr[6];
    E_Na = arr[7];
    E_K = arr[8];
    E_L = arr[9];
    C_mem = arr[10];
    n = arr[11];
    m = arr[12];
    h = arr[13];
    divider = arr[14];
    exc = arr[15];
    return 0;
}

int NeuronHodgkinHuxleyRK::numEssentialVariables(){
    return 16;
}

/// Prototype of new neuron class. Do not alter. Copy, then change.
std::string NeuronPrototype::neurotype = "Insert_name_here";

std::string NeuronPrototype::getName(){
    return neurotype;
}

void NeuronPrototype::setExcitatory(int f){
    I = 0.0;
    if(f){
        exc = 1;
    } else {
        exc = 0;
    }
}

int NeuronPrototype::evolve(double dt, double time){
    // Evolve parameters of your neurons somehow.
    // If the neuron spikes, return 1, else return 0;
    // Befor returning, set I to 0

    I = 0.0;
    return 0;
}
