#ifndef NEURON_H
#define NEURON_H

#include <iostream>
#include <cstdio>

#include "vfdiscrete.h"
#include "vfdistributions.h"
#include "vffile.h"
#include "cad.h"

class Neuron{
public:
    /// Commonly used variables:
    double V;           //current potential, mV
    double Vth;         //threshold potential, mV
    double Vrest;       //rest potential, mV
    double I;           //current, nA
    double x;           //coordinates
    double y;
    double last_spiked; //time of last spike
    int exc;            //1 if excitatory, 0 if inhibitory

    /// MUST overwrite functions:
    // Evolves a neuron for "dt" in time "time"
    virtual int evolve(double dt, double time) = 0;
    // Sets as excitatory and (optionally) some other parameters
    virtual void setExcitatory(int f) = 0;
    /// OPTIONALLY overwrite functions:
    // Specifies how the neuron accepts current
    virtual void addCurrent(double a);
    // Gets the name for this type of neurons:
    virtual std::string getName();
    // Export all data to reproduce.
    virtual double* exportData();
    // Import data from file.
    virtual int importData(double *arr);
    virtual int initNeurons();
    virtual int numEssentialVariables();
    virtual int excite(double time);
    int setCoordinates(double xo, double yo);
};

class NeuronLeakyIAF: public Neuron{
    // Euler-style
public:
    static std::string neurotype;
    double Rin;         //membrance resistence Ohm
    double tau_m;       //membrane time constant msec
    double tau_ref_abs; //absolute refractory period
    double Vreset;      //for i-a-f to work

    static double init_Rin, init_tau_m, init_tau_ref_abs_exc;
    static double init_tau_ref_abs_inh, init_Vreset;
    static double init_V, init_Vth, init_Vrest;

    /// Overwriting abstract functions:
    void setExcitatory(int f);
    int evolve(double dt, double time);
    std::string getName();
    double* exportData();
    int importData(double *arr);
    int initNeurons();
    static int initNeuronsLocal();
    int numEssentialVariables();
};

class NeuronHodgkinHuxleyRK: public Neuron{
    // RK4 style
public:
    static std::string neurotype;
    double a, b;  //buff used in m,n,h func
    double m, n, h;
    double g_Na, g_K, g_L;
    double E_Na, E_K, E_L;
    double k_1_n, k_1_m, k_1_h, k_1_V;
    double k_2_n, k_2_m, k_2_h, k_2_V;
    double k_3_n, k_3_m, k_3_h, k_3_V;
    double k_4_n, k_4_m, k_4_h, k_4_V;
    double C_mem;
    double tau_spike;
    int Ndt;
    double divider;

    double mFRK(double V, double m, double h, double n, double dt);
    double hFRK(double V, double m, double h, double n, double dt);
    double nFRK(double V, double m, double h, double n, double dt);
    double VFRK(double V, double m, double h, double n, double dt);
    int subEvolve(double t, double dt);

    static double init_tau_spike, init_C_mem, init_g_Na, init_g_K, init_g_L;
    static double init_E_Na, init_E_K, init_E_L, init_Vreset;
    static double init_n, init_m, init_h;
    static double init_V, init_Vth, init_Vrest, init_divider;

    /// Overwriting abstract functions:
    void setExcitatory(int f);
    void addCurrent(double a);
    int evolve(double dt, double time);
    std::string getName();
    int initNeurons();
    static int initNeuronsLocal();
    double* exportData();
    int importData(double *arr);
    int numEssentialVariables();
};

#endif // NEURON_H
