#ifndef SYNAPSE_H
#define SYNAPSE_H

#include <iostream>
#include <cstdio>

#include "cad.h"
#include "vfdistributions.h"
#include "vfdiscrete.h"

class Synapse{
public:
    /// Commonly used variables
    double weight;      //0..1
    double delay;       //in msec
    int presynaptic;    //number of post- and presynaptic neuron in array
    int postsynaptic;       // in storage
    double last_spiked;     //time of last spike of synapse(presyn neuron)
    double last_spiked_post;    //time of last spike of postsynaptic neuron

    /// Tech variables:
    int working;        //1 if exists, 0 if not
    double out_current;     //stores current to send
    double* delivery;
    int number_of_possible_transfers;
    double buf;
    static double* innerDataArr;

    /// Everyones' functions
    void incSpike(double timen);
    void incAfterSpike(double timen);
    double moveDeliveries();
    void setDeliveries(double dt);
    /// OPTIONALLY overwritable functions
    virtual double* getInnerData();
    virtual std::string getName();
    virtual int importData(double *arr);
    virtual double* exportData();
    virtual int initSynapses();
    virtual int numEssentialVariables();
    /// MUST overwritable functions
    // Evolve synapse for a short period of time and set out_current that will
    // be delivered after delay to neuron
    virtual double evolve(double dt, double time, double Vpre, double Vpost) = 0;
    // Set some data-parameters to synapse
    virtual void setData (int pre, int pos, int preex, int posex, double dt) = 0;
};

class SynapseStatic:  public Synapse{
    // Nothing changes
public:
    static std::string synapsetype;

    static double init_weight;

    SynapseStatic ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
    double* exportData();
    int initSynapses();
    static int initSynapsesLocal();
};

class SynapseTsodyksMarkram: public Synapse{
    // Tsodyks-Markram euler style
public:
    static std::string synapsetype;
    double x,y,z,u;
    double xo, yo, zo, uo, dd;  //buffer
    double tau_one, tau_rec, tau_facil, U, A;
    int exc;
    static double xav, yav, zav, uav;
    static int number_of_synapses;

    SynapseTsodyksMarkram ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
    double* getInnerData();
};

class SynapseTsodyksMarkramRK: public Synapse{
    // Tsodyks-Markram rk4 style
public:
    static std::string synapsetype;
    double x, y, z, u;
    double tau_one, tau_rec, tau_facil, U, A;
    int exc;
    // Runge-Kutta stuff:
    double k_1_x, k_1_y, k_1_z, k_1_u;
    double k_2_x, k_2_y, k_2_z, k_2_u;
    double k_3_x, k_3_y, k_3_z, k_3_u;
    double k_4_x, k_4_y, k_4_z, k_4_u;
    double Xr(double x1, double y1, double z1, double u1, double t1, double dt);
    double Yr(double x1, double y1, double z1, double u1, double t1, double dt);
    double Zr(double x1, double y1, double z1, double u1, double t1, double dt);
    double Ur(double x1, double y1, double z1, double u1, double t1, double dt);

    static double xav, yav, zav, uav;
    static int number_of_synapses;

    static double init_tau_one, init_x, init_y, init_z, init_u;
    static double init_tau_recee, init_tau_facilee, init_Uee, init_Aee;
    static double init_tau_recei, init_tau_facilei, init_Uei, init_Aei;
    static double init_tau_recie, init_tau_facilie, init_Uie, init_Aie;
    static double init_tau_recii, init_tau_facilii, init_Uii, init_Aii;

    SynapseTsodyksMarkramRK ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
    double* getInnerData();
    double* exportData();
    int initSynapses();
    static int initSynapsesLocal();
    int importData(double *arr);
    int numEssentialVariables();
};

class SynapseTsodyksMarkramAn: public Synapse{
    // Tsodyks-Markram analytical
public:
    static std::string synapsetype;
    double x, y, z, u;
    double xsp, ysp, zsp, usp;
    double tau_one, tau_rec, tau_facil, U, A;
    double efacil, eone, erec;
    int exc;
    int buf00;

    static double init_tau_one, init_x, init_y, init_z, init_u;
    static double init_tau_recee, init_tau_facilee, init_Uee, init_Aee;
    static double init_tau_recei, init_tau_facilei, init_Uei, init_Aei;
    static double init_tau_recie, init_tau_facilie, init_Uie, init_Aie;
    static double init_tau_recii, init_tau_facilii, init_Uii, init_Aii;

    SynapseTsodyksMarkramAn ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
    double* exportData();
    int initSynapses();
    static int initSynapsesLocal();
    int importData(double *arr);
    int numEssentialVariables();
};

class SynapseTsodyksMarkramRKNest: public Synapse{
    // Tsodyks-Markram rk4 style
public:
    static std::string synapsetype;
    double xo, yo, zo, uo, h;
    double x1, y1, u1, u2;
    double tau_one, tau_rec, tau_facil, U, A;
    double g, tau_g, E;
    double lsp2;
    int exc;

    static double init_tau_one, init_x, init_y, init_z, init_u;
    static double init_tau_recee, init_tau_facilee, init_Uee, init_Aee;
    static double init_tau_recei, init_tau_facilei, init_Uei, init_Aei;
    static double init_tau_recie, init_tau_facilie, init_Uie, init_Aie;
    static double init_tau_recii, init_tau_facilii, init_Uii, init_Aii;

    SynapseTsodyksMarkramRKNest ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
    double* exportData();
    int initSynapses();
    static int initSynapsesLocal();
    int importData(double *arr);
    int numEssentialVariables();
};

class SynapseSTDP: public Synapse{
public:
    static std::string synapsetype;
    double lambda, alpha, tau_corr;

    SynapseSTDP ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int post, int preex, int posex, double dt);
};

class SynapseSTDPG: public Synapse{
    // STDP+G
public:
    static std::string synapsetype;
    double lambda, alpha, tau_corr;
    double g, tau_s;

    static double init_lambda, init_alpha, init_tau_corr, init_g, init_tau_s;
    static double init_weight;

    SynapseSTDPG ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int post, int preex, int posex, double dt);
    double* exportData();
    int initSynapses();
    static int initSynapsesLocal();
    int importData(double *arr);
    int numEssentialVariables();
};

class SynapseTMSTDP: public Synapse{
    // Tsodyks-Markram + STDP
public:
    static std::string synapsetype;
    double lambda, alpha, tau_corr;
    double x, y, z, u;
    double tau_one, tau_rec, tau_facil, U, A;
    double t_start;
    int exc;

    // Runge-Kutta stuff:
    double k_1_x, k_1_y, k_1_z, k_1_u;
    double k_2_x, k_2_y, k_2_z, k_2_u;
    double k_3_x, k_3_y, k_3_z, k_3_u;
    double k_4_x, k_4_y, k_4_z, k_4_u;
    double Xr(double x1, double y1, double z1, double u1, double t1, double dt);
    double Yr(double x1, double y1, double z1, double u1, double t1, double dt);
    double Zr(double x1, double y1, double z1, double u1, double t1, double dt);
    double Ur(double x1, double y1, double z1, double u1, double t1, double dt);

    static double init_tau_one, init_x, init_y, init_z, init_u;
    static double init_tau_recee, init_tau_facilee, init_Uee, init_Aee;
    static double init_tau_recei, init_tau_facilei, init_Uei, init_Aei;
    static double init_tau_recie, init_tau_facilie, init_Uie, init_Aie;
    static double init_tau_recii, init_tau_facilii, init_Uii, init_Aii;
    static double init_lambda, init_alpha, init_tau_corr;
    static double init_weight, init_t_start;
    static int init_type_of_weight;

    SynapseTMSTDP ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int post, int preex, int posex, double dt);
    double* exportData();
    int importData(double *arr);
    int initSynapses();
    static int initSynapsesLocal();
    int numEssentialVariables();
};

class SynapseTMSTDPAsymmetrical: public Synapse{
    // Tsodyks-Markram + STDP (asymmetrical version)
public:
    static std::string synapsetype;
    //STDP
    double A_plus, A_minus, tau_corr_plus, tau_corr_minus;
    double t_start;
    //TM
    double x, y, z, u;
    double tau_one, tau_rec, tau_facil, U, A;
    int exc;

    // Runge-Kutta stuff:
    double k_1_x, k_1_y, k_1_z, k_1_u;
    double k_2_x, k_2_y, k_2_z, k_2_u;
    double k_3_x, k_3_y, k_3_z, k_3_u;
    double k_4_x, k_4_y, k_4_z, k_4_u;
    double Xr(double x1, double y1, double z1, double u1, double t1, double dt);
    double Yr(double x1, double y1, double z1, double u1, double t1, double dt);
    double Zr(double x1, double y1, double z1, double u1, double t1, double dt);
    double Ur(double x1, double y1, double z1, double u1, double t1, double dt);

    static double init_tau_one, init_x, init_y, init_z, init_u;
    static double init_tau_recee, init_tau_facilee, init_Uee, init_Aee;
    static double init_tau_recei, init_tau_facilei, init_Uei, init_Aei;
    static double init_tau_recie, init_tau_facilie, init_Uie, init_Aie;
    static double init_tau_recii, init_tau_facilii, init_Uii, init_Aii;
    static double init_A_plus, init_A_minus;
    static double init_tau_corr_plus, init_tau_corr_minus;
    static double init_weight, init_t_start;
    static int init_type_of_weight;

    SynapseTMSTDPAsymmetrical ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int post, int preex, int posex, double dt);
    double* exportData();
    int importData(double *arr);
    int initSynapses();
    static int initSynapsesLocal();
    int numEssentialVariables();
};

class SynapseGFirstType: public Synapse{
    // Differential eq, no cut
public:
    static std::string synapsetype;
    double g, gs;
    double tau_s;
    double E;   //reverse potential

    static double init_g, init_gs, init_tau_se, init_tau_si, init_Ee, init_Ei;

    SynapseGFirstType ();
    std::string getName();
    double evolve(double dt, double timen, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
    double* exportData();
    int importData(double *arr);
    int initSynapses();
    static int initSynapsesLocal();
    int numEssentialVariables();
};

class SynapseGFirstTypeWCUT: public Synapse{
    // Differential eq, with cut
public:
    static std::string synapsetype;
    double g, gs;
    double tau_s;
    double E;   //reverse potential

    static double init_g, init_gs, init_tau_se, init_tau_si, init_Ee, init_Ei;

    SynapseGFirstTypeWCUT ();
    std::string getName();
    double evolve(double dt, double timen, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
    double* exportData();
    int importData(double *arr);
    int initSynapses();
    static int initSynapsesLocal();
    int numEssentialVariables();
};

class SynapseGSecondType: public Synapse{
    // Linear eq, no cut
public:
    static std::string synapsetype;
    double g, gs;
    double tau_s;
    int max_spikes;
    double *spikerow;
    int now_spikes;
    double E;   //reverse potential

    static double init_g, init_gs, init_tau_se, init_tau_si, init_Ee, init_Ei;

    SynapseGSecondType ();
    std::string getName();
    double alphaFunction(double x);
    double evolve(double dt, double timen, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
    double* exportData();
    int importData(double *arr);
    int initSynapses();
    static int initSynapsesLocal();
    int numEssentialVariables();
};

class SynapseGSecondTypeWCUT: public Synapse{
    // Linear eq, cut
public:
    static std::string synapsetype;
    double g, gs;
    double tau_s;
    int max_spikes;
    double *spikerow;
    int now_spikes;
    double E;   //reverse potential

    static double init_g, init_gs, init_tau_se, init_tau_si, init_Ee, init_Ei;

    SynapseGSecondTypeWCUT ();
    std::string getName();
    double alphaFunction(double x);
    double evolve(double dt, double timen, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
    double* exportData();
    int importData(double *arr);
    int initSynapses();
    static int initSynapsesLocal();
    int numEssentialVariables();
};

class SynapsePrototype:  public Synapse{
    // Use this as a prototype. Do not change, copy, then alter.
public:
    static std::string synapsetype;

    SynapsePrototype ();
    std::string getName();
    double evolve(double dt, double time, double Vpre, double Vpost);
    void setData(int pre, int pos, int preex, int posex, double dt);
};

#endif // SYNAPSE_H
