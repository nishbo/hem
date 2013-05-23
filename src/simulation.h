#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>

#include "cad.h"
#include "neuron.h"
#include "synapse.h"
#include "topology.h"
#include "vfdistributions.h"
#include "vfdiscrete.h"
#include "vffile.h"

using namespace std;

class SimulationSingleton{
public:
    /// Parameters.
    // Main parameters:
    int neurons_in_simulation;      // total amount of neurons in simulation
    int N;                          // copy of neurons_in_simulation
    int amount_of_inh_neurons;
    double length_of_simulation;    //in msec
    double dt;                      //in msec
    double time_now;                // what time is it?

    Neuron** neuron_array;      //array of all neurons
    Synapse** synapse_array;    //array of all synapses
    int* connectivity_matrix; //it is just an array. new line every N elements.
    int M;                      //total number of synapses
    int** outgoing_synapses;    //lists of outgoing synapses per neuron
    int** incoming_synapses;    //lists of incoming synapses per neuron
    int** outgoing_synapses_to;
    int** incoming_synapses_from;

    //System parameters:
    int number_of_threads;
    int enable_test;                //if 1 = testing

    // Type parameters:
    int type_of_neuron;
    int type_of_synapse;
    int type_of_stimulation;//0 for uniform, 1 for normal

    // Stimulation parameters:
    double* Inoise;         //Stimuli current
    double* Inoise2;        //Stimuli current buffer
    double Imax, Imin, Imean, Isd;
    double tau_stim;        //every tau_stim stimulation current changes
    double* syn_noise_period; //period of noise excitating of synapses (in msec)
    int syn_noise_freq_mean;  //mean in poisson distr
    double stim_start;      //time of start of stimulation
    double* stim_start_pers;//personal times of start of stimulation

    // Output parameters:
    double time_between_weight_exports;
    double time_between_exports;    //in msec
    double time_between_vi_exports; //V, I, synapse
    FILE *storage_file;
    FILE *spike_file;
    FILE *synapse_file;
    FILE *weight_file;
    FILE *current_file;
    std::string one_line;  //less accessing the file
    double *spikes_time;
    int *spikes_numbers;
    int spikes_resent;
    double *buffer;

    // Import parameters:
    int import_neurons;
    int import_synapses;

    // Buffer and other parameters:
    int buf0, buf00, buf000;
    float buf01;
    double buf1, buf2;
    string buf31;
    char* buf30;
    int error_number;     //number of error that occured. 0 if OK

    //Topology getting:
    double *x, *y, *delays;

    /// Functions.
    static SimulationSingleton* instance();
    // Creating network:
    int createNetwork();

    int loadParameters();
    int loadParametersFromFile();

    int createTopology();

    int createStimulation();
    int createSynapseStimulation();
    int createNeuronStimulation();

    int importNeurons();
    int createNeurons();
    int createNamedNeurons();

    int importSynapses();
    int createSynapses();
    int createNamedSynapses();
    int createSynapseLists();


    int errorReport();
    // Simulating:
    void sendNeuralNoise();
    void sendSynapseNoise();
    void saveSpike(int n);
    int evolveAllNeurons();
    int evolveAllSynapses();

    // Output, etc. In "inout.cpp" are bodies.
    int setOutputFile();
    void closeOutputFile();
    // export data to file
    void outputChangingData();
    int outputParametersInFile();      //delays, Vrest, Vth
    int outputConnectivityMatrixInFile();
    void outputChangingDataInFile(double time);      //potentials, weights
    void outputSynapseDataInFile(double time);
    void outputSpikesInFile();
    int outputWeightsInFile(double time);
    int outputCurrentsInFile(double time);
    int exportNeurons();
    int exportSynapses();

    void test();
private:
    static SimulationSingleton* self;
    SimulationSingleton();
    SimulationSingleton(SimulationSingleton&p){}
    SimulationSingleton& operator=(SimulationSingleton&){}
};

#endif // SIMULATION_H
