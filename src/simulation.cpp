#include "simulation.h"

SimulationSingleton* SimulationSingleton::self=NULL;

SimulationSingleton* SimulationSingleton::instance(){
    if (!self) self = new SimulationSingleton;
    return self;
}

SimulationSingleton::SimulationSingleton(){
    // DO NOT CHANGE!
    // Reset default parameters in loadDefaultParameters() function
    neurons_in_simulation = 100;
    N = neurons_in_simulation;
    amount_of_inh_neurons = 0.2 * neurons_in_simulation;
    probability_of_connection = 0.1 * 500 / neurons_in_simulation;
    length_of_simulation = 1000;
    dt = 0.1;
    time_between_exports = 1;
    type_of_delay = 2;
    type_of_neuron = 0;
    type_of_synapse = 7;
    type_of_topology = 0;
    smw_beta = 0.2;
    smw_local = 3;
    syn_noise_freq_mean = 0;
    tau_stim = 0;
    time_between_weight_exports = 100;
    import_neurons = 0;
    import_synapses = 0;
    stim_start = 0;
    time_between_vi_exports = 1000;

    // System variables
    error_number = 0;
    time_now = 0.0;
}

int SimulationSingleton::errorReport(){
    // reports an error to terminal and returns it's number.
    //otherwise returns 0;
    switch (error_number){
    case 0:
        break;
    case 1:
        cout<<"\nError: Invalid input.\n";
        exit(error_number);
        break;
    case 2:
        cout<<"\nError: Problems opening a storage file.\n";
        exit(error_number);
        break;
    case 3:
        cout<<"\nError: Synapse type undefined.\n";
        exit(error_number);
        break;
    case 4:
        cout<<"\nError: Neuron type undefined.\n";
        exit(error_number);
        break;
    case 5:
        cout<<"\nError: Computational error.\n";
        exit(error_number);
        break;
    case 6:
        cout<<"\nError: Wrong topology.\n";
        exit(error_number);
        break;
    case 7:
        cout<<"\nError: Error opening ini file.\n";
        exit(error_number);
        break;
    default:
        cout<<"\nError: Unknown error!\n";
        exit(error_number);
    }

    return error_number;
}

int SimulationSingleton::createNeurons(){
    switch(type_of_neuron){
    case 0:
        for(int i=0; i<N; i++)
            neuron_array[i] = new NeuronLeakyIAF;
        break;
    case 1:
        for(int i=0; i<N; i++)
            neuron_array[i] = new NeuronHodgkinHuxleyRK;
        break;
    default:
        error_number = 4;
    }

    neuron_array[0]->initNeurons();

    for(int i=0; i<N; i++){
        if(i<amount_of_inh_neurons)
            neuron_array[i]->setExcitatory(0);
        else
            neuron_array[i]->setExcitatory(1);
    }

    if(import_neurons){
        importNeurons();
    }

    return 0;
}

int SimulationSingleton::createSynapses(){
    switch (type_of_synapse){
    case 0:
        for(int i=0; i<N*N; i++)
            synapse_array[i] = new SynapseStatic;
        break;
    case 1:
        for(int i=0; i<N*N; i++)
            synapse_array[i] = new SynapseTsodyksMarkramRK;
        break;
    case 2:
        for(int i=0; i<N*N; i++)
            synapse_array[i] = new SynapseGFirstType;
        break;
    case 3:
        for(int i=0; i<N*N; i++)
            synapse_array[i] = new SynapseGSecondType;
        break;
    case 4:
        for(int i=0; i<N*N; i++)
            synapse_array[i] = new SynapseGFirstTypeWCUT;
        break;
    case 5:
        for(int i=0; i<N*N; i++)
            synapse_array[i] = new SynapseGSecondTypeWCUT;
        break;
    case 6:
        for(int i=0; i<N*N; i++)
            synapse_array[i] = new SynapseSTDPG;
        break;
    case 7:
        for(int i=0; i<N*N; i++)
            synapse_array[i] = new SynapseTMSTDP;
        break;
    case 8:
        for(int i=0; i<N*N; i++)
            synapse_array[i] = new SynapseTMSTDPAsymmetrical;
        break;
    default:
        error_number = 4;
    }

    synapse_array[0]->initSynapses();

    /// IF WE IMPORT SYNAPSES, WE DO NOT CREATE ANYTHING ELSE HERE
    if(import_synapses){
        importSynapses();
        return 0;
    }

    /// Creating TOPOLOGY (CONNECTIONS and DELAYS)

    switch(type_of_topology){
    case 0:
        Topology::randomConnections(N, connectivity_matrix, \
                                    probability_of_connection);
        break;
    case 1:
        Topology::smallWorld(N, connectivity_matrix, smw_local, smw_beta);
        break;
    default:
        error_number = 6;
    }

    buf0 = 0;  // for running through connectivity_matrix
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            buf0 = i * N + j;
            if(connectivity_matrix[buf0]){ //if topology said that conn exists
                connectivity_matrix[buf0] = buf0;    //connection create
                synapse_array[buf0]->setData(i, j,
                                            neuron_array[i]->exc, \
                                            neuron_array[j]->exc, dt);
                synapse_array[buf0]->delay = Topology::setDelay(\
                            neuron_array[i]->x, neuron_array[i]->y, \
                            neuron_array[j]->x, neuron_array[j]->y, \
                            type_of_delay, spike_velocity, \
                            delay_max, dt);
                synapse_array[buf0]->setDeliveries(dt);
            } else {
                connectivity_matrix[buf0] = -1; // means NO CONNECTION
            }
        }
    }

    return 0;
}

int SimulationSingleton::createStimulation(){
    /// NEURON STIMULATION:
    switch (type_of_stimulation){
    case 0:
        for(int i=0; i<N; i++){
            Inoise[i] = VFDistributions::uniform(Imin, Imax);
            Inoise2[i] = Inoise[i];
        }
        break;
    case 1:
        for(int i=0; i<N; i++){
            Inoise[i] = VFDistributions::normal(Imean, Isd, Imin, Imax);
            Inoise2[i] = Inoise[i];
        }
        break;
    default:
        error_number = 1;
    }
    for(int i=0; i<N; i++){
        stim_start_pers[i] = VFDistributions::uniform(0, stim_start);
    }

    /// SYNAPSE STIMUALTION:
    // Calculating frequence (per second):
    switch (type_of_synapse){
    case 2:
    case 3:
    case 4:
    case 5:
//        for(int i=0; i< neurons_in_simulation*neurons_in_simulation; i++){
//            if(connectivity_matrix[i] > -1){
//                syn_noise_period[i] = VFDistributions::normal(\
//                            syn_noise_freq_mean, syn_noise_freq_mean,\
//                            syn_noise_freq_mean/10, syn_noise_freq_mean*2);
//            } else {
//                syn_noise_period[i] = 0;
//            }
//        }
        break;
    default:
        for(int i=0; i < N*N; i++)
            syn_noise_period[i] = 0;
    }
    // Calculating period:
    for(int i=0; i < N*N; i++)
        if(syn_noise_period[i]>0){
            syn_noise_period[i] = 1000 / syn_noise_period[i];
        }
}

void SimulationSingleton::createNetwork(){
    neuron_array = Malloc(N, Neuron*);
    synapse_array = Malloc(N*N, Synapse*);
    connectivity_matrix = new int[N*N];
    Inoise = new double[N];
    Inoise2 = new double[N];
    syn_noise_period = new double[N*N];
    stim_start_pers = new double[N];

    createStimulation();
    createNeurons();
    createSynapses();


    /// Set BUFFER for SPIKES:
    buf0 = (floor(time_between_exports/dt)+2)*N;
    spikes_time = new double[buf0];
    spikes_numbers = new int[buf0];
    spikes_resent = 0;
}

void SimulationSingleton::sendNeuralNoise(){
    // Adds noise current to all neurons according to ini

    if (tau_stim > 0)
        if (VFDiscrete::inBetween(time_now, tau_stim, dt)){
            for(int i=0; i<N; i++)
                Inoise2[i] = VFDistributions::normal(Inoise[i], Inoise[i]/10, \
                                                     Imin, Imax);
        }

    for(int i=0; i<N; i++){
        if(time_now>stim_start_pers[i])
            neuron_array[i]->addCurrent(Inoise2[i]);
    }
}

void SimulationSingleton::sendSynapseNoise(){
    // excitates some synapses
    for(int i=0; i<N*N; i++)
        if(syn_noise_period[i]>0){
            if(VFDiscrete::inBetween(time_now, syn_noise_period[i], dt)){
                synapse_array[i]->incSpike(time_now);
            }
        }
}

int SimulationSingleton::evolveAllNeurons(){
    for(int i=0; i < N; i++){
        if((neuron_array[i])->evolve(dt, time_now)){
            // if neuron sent a spike
            saveSpike(i);

            // send signal to POSTsynapses
            for(int j=0; j < N; j++){
                //check all possible connections from this neuron
                if(connectivity_matrix[i*N + j] > -1){
                    //if a synapse exists - send it a message about spike
                    synapse_array[i*N + j]->incSpike(time_now);
                }
            }

            // send signal to PREsynapses
            for(int j=0; j < N; j++){
                //check all possible connections to this neuron
                if(connectivity_matrix[j*N + i] > -1){
                    //if a synapse exists - send it a message about post-spike
                    synapse_array[j*N + i]->incAfterSpike(time_now);
                }
            }
        }
    }

    return 0;
}

int SimulationSingleton::evolveAllSynapses(){
    for(int i=0; i < N; i++){
        for(int j=0; j<N; j++){
            if(connectivity_matrix[i*N + j] > -1){
                // if synapse exists - evolve it
                buf1 = synapse_array[i*N + j]->evolve(dt, time_now, \
                                                      neuron_array[i]->V, \
                                                      neuron_array[j]->V);
                if(buf1 > 1e-030)
                    neuron_array[j]->addCurrent(buf1);
            }
        }
    }

    return 0;
}

void SimulationSingleton::saveSpike(int n){
    spikes_time[spikes_resent] = time_now;
    spikes_numbers[spikes_resent] = n;
    spikes_resent++;
}
