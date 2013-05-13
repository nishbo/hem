#include "simulation.h"

SimulationSingleton* SimulationSingleton::self=NULL;

SimulationSingleton* SimulationSingleton::instance(){
    if (!self) self = new SimulationSingleton;
    return self;
}

SimulationSingleton::SimulationSingleton(){
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

int SimulationSingleton::createNamedNeurons(){
    neuron_array = Malloc(N, Neuron*);

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
        return error_number;
    }

    return 0;
}

int SimulationSingleton::createNeurons(){
    if(import_neurons)
        cout<<"\nImporting neurons...\r";
    else
        cout<<"\nCreating neurons...\r";

    createNamedNeurons();
    neuron_array[0]->initNeurons();

    if(import_neurons)
        importNeurons();
    else
        for(int i=0; i<N; i++){
            if(i<amount_of_inh_neurons)
                neuron_array[i]->setExcitatory(0);
            else
                neuron_array[i]->setExcitatory(1);
            neuron_array[i]->setCoordinates(x[i], y[i]);
        }

    if(import_neurons)
        cout<<N<<" Neurons imported.       ";
    else
        cout<<N<<" Neurons created.       ";

    return 0;
}

int SimulationSingleton::createNamedSynapses(){
    if(M<2){
        M = 1;
        synapse_array = Malloc(1, Synapse*);
    } else
        synapse_array = Malloc(M, Synapse*);

    switch (type_of_synapse){
    case 0:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseStatic;
        break;
    case 1:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseTsodyksMarkramRK;
        break;
    case 2:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseSTDPG;
        break;
    case 3:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseTMSTDP;
        break;
    case 4:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseTsodyksMarkramRKNest;
        break;
    case 5:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseTMexcSTDP;
        break;
    case 51:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseGFirstType;
        break;
    case 52:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseGSecondType;
        break;
    case 53:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseGFirstTypeWCUT;
        break;
    case 54:
        for(int i=0; i<M; i++)
            synapse_array[i] = new SynapseGSecondTypeWCUT;
        break;
    default:
        error_number = 4;
        return error_number;
    }
    return 0;
}

int SimulationSingleton::createSynapseLists(){
    outgoing_synapses = Malloc(N, int*);
    incoming_synapses = Malloc(N, int*);
    outgoing_synapses_to = Malloc(N, int*);
    incoming_synapses_from = Malloc(N, int*);

    for(int i=0; i<N; i++){
        buf0 = 0;
        buf00 = 0;
        for(int j=0; j<N; j++){
            if(connectivity_matrix[i*N+j]>-1)
                buf0++;
            if(connectivity_matrix[j*N+i]>-1)
                buf00++;
        }
        outgoing_synapses[i] = new int[buf0 + 1];
        outgoing_synapses[i][0] = buf0;
        outgoing_synapses_to[i] = new int[buf0 + 1];
        outgoing_synapses_to[i][0] = buf0;
        incoming_synapses[i] = new int[buf00 + 1];
        incoming_synapses[i][0] = buf00;
        incoming_synapses_from[i] = new int[buf00 + 1];
        incoming_synapses_from[i][0] = buf00;
    }

    for(int i=0; i<N; i++){
        buf0 = 1;
        buf00 = 1;
        for(int j=0; j<N; j++){
            if(connectivity_matrix[i*N+j]>-1){
                outgoing_synapses[i][buf0] = connectivity_matrix[i*N+j];
                outgoing_synapses_to[i][buf0] = j;
                buf0++;
            }
            if(connectivity_matrix[j*N+i]>-1){
                incoming_synapses[i][buf00] = connectivity_matrix[j*N+i];
                incoming_synapses_from[i][buf00] = j;
                buf00++;
            }
        }
    }

    free(connectivity_matrix);

    return 0;
}

int SimulationSingleton::createSynapses(){
    int *arr = new int[N];
    for(int i=0; i < N; i++){
        arr[i] = 0;
    }

    createNamedSynapses();
    synapse_array[0]->initSynapses();

    /// Filling synapses with data
    buf00 = 0;  //free synapse to create;
    for(int i=0; i < N; i++){
        outgoing_synapses[i] = new int[outgoing_synapses_to[i][0]+1];
        outgoing_synapses[i][0] = outgoing_synapses_to[i][0];
        for(int j=0; j < outgoing_synapses_to[i][0]; j++){
            outgoing_synapses[i][j+1] = buf00;
            synapse_array[buf00]->setData(i, outgoing_synapses_to[i][j+1], \
                neuron_array[i]->exc, \
                neuron_array[outgoing_synapses_to[i][j+1]]->exc, dt);
            synapse_array[buf00]->delay = delays[buf00];
            synapse_array[buf00]->setDeliveries(dt);

            arr[outgoing_synapses_to[i][j+1]]++;

            buf00++;
        }
    }

    for(int i=0; i < N; i++){
        incoming_synapses[i] = new int[arr[i]+1];
        incoming_synapses[i][0] = 0;
        incoming_synapses_from[i] = new int[arr[i]+1];
        incoming_synapses_from[i][0] = 0;
    }
    for(int i=0; i < N; i++){
        for(int j=0; j < outgoing_synapses[i][0]; j++){
            buf0 = outgoing_synapses_to[i][j+1];
            incoming_synapses_from[buf0][++incoming_synapses_from[buf0][0]] = i;
            incoming_synapses[buf0][++incoming_synapses[buf0][0]] = outgoing_synapses[i][j+1];
        }
    }

    free(arr);
    return 0;
}

int SimulationSingleton::createNeuronStimulation(){
    using namespace vf_distributions;

    switch (type_of_stimulation){
    case 0:
        for(int i=0; i<N; i++){
            Inoise[i] = uniform(Imin, Imax);
            Inoise2[i] = Inoise[i];
        }
        break;
    case 1:
        for(int i=0; i<N; i++){
            Inoise[i] = normal(Imean, Isd, Imin, Imax);
            Inoise2[i] = Inoise[i];
        }
        break;
    case 2:
        for(int i=0; i<N; i++){
            Inoise[i] = normal(Imean, Isd, Imin, Imax);
            if(Inoise[i]<15)
                Inoise[i] = 0;
            Inoise2[i] = Inoise[i];
        }
        break;
    default:
        error_number = 1;
    }

    for(int i=0; i<N; i++){
        stim_start_pers[i] = uniform(0, stim_start);
    }

    return 0;
}

int SimulationSingleton::createSynapseStimulation(){
    // Calculating frequence (per second):
    for(int i=0; i < M; i++)
        syn_noise_period[i] = vf_distributions::normal(\
                    syn_noise_freq_mean, syn_noise_freq_mean,\
                    syn_noise_freq_mean/10, syn_noise_freq_mean*2);
    
    // Calculating period:
    for(int i=0; i < M; i++)
        if(syn_noise_period[i]>0)
            syn_noise_period[i] = 1000 / syn_noise_period[i];

    return 0;
}

int SimulationSingleton::createStimulation(){
    cout<<"\nCreating stimulation...\r";

    if(!import_neurons)
        createNeuronStimulation();
    createSynapseStimulation();

    cout<<"Stimulation created.   ";
    return errorReport();
}

int SimulationSingleton::createTopology(){
    x = new double[N];
    y = new double[N];
    outgoing_synapses_to = Malloc(N, int*);
    outgoing_synapses = Malloc(N, int*);
    incoming_synapses = Malloc(N, int*);
    incoming_synapses_from = Malloc(N, int*);

    int *Ms = Malloc(1, int);
    Topology::setTopology(N, Ms, x, y, outgoing_synapses_to, &delays);
    M = *Ms;
    free(Ms);

    return 0;
}

int SimulationSingleton::createNetwork(){
    loadParameters();
    createTopology();

    Inoise = new double[N];
    Inoise2 = new double[N];
    stim_start_pers = new double[N];
    createNeurons();
    if(import_synapses){
        cout<<"\nImporting synapses...\r";
        importSynapses();
        createSynapseLists();
        cout<<M<<" Synapses imported.     ";
    } else {
        cout<<"\nCreating synapses...\r";
        createSynapses();
        cout<<M<<" Synapses created.     ";
    }
    syn_noise_period = new double[M];
    createStimulation();

    /// Set BUFFER for SPIKES:
    buf0 = (floor(time_between_exports/dt)+2)*N;
    spikes_time = new double[buf0];
    spikes_numbers = new int[buf0];
    spikes_resent = 0;

    return errorReport();
}

void SimulationSingleton::sendNeuralNoise(){
    // Adds noise current to all neurons according to ini

    if (tau_stim > 0)
        if (vf_discrete::inBetween(time_now, tau_stim, dt)){
            for(int i=0; i<N; i++)
                Inoise2[i] = vf_distributions::normal(Inoise[i], Inoise[i]/10, \
                                                     Imin, Imax);
        }

    for(int i=0; i<N; i++){
        if(time_now>stim_start_pers[i])
            neuron_array[i]->addCurrent(Inoise2[i]);
    }
}

void SimulationSingleton::sendSynapseNoise(){
    // excites some synapses
    for(int i=0; i<M; i++)
        if(syn_noise_period[i] > 0)
            if(vf_discrete::inBetween(time_now, syn_noise_period[i], dt))
                synapse_array[i]->incSpike(time_now);
}

int SimulationSingleton::evolveAllNeurons(){
    for(int i=0; i < N; i++){
        if((neuron_array[i])->evolve(dt, time_now)){
            // if neuron sent a spike
            saveSpike(i);

            // send signal to POSTsynapses
            for(int j=0; j < outgoing_synapses[i][0]; j++)
                synapse_array[outgoing_synapses[i][j+1]]->incSpike(time_now);

            // send signal to PREsynapses
            for(int j=0; j < incoming_synapses[i][0]; j++)
                synapse_array[incoming_synapses[i][j+1]]->incAfterSpike(time_now);
        }
    }

    return 0;
}

int SimulationSingleton::evolveAllSynapses(){
    for(int i=0; i < M; i++){
        neuron_array[synapse_array[i]->to()]->addCurrent(\
            synapse_array[i]->evolve(dt, time_now, \
                neuron_array[synapse_array[i]->from()]->V, \
                neuron_array[synapse_array[i]->to()]->V));
    }

    return 0;
}

void SimulationSingleton::saveSpike(int n){
    spikes_time[spikes_resent] = time_now;
    spikes_numbers[spikes_resent] = n;
    spikes_resent++;
}
