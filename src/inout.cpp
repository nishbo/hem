#include "simulation.h"

void SimulationSingleton::loadParameters(){
    // checks USE_DEFAULT_PARAMETERS definition and if 0 asks users various
    // questions
    if(USE_DEFAULT_PARAMETERS == 1){
        cout<<"You decided to use default parameters.\n";
        loadDefaultParameters();
    } else if(USE_DEFAULT_PARAMETERS == 0){
        cout<<"Please, set simulation parameters:\n\n";
        cout<<"Function removed from current release."<<\
              " Please use other options.";
        exit(76);
        loadParametersFromTerminal();
    } else if(USE_DEFAULT_PARAMETERS == 2){
        cout<<"You decided to load parameters from file.\n";
        loadParametersFromFile();
    }
}

void SimulationSingleton::loadDefaultParameters(){
    // You can reset default parameters here.
    // DO NOT change them in creation function.

}

void SimulationSingleton::loadParametersFromTerminal(){
//    cout<<"Type of neurons (0 for i-a-f, 1 for h-h, 2 for g):";
//    type_of_neuron = VFFile::checkNumberFromCin();
//    cout<<"Type of synapse (0 for static, 1 for ts-m, 2 for g): ";
//    type_of_synapse = VFFile::checkNumberFromCin();
//    cout<<"Amount of neurons in simulation: ";
//    neurons_in_simulation = VFFile::checkNumberFromCin();
//    cout<<"Amount of inhibitory neurons in simulation: ";
//    amount_of_inh_neurons = VFFile::checkNumberFromCin();
//    cout<<"Time-length of simulation (msec): ";
//    length_of_simulation = VFFile::checkNumberFromCin();
//    cout<<"Choose type of topology. Type \"0\" for random, \"1\" for"<<\
//          " small-world: ";
//    type_of_topology = VFFile::checkNumberFromCin();
//    if(type_of_topology == 0){
//        cout<<"\tProbability of connection between neurons: ";
//        probability_of_connection = VFFile::checkNumberFromCin();
//    } else if(type_of_topology == 1){
//        cout<<"\tAmount of local nodes connections divided by 2 (for 4 conn"<<\
//              " type 2, for 6 type 3): ";
//        smw_local = VFFile::checkNumberFromCin();
//        cout<<"\tProbability of rewriting of connection: ";
//        smw_beta = VFFile::checkNumberFromCin();
//    }
//    cout<<"Choose type of delay. Type \"0\" for random, \"1\" for"<<\
//          " length-dependent: ";
//    type_of_delay = VFFile::checkNumberFromCin();
}

int SimulationSingleton::loadParametersFromFile(){
    FILE* fid = fopen(FILE_IMP_PARAMETERS, "r");
    if(!fid){
        error_number = 7;
        fclose(fid);
        return error_number;
    }

//    string* param_names = new string[16];
//    param_names[0] = "--";
//    param_names[1] = "Amount_of_neurons";
//    param_names[2] = "Amount_of_inhibitory_neurons";
//    param_names[3] = "Probability_of_connection";
//    param_names[4] = "Length_of_simulation";
//    param_names[5] = "dt";
//    param_names[6] = "Time_between_exports";
//    param_names[7] = "Type_of_delay";
//    param_names[8] = "Type_of_neuron";
//    param_names[9] = "Type_of_synapse";
//    param_names[10] = "Type_of_topology";
//    param_names[11] = "smw_beta";
//    param_names[12] = "smw_local";
//    param_names[13] = "Synaptic_noise_mean_frequency";
//    param_names[14] = "Tau_stimulation";
//    param_names[15] = "Time_between_weight_exports";

    buf01 = 0;
    fscanf(fid, "Amount_of_neurons = %d\n", &neurons_in_simulation);
    N = neurons_in_simulation;
    fscanf(fid, "Amount_of_inhibitory_neurons = %d\n", &amount_of_inh_neurons);
    fscanf(fid, "Probability_of_connection = %f\n", &buf01);
    probability_of_connection = buf01;
    fscanf(fid, "Length_of_simulation = %f\n", &buf01);
    length_of_simulation = buf01;
    fscanf(fid, "dt = %f\n", &buf01);
    dt = buf01;
    fscanf(fid, "Time_between_exports = %f\n", &buf01);
    time_between_exports = buf01;
    fscanf(fid, "Time_between_vi_exports = %f\n", &buf01);
    time_between_vi_exports = buf01;
    fscanf(fid, "Type_of_delay = %d\n", &type_of_delay);
    fscanf(fid, "Type_of_neuron = %d\n", &type_of_neuron);
    fscanf(fid, "Type_of_synapse = %d\n", &type_of_synapse);
    fscanf(fid, "Type_of_topology = %d\n", &type_of_topology);
    fscanf(fid, "Type_of_stimulation = %d\n", &type_of_stimulation);
    fscanf(fid, "smw_beta = %f\n", &buf01);
    smw_beta = buf01;
    fscanf(fid, "smw_local = %d\n", &smw_local);
    fscanf(fid, "Synaptic_noise_mean_frequency = %d\n", &syn_noise_freq_mean);
    fscanf(fid, "Tau_stimulation = %f\n", &buf01);
    tau_stim = buf01;
    fscanf(fid, "Start_of_stimulation = %f\n", &buf01);
    stim_start = buf01;
    fscanf(fid, "Time_between_weight_exports = %f\n", &buf01);
    time_between_weight_exports = buf01;
    fscanf(fid, "Min_noise = %f\n", &buf01);
    Imin = buf01;
    fscanf(fid, "Max_noise = %f\n", &buf01);
    Imax = buf01;
    fscanf(fid, "Mean_stimulation = %f\n", &buf01);
    Imean = buf01;
    fscanf(fid, "Sigma_stimulation = %f\n", &buf01);
    Isd = buf01;
    fscanf(fid, "Border_length_of_box = %f\n", &buf01);
    border_length_of_box = buf01;
    fscanf(fid, "Spike_velocity = %f\n", &buf01);
    spike_velocity = buf01;
    fscanf(fid, "Delay_max = %f\n", &buf01);
    delay_max = buf01;
    fscanf(fid, "Import_network_neurons = %d\n", &import_neurons);
    fscanf(fid, "Import_network_synapses = %d\n", &import_synapses);

    fclose(fid);
    return 0;
}

void SimulationSingleton::outputChangingData(){
    if(VFDiscrete::inBetween(time_now, time_between_vi_exports, dt)){
        outputChangingDataInFile(time_now);
        outputSynapseDataInFile(time_now);
        outputCurrentsInFile(time_now);
    }
    if(VFDiscrete::inBetween(time_now, time_between_exports, dt)){
        outputSpikesInFile();
    }
    if(VFDiscrete::inBetween(time_now, time_between_weight_exports, dt)){
        outputWeightsInFile(time_now);
    }
}

int SimulationSingleton::setOutputFile(){
    if(!(VFFile::tryFile(FILE_EXP_MAIN_OUTPUT)) || \
       !(VFFile::tryFile(FILE_EXP_SPIKES)) || \
       !(VFFile::tryFile(FILE_EXP_SYNAPSE)) || \
       !(VFFile::tryFile(FILE_EXP_WEIGHTS)) || \
       !(VFFile::tryFile(FILE_EXP_CURRENT)) ){
        error_number = 2;
        return error_number;
    }

    storage_file = fopen (FILE_EXP_MAIN_OUTPUT, "w");
    spike_file = fopen (FILE_EXP_SPIKES, "w");
    synapse_file = fopen (FILE_EXP_SYNAPSE, "w");
    weight_file = fopen (FILE_EXP_WEIGHTS, "w");
    current_file = fopen (FILE_EXP_CURRENT, "w");

    if(!storage_file || !synapse_file || !spike_file || !current_file){
        error_number = 2;
        return 2;
    }

    return 0;
}

void SimulationSingleton::closeOutputFile(){
    fclose(storage_file);
    fclose(spike_file);
    fclose(synapse_file);
    fclose(weight_file);
    fclose(current_file);
}

int SimulationSingleton::outputParametersInFile(){
    // Saves parameters in parameters.txt
    if(!(VFFile::tryFile(FILE_EXP_PARAMETERS)) || \
       !(VFFile::tryFile(FILE_EXP_NEURON_PARAMS)) ){
        error_number = 2;
        return error_number;
    }
    FILE* param_file = fopen (FILE_EXP_PARAMETERS, "w");
    if(!param_file){
        error_number = 2;
        return 2;
    }

    fprintf(param_file, "________Simulation parameters:\n");
    fprintf(param_file, "Type of neurons _=_ ");
    fprintf(param_file, neuron_array[0]->getName().c_str());
    fprintf(param_file, ";\n");

    fprintf(param_file, "Type of synapses _=_ ");
    fprintf(param_file, synapse_array[0]->getName().c_str());
    fprintf(param_file, ";\n");

    fprintf(param_file, "Neurons in simulation _=_ %d;\n", N);
    fprintf(param_file, "Probability of connection _=_ %.2f;\n", \
                        probability_of_connection);
    fprintf(param_file, "Length of simulation (msec) _=_ %.2f;\n", \
                        length_of_simulation);
    fprintf(param_file, "Time between exports (msec) _=_ %.2f;\n", \
                        time_between_exports);
    fprintf(param_file, "Time-step (msec) _=_ %.4f;\n", \
                        dt);
    fprintf(param_file, "\n");

    fclose(param_file);

    param_file = fopen (FILE_EXP_NEURON_PARAMS, "w");
    if(!param_file){
        error_number = 2;
        return 2;
    }

    fprintf(param_file, "________Neurons potential threshold:\n");
    buf31 = "";
    for(int i=0; i<N; i++)
        buf31+= VFFile::convertDoubleToString(neuron_array[i]->Vth) + " ";
    fprintf(param_file, buf31.c_str());
    fprintf(param_file, "\n");

    fprintf(param_file, "________Neurons Vrest:\n");
    buf31 = "";
    for(int i=0; i<N; i++)
        buf31+= VFFile::convertDoubleToString(neuron_array[i]->Vrest) + " ";
    fprintf(param_file, buf31.c_str());
    fprintf(param_file, "\n");

    buf0 = 0;
    fprintf(param_file, "________Synapse delays:\n");
    for(int i=0; i<N; i++){
        buf31 = "";
        for(int j=0; j<N; j++){
            buf0 = i * N + j;
            if(this->connectivity_matrix[buf0]>-1)
                buf1 = synapse_array[buf0]->delay;
            else
                buf1 = 0;
            buf1 = buf1 - (buf1 * 100 - (int)(buf1 * 100))/100;
            buf31 += VFFile::convertDoubleToString(buf1) + " ";
        }
        fprintf(param_file, buf31.c_str());
        fprintf(param_file, "\n");
    }
    fclose(param_file);
    return 0;
}

int SimulationSingleton::outputConnectivityMatrixInFile(){
    // Saves parameters in parameters.txt
    if(!(VFFile::tryFile(FILE_EXP_CONN_MATR))){
        error_number = 2;
        return error_number;
    }
    FILE *param_file = fopen (FILE_EXP_CONN_MATR, "w");
    if(!param_file){
        error_number = 2;
        return error_number;
    }

    fprintf(param_file, "________Connectivity matrix:\n");
    for(int i=0; i<N; i++){
        buf31 = "";
        for(int j=0; j<N; j++){
            buf0 = i * N + j;
            if(connectivity_matrix[buf0]>-1)
                buf31+= "+ ";
            else
                buf31 += "- ";
        }
        fprintf(param_file, buf31.c_str());
        fprintf(param_file, "\n");
    }
    fclose(param_file);
    return 0;
}

void SimulationSingleton::outputChangingDataInFile(double time){
    fprintf(storage_file, "_________Time = %.2f\n", time);

    fprintf(storage_file, "____Potentials:\n");
    one_line = "";
    for(int i=0; i<neurons_in_simulation; i++)
        one_line += VFFile::convertDoubleToString(neuron_array[i]->V) + " ";
    fprintf(storage_file, one_line.c_str());
    fprintf(storage_file, "\n");
}

int SimulationSingleton::outputCurrentsInFile(double time){
    fprintf(current_file, "_________Time = %.2f\n", time);

    one_line = "";
    for(int i=0; i<neurons_in_simulation; i++)
        one_line += VFFile::convertDoubleToString(neuron_array[i]->I) + " ";
    fprintf(current_file, one_line.c_str());
    fprintf(current_file, "\n");
}

void SimulationSingleton::outputSpikesInFile(){
    one_line = "";
    for(int i=0; i<spikes_resent; i++){
        one_line += VFFile::convertDoubleToString(spikes_time[i]) + " " + \
                VFFile::convertDoubleToString(spikes_numbers[i]) + "\n";
    }
    fprintf(spike_file, one_line.c_str());

    spikes_resent = 0;
}

void SimulationSingleton::outputSynapseDataInFile(double time){
    buffer = synapse_array[0]->getInnerData();
    if(buffer){
        fprintf(synapse_file, "_________Time = %.2f\n", time);
        one_line = "";
        for(int i=1; i<buffer[0]; i++){
            one_line += VFFile::convertDoubleToString(buffer[i]) + " ";
        }
        fprintf(synapse_file, one_line.c_str());
        fprintf(synapse_file, "\n");
    }
}

int SimulationSingleton::outputWeightsInFile(double time){
    fprintf(weight_file, "_________Time = %.2f\n", time);
    for(int i=0; i < N; i++){
        one_line = "";
        for(int j=0; j < N; j++){
            buf0 = i * N + j;
            if(connectivity_matrix[buf0] > -1)
                buf1 = synapse_array[buf0]->weight;
            else
                buf1 = 0;
            one_line+= VFFile::convertDoubleToString(buf1) + " ";
        }
        fprintf(weight_file, one_line.c_str());
        fprintf(weight_file, "\n");
    }

    return 0;
}

int SimulationSingleton::exportNeurons(){
    if(!(neuron_array[0]->exportData()))
        return 1;

    FILE* fid = fopen(FILE_EXP_NEURONS, "w");
    double *arr;
    for(int i=0; i<N; i++){
        one_line = VFFile::convertDoubleToString(i) + ": ";
        arr = neuron_array[i]->exportData();
        for(int j=1;j<arr[0];j++)
            one_line += VFFile::convertDoubleToString(arr[j]) + " ";
        one_line += VFFile::convertDoubleToString(Inoise[i]) + "\n";
        fprintf(fid, one_line.c_str());
    }
    fclose(fid);
    return 0;
}

int SimulationSingleton::exportSynapses(){
    if(!(synapse_array[0]->exportData()))
        return 1;

    FILE* fid = fopen(FILE_EXP_SYNAPSES, "w");
    double *arr;

    buf0 = 0;
    for(int i=0; i<N*N; i++){
        if(connectivity_matrix[i]>-1)
            buf0++;
    }
    one_line = "Number of synapses = " + \
            VFFile::convertDoubleToString(buf0) + "\n";
    fprintf(fid, one_line.c_str());

    for(int i=0; i<N*N;i++)
        if(connectivity_matrix[i]>-1){
            one_line = VFFile::convertDoubleToString(i)+ ": ";
            arr = synapse_array[i]->exportData();
            one_line += "from " + VFFile::convertDoubleToString(arr[1]) + \
                    " to " + VFFile::convertDoubleToString(arr[2]) + "; ";
            for(int j=3; j<arr[0];j++)
                one_line += VFFile::convertDoubleToString(arr[j]) + " ";
            one_line += "\n";
            fprintf(fid, one_line.c_str());
        }
    fclose(fid);
    return 0;
}

int SimulationSingleton::importNeurons(){
    FILE* fid = fopen(FILE_IMP_NEURONS, "r");
    if(!fid){
        error_number = 7;
        fclose(fid);
        return error_number;
    }
    double* arr;

    buf0 = 0;
    buf01 = 0;
    buf00 = neuron_array[0]->numEssentialVariables();
    arr = new double[buf00];
    for(int i=0; i<N; i++){
        fscanf(fid, "%d: ", &buf0);
        for(int j=0; j<buf00; j++){
            fscanf(fid, "%f", &buf01);
            arr[j] = buf01;
        }
        neuron_array[i]->importData(arr);

        fscanf(fid, "%f", &buf01);
        Inoise[i] = buf01;
        Inoise2[i] = buf01;
    }

    fclose(fid);
    return 0;
}

int SimulationSingleton::importSynapses(){
    FILE* fid = fopen(FILE_IMP_SYNAPSES, "r");
    if(!fid){
        error_number = 7;
        fclose(fid);
        return error_number;
    }
    double* arr;

    for(int i=0; i<N*N;i++){
        connectivity_matrix[i] = -1;
    }

    buf0 = 0;
    buf00 = 0;
    buf01 = 0;
    buf000 = synapse_array[0]->numEssentialVariables();
    arr = new double[buf000];
    fscanf(fid, "Number of synapses = %d", &buf0);
    for(int i=0; i<buf0; i++){
        fscanf(fid, "%d: ", &buf00);
        fscanf(fid, "from %d ", &buf00);
        arr[0] = buf00;
        fscanf(fid, "to %d; ", &buf00);
        arr[1] = buf00;
        for(int j=2; j<buf000; j++){
            fscanf(fid, "%f ", &buf01);
            arr[j] = buf01;
        }
        synapse_array[(int)arr[0]*N+(int)arr[1]]->importData(arr);
        connectivity_matrix[(int)arr[0]*N+(int)arr[1]] = \
                (int)arr[0]*N+(int)arr[1];
        synapse_array[(int)arr[0]*N+(int)arr[1]]->setDeliveries(dt);
    }

    fclose(fid);
    return 0;
}
