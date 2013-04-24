#include "simulation.h"

int SimulationSingleton::loadParameters(){
    loadParametersFromFile();
    return errorReport();
}

int SimulationSingleton::loadParametersFromFile(){
    cout<<"\nLoading main parameters...\r";
    using namespace vf_file;

    buf31 = loadFileToString(FILE_IMP_PARAMETERS, DATAFILES);

    neurons_in_simulation = getParameterIni("Amount_of_neurons", buf31);
    N = neurons_in_simulation;
    amount_of_inh_neurons = getParameterIni(\
        "Amount_of_inhibitory_neurons", buf31);
    length_of_simulation = getParameterIni("Length_of_simulation", \
        buf31);
    dt = getParameterIni("dt", buf31);
    time_between_exports = getParameterIni("Time_between_exports", \
        buf31);
    time_between_vi_exports = getParameterIni(\
        "Time_between_vi_exports", buf31);
    type_of_neuron = getParameterIni("Type_of_neuron", buf31);
    type_of_synapse = getParameterIni("Type_of_synapse", buf31);
    type_of_stimulation = getParameterIni("Type_of_stimulation", buf31);
    syn_noise_freq_mean = getParameterIni(\
        "Synaptic_noise_mean_frequency", buf31);
    tau_stim = getParameterIni("Tau_stimulation", buf31);
    stim_start = getParameterIni("Start_of_stimulation", buf31);
    time_between_weight_exports = getParameterIni(\
        "Time_between_weight_exports", buf31);
    Imin = getParameterIni("Min_noise", buf31);
    Imax = getParameterIni("Max_noise", buf31);
    Imean = getParameterIni("Mean_stimulation", buf31);
    Isd = getParameterIni("Sigma_stimulation", buf31);
    import_neurons = getParameterIni("Import_network_neurons", buf31);
    import_synapses = getParameterIni("Import_network_synapses", buf31);


    buf31 = loadFileToString(FILE_SYS_DATA, DATAFILES);

    number_of_threads = getParameterIni("Number_of_threads", buf31);
    enable_test = getParameterIni("Test", buf31);

    cout<<"Main parameters loaded.   ";
    return 0;
}

void SimulationSingleton::outputChangingData(){
    using namespace vf_discrete;
    if(inBetween(time_now, time_between_vi_exports, dt)){
        outputChangingDataInFile(time_now);
        outputSynapseDataInFile(time_now);
        outputCurrentsInFile(time_now);
    }
    if(inBetween(time_now, time_between_exports, dt)){
        outputSpikesInFile();
    }
    if(inBetween(time_now, time_between_weight_exports, dt)){
        outputWeightsInFile(time_now);
    }
}

int SimulationSingleton::setOutputFile(){
    cout<<"\nOpening output files...\r";
    using namespace vf_file;
    if(!enable_test)
        if(!(tryFile(getFilenameFromIni(DATAFILES, FILE_EXP_MAIN_OUTPUT))) || \
           !(tryFile(getFilenameFromIni(DATAFILES, FILE_EXP_SPIKES))) || \
           !(tryFile(getFilenameFromIni(DATAFILES, FILE_EXP_SYNAPSE))) || \
           !(tryFile(getFilenameFromIni(DATAFILES, FILE_EXP_WEIGHTS))) || \
           !(tryFile(getFilenameFromIni(DATAFILES, FILE_EXP_CURRENT))) ){
            error_number = 2;
            return error_number;
        }

    storage_file = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_MAIN_OUTPUT).c_str(), "w");
    spike_file = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_SPIKES).c_str(), "w");
    synapse_file = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_SYNAPSE).c_str(), "w");
    weight_file = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_WEIGHTS).c_str(), "w");
    current_file = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_CURRENT).c_str(), "w");

    if(!storage_file || !synapse_file || !spike_file || !current_file){
        error_number = 2;
        return 2;
    }

    cout<<"Output files opened.   ";
    return errorReport();
}

void SimulationSingleton::closeOutputFile(){
    fclose(storage_file);
    fclose(spike_file);
    fclose(synapse_file);
    fclose(weight_file);
    fclose(current_file);
}

int SimulationSingleton::outputParametersInFile(){
    using namespace vf_file;
    if(!enable_test)
        if(!(tryFile(getFilenameFromIni(DATAFILES, FILE_EXP_PARAMETERS))) || \
           !(tryFile(getFilenameFromIni(DATAFILES, FILE_EXP_NEURON_PARAMS))) ){
            error_number = 2;
            return error_number;
        }
    FILE* param_file = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_PARAMETERS).c_str(), "w");
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
    fprintf(param_file, "Length of simulation (msec) _=_ %.2f;\n", \
                        length_of_simulation);
    fprintf(param_file, "Time between exports (msec) _=_ %.2f;\n", \
                        time_between_exports);
    fprintf(param_file, "Time between I/V exports (msec) _=_ %.2f;\n", \
                        time_between_vi_exports);
    fprintf(param_file, "Time-step (msec) _=_ %.4f;\n", \
                        dt);
    fprintf(param_file, "\n");

    fclose(param_file);

    param_file = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_NEURON_PARAMS).c_str(), "w");
    if(!param_file){
        error_number = 2;
        return 2;
    }

    fprintf(param_file, "________Neurons potential threshold:\n");
    buf31 = "";
    for(int i=0; i<N; i++)
        buf31+= convertDoubleToString(neuron_array[i]->Vth) + " ";
    fprintf(param_file, buf31.c_str());
    fprintf(param_file, "\n");

    fprintf(param_file, "________Neurons Vrest:\n");
    buf31 = "";
    for(int i=0; i<N; i++)
        buf31+= convertDoubleToString(neuron_array[i]->Vrest) + " ";
    fprintf(param_file, buf31.c_str());
    fprintf(param_file, "\n");

    fprintf(param_file, "________Neurons Ibg:\n");
    buf31 = "";
    for(int i=0; i<N; i++)
        buf31+= convertDoubleToString(Inoise[i]) + " ";
    fprintf(param_file, buf31.c_str());
    fprintf(param_file, "\n");

    buf0 = 0;
    fprintf(param_file, "________Synapse delays:\n");
    for(int i=0; i<N; i++){
        buf31 = "";
        for(int j=0; j<N; j++){
            buf0 = i * N + j;
            if(connectivity_matrix[buf0]>-1)
                buf1 = synapse_array[connectivity_matrix[buf0]]->delay;
            else
                buf1 = 0;
            buf1 = buf1 - (buf1 * 100 - (int)(buf1 * 100))/100;
            buf31 += convertDoubleToString(buf1) + " ";
        }
        fprintf(param_file, buf31.c_str());
        fprintf(param_file, "\n");
    }
    fclose(param_file);

    return errorReport();
}

int SimulationSingleton::outputConnectivityMatrixInFile(){
    // Saves parameters in parameters.txt
    using namespace vf_file;
    if(!enable_test)
        if(!(tryFile(getFilenameFromIni(DATAFILES, FILE_EXP_CONN_MATR)))){
            error_number = 2;
            return error_number;
        }
    FILE *param_file = fopen (getFilenameFromIni(DATAFILES, FILE_EXP_CONN_MATR).c_str(), "w");
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
    return errorReport();
}

void SimulationSingleton::outputChangingDataInFile(double time){
    fprintf(storage_file, "_________Time = %.2f\n", time);

    fprintf(storage_file, "____Potentials:\n");
    one_line = "";
    for(int i=0; i<neurons_in_simulation; i++)
        one_line += vf_file::convertDoubleToString(neuron_array[i]->V) + " ";
    fprintf(storage_file, one_line.c_str());
    fprintf(storage_file, "\n");
}

int SimulationSingleton::outputCurrentsInFile(double time){
    fprintf(current_file, "_________Time = %.2f\n", time);

    one_line = "";
    for(int i=0; i<neurons_in_simulation; i++)
        one_line += vf_file::convertDoubleToString(neuron_array[i]->I) + " ";
    fprintf(current_file, one_line.c_str());
    fprintf(current_file, "\n");
    return 0;
}

void SimulationSingleton::outputSpikesInFile(){
    one_line = "";
    for(int i=0; i<spikes_resent; i++){
        one_line += vf_file::convertDoubleToString(spikes_time[i]) + " " + \
                vf_file::convertDoubleToString(spikes_numbers[i]) + "\n";
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
            one_line += vf_file::convertDoubleToString(buffer[i]) + " ";
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
                buf1 = synapse_array[connectivity_matrix[buf0]]->weight;
            else
                buf1 = 0;
            one_line+= vf_file::convertDoubleToString(buf1) + " ";
        }
        fprintf(weight_file, one_line.c_str());
        fprintf(weight_file, "\n");
    }

    return 0;
}

int SimulationSingleton::exportNeurons(){
    cout<<"\nEporting neurons...\r";
    using namespace vf_file;
    if(!(neuron_array[0]->exportData()))
        return 1;

    FILE* fid = fopen(getFilenameFromIni(DATAFILES, FILE_EXP_NEURONS).c_str(), "w");
    double *arr;
    for(int i=0; i<N; i++){
        one_line = convertDoubleToString(i) + ": ";
        arr = neuron_array[i]->exportData();
        for(int j=1;j<arr[0];j++)
            one_line += convertDoubleToString(arr[j]) + " ";
        one_line += convertDoubleToString(Inoise[i]) + "\n";
        fprintf(fid, one_line.c_str());
    }
    fclose(fid);
    cout<<"Neurons exported.  ";
    return 0;
}

int SimulationSingleton::exportSynapses(){
    cout<<"\nExport synapses...\r";
    using namespace vf_file;
    if(!(synapse_array[0]->exportData()))
        return 1;

    FILE* fid = fopen(getFilenameFromIni(DATAFILES, FILE_EXP_SYNAPSES).c_str(), "w");
    double *arr;

    one_line = "Number of synapses = " + \
            convertDoubleToString(M) + "\n";
    fprintf(fid, one_line.c_str());

    for(int i=0; i<N*N;i++)
        if(connectivity_matrix[i]>-1){
            one_line = convertDoubleToString(connectivity_matrix[i])+ ": ";
            arr = synapse_array[connectivity_matrix[i]]->exportData();
            one_line += "from " + convertDoubleToString(arr[1]) + \
                    " to " + convertDoubleToString(arr[2]) + "; ";
            for(int j=3; j<arr[0];j++)
                one_line += convertDoubleToString(arr[j]) + " ";
            one_line += "\n";
            fprintf(fid, one_line.c_str());
        }
    fclose(fid);
    cout<<"Synapses exported.   ";
    return 0;
}

int SimulationSingleton::importNeurons(){
    FILE* fid = fopen(vf_file::getFilenameFromIni(DATAFILES, FILE_IMP_NEURONS).c_str(), "r");
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
    FILE* fid = fopen(vf_file::getFilenameFromIni(DATAFILES, FILE_IMP_SYNAPSES).c_str(), "r");
    if(!fid){
        error_number = 7;
        fclose(fid);
        return error_number;
    }
    double* arr;

    for(int i=0; i<N*N;i++){
        connectivity_matrix[i] = -1;
    }

    M = 0;
    buf0 = 0;
    buf00 = 0;
    buf01 = 0;
    fscanf(fid, "Number of synapses = %d", &M);
    createNamedSynapses();
    buf000 = synapse_array[0]->numEssentialVariables();
    arr = new double[buf000];
    for(int i=0; i<M; i++){
        fscanf(fid, "%d: ", &buf00);
        fscanf(fid, "from %d ", &buf00);
        arr[0] = buf00;
        fscanf(fid, "to %d; ", &buf00);
        arr[1] = buf00;
        for(int j=2; j<buf000; j++){
            fscanf(fid, "%f ", &buf01);
            arr[j] = buf01;
        }
        connectivity_matrix[(int)arr[0]*N+(int)arr[1]] = buf0;
        synapse_array[buf0]->importData(arr);
        synapse_array[buf0]->setDeliveries(dt);
        buf0++;
    }
    fclose(fid);
    return 0;
}
