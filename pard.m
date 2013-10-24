%% Print Analysed Raw Data
function pard
global pi
pi = 3.14159265359;
    fprintf('\n\t\tNew instance %4.0f\n', random('Uniform', 1, 99999));

    LOAD_NEURON_PARAMETERS = 0;
    LOAD_POTENTIALS = 0;
    LOAD_SYNAPSE_DATA = 0;
    LOAD_SPIKES = 1;
    LOAD_SYNAPSES = 0;
    LOAD_CURRENTS = 0;
    LOAD_STIM_SYNAPSE = 0;
    LOAD_STIM_NEURON = 0;
    
    CALCULATE_AVERAGE_POTENTIAL = LOAD_POTENTIALS * 1;              %Plots an average potential vs time
    CALCULATE_PSD_POTENTIAL = CALCULATE_AVERAGE_POTENTIAL * 0;      %Plots a PSD (power spectral density) based on neuron potentials
    CALCULATE_AVERAGE_CURRENT = LOAD_CURRENTS * 1;                  %Plots average synaptic current (with SD) vs time
    CALCULATE_PSD_CURRENT = CALCULATE_AVERAGE_CURRENT * 1;          %Plots a PSD (power spectral density) based on synaptic currents
    CALCULATE_BG_CURRENT_HISTOGRAM = LOAD_NEURON_PARAMETERS * 1;    %Plots a histogram of background currents of neurons
    CALCULATE_STIM_NEURON_HISTOGRAM = LOAD_STIM_NEURON * 1;         %Plots a histogram of background currents of neurons
    CALCULATE_STIM_SYNAPSE_HISTOGRAM = LOAD_STIM_SYNAPSE * 1;         %Plots a histogram of background currents of neurons
    
    CALCULATE_ACTIVITY_HISTOGRAM = LOAD_SPIKES * 0;                 %Plots an activity histogram (activity of network vs time)
    CALCULATE_AVERAGE_ISI = LOAD_SPIKES * 0;                        %Plots average InterSpike Interval vs time
    CALCULATE_ISI_HISTOGRAM = LOAD_SPIKES * 0;                      %Plots a histogram of ISI
    
    FIND_BURSTS = CALCULATE_ACTIVITY_HISTOGRAM * 0;                 %Finds bursts and marks them on other plots
    BURSTS_STATISTICS = FIND_BURSTS * 1;                            %Calculates and outputs some statistics about bursts, like length, amount os spikes, etc.
    CALCULATE_IBI = FIND_BURSTS * 0;    PLOT_IBI = CALCULATE_IBI * 0;%Calculates (and Plots resp.) time since last burst vs time.
    CALCULATE_IBI_HISTOGRAM = CALCULATE_IBI * 1;                    %Plots InterBurst Interval histogram vs time.
    FULL_BURST_AMPLITUDE = FIND_BURSTS * LOAD_SYNAPSES * 1;         %Calculates amplitude of bursts in terms of active synapses. Very slow.

%%%%%%%%%%%%%%%%%%%%%%%%%%%LOADING FROM FILES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Analyzing parameters.txt
    fid = fopen('data/parameters.txt');
    fscanf(fid,'________Simulation parameters:\n');
    ton = fscanf(fid,'Type of neurons = %s');
    tos = fscanf(fid,'\nType of synapses = %s');
    N = fscanf(fid,'Neurons in simulation = %f;');
    p = fscanf(fid,'\nProbability of connection = %f;');
    time_length = fscanf(fid,'\nLength of simulation (msec) = %f;');
    tbe = fscanf(fid,'\nTime between exports (msec) = %f');
    fscanf(fid,'%s', 1);
    tbvie = fscanf(fid,'\nTime between I/V exports (msec) = %f');
    dt = fscanf(fid,';\nTime-step (msec) = %f;');
    fscanf(fid,'\n');
    tbwe = fscanf(fid,'\nTime between weight exports = %f;');
    aoin = fscanf(fid,'\nAmount of inhibitory neurons = %f;');
    tostim = fscanf(fid,'\nType of stimulation = %f;');
    syn_noise_freq = fscanf(fid,'\nSynaptic noise frequency = %f;');
    tau_stim = fscanf(fid,'\nTau stimulation = %f;');
    min_noise = fscanf(fid,'\nMin noise = %f;');
    max_noise = fscanf(fid,'ax noise = %f;');
    mean_stim = fscanf(fid,'\nMean stimulation = %f;');
    sd_stim = fscanf(fid,'\nSigma stimulation = %f;');
    import_neurons = fscanf(fid,'\nImport neurons = %d;');
    import_synapses = fscanf(fid,'synapses = %d;');
    fclose(fid);

    number_of_exports = floor(time_length / tbvie) + 1;
    fprintf('Parameters loaded.\n');

    %% Analyzing output.txt%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOGLER = 1;
    if LOAD_POTENTIALS == 1
        if TOGLER == 0
            fid = fopen('data/output.txt');
        else
            fid = fopen('data/output2.txt');
        end
        V = zeros(number_of_exports, N);
        time13 = zeros(1, number_of_exports);
        for i=1 : 1 : number_of_exports
            time13(i) = fscanf(fid, '_________Time = %f');
            fscanf(fid, '\n____Potentials:\n');
            asd = fscanf(fid, '%f ', [1 N]);
            for j=1 : 1 : N
                V(i, j) = asd(j);
            end
            fscanf(fid, '\n');
        end
        fclose(fid);
        fprintf('Output loaded.\n');
    end

    %% Analyzing spikes.txt
    if LOAD_SPIKES == 1
        fid = fopen('data/spikes.txt');
        rastr1 = fscanf(fid, '%f %f', [2 inf]);
        fclose(fid);
        rastr = rastr1(2,:);
        rastr_time = rastr1(1, :);
        rastr = rastr.';
        rastr_time = rastr_time.';
        num = size(rastr, 1);
        max_rastr_time = max(rastr_time);
        fprintf('Spikes loaded.\n');
        disp(size(rastr));
    end

    %% Analysing current.txt
    if LOAD_CURRENTS == 1
        fid = fopen('data/current.txt');
        I = zeros(number_of_exports, N);
        I_time = zeros(1, number_of_exports);
        for i=1 : 1 : number_of_exports
            I_time(i) = fscanf(fid, '_________Time = %f');
            asd = fscanf(fid, '%f ', [1 N]);
            for j=1 : 1 : N
                I(i, j) = asd(j);
            end
            fscanf(fid, '\n');
        end
        fclose(fid);
        fprintf('Currents loaded.\n');
    end
    
    fprintf('\tSimulation loaded.\n\n');
%%%%%%%%%%%%%%%%%%%%%%%%CALCULATING DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate average potential
    if CALCULATE_AVERAGE_POTENTIAL == 1
        V_average = zeros(1, number_of_exports);
        for i=1 : 1 : number_of_exports
            for j=1 : 1 : N
                V_average(i) =  V(i, j);
            end
            V_average(i) = V_average(i) ;
        end
        fprintf('Average potential calculated.\n');
    end

    %% Activity histogram
    if CALCULATE_ACTIVITY_HISTOGRAM == 1
        time_int_spikes = 1;
        integrated_spikes_amount = floor(time_length / time_int_spikes)+1;
        integrated_spikes = zeros(1, integrated_spikes_amount);
        integrated_spikes_time = zeros(1, integrated_spikes_amount);
        for i=1 : 1 : num
            integrated_spikes(floor(rastr_time(i) / time_int_spikes)+1) = ...
              integrated_spikes(floor(rastr_time(i) / time_int_spikes)+1) + 1;
        end
        integrated_spikes_time = 0 : time_int_spikes : ...
            (integrated_spikes_amount - 1) * time_int_spikes;
        integrated_spikes = integrated_spikes / N;
        fprintf('Activity histogram calculated.\n');
    end

    fprintf('\t\tData analysed.\n');
%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotnum = 1;    %number of the free figure
    %% Average potential%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if CALCULATE_AVERAGE_POTENTIAL == 1
        figure(plotnum);
        hold on
        if TOGLER == 0
            plot(time13, V_average / 15);
        else
            plot(time13, V_average / 15, 'k');
        end
        axis([0 max(max(time13)) 0 1.2]);
        title('Average potential');
        xlabel('Time, ms');
        ylabel('Potential / Vth');
        plotnum = plotnum + 1;
    end

    %% RASTER
    if LOAD_SPIKES == 1
        figure(plotnum);
        hold on;
        plot(rastr_time, rastr, 'k.','MarkerSize',1);
        if FIND_BURSTS == 1
            for i=1 : 1 : amount_of_bursts
                plot([burst_time(i) burst_time(i)], [0 N], 'r');
            end
        end
        plot([max_rastr_time max_rastr_time], [0 N], 'g.-');
        title('Raster plot of spiking neurons');
        axis([0 time_length 0 N]);
        xlabel('Time, ms');
        ylabel('Number of neuron');
        hold off;
        plotnum = plotnum + 1;
    end

    %% Activity
    if CALCULATE_ACTIVITY_HISTOGRAM == 1
        figure(plotnum);
        hold on;
        plot(integrated_spikes_time, integrated_spikes);
        plot([max_rastr_time max_rastr_time], [0 N], 'g.-');
        if FIND_BURSTS == 1
            for i=1 : 1 : amount_of_bursts
                plot([burst_time(i) burst_time(i)], [0 burst_amps(i)/N], 'r');
            end
        end
        title('Fraction of neurons spiking per time');
        axis([0 time_length 0 1]);
        xlabel('Time, ms');
        ylabel('Fraction of network');
        hold off;
        plotnum = plotnum + 1;
    end
    

    fprintf('\n');
    
end

function answ = normalDistrTheor(average, sd, index)
global pi
    answ = 1 ./ (sqrt(2.*pi).*sd) .* exp(- (index - average).^2 ./ (2 .* sd.^2));
end

function [index answ] = createDefiniteHistogram(data, min, max, step, N)
    if step > 0
        index = min : step : max;
        N = length(index);
    else
        step = (max - min) / (N - 1);
        index = min : step : max;
    end
    answ = zeros(1, N);
    for i=1 : 1 : length(data)
        buf = floor( (data(i) - min) / step) + 1;
        answ(buf) = answ(buf) + 1;
    end
end