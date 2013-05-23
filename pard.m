%% Print Analysed Raw Data
function pard
global pi
pi = 3.14159265359;
    fprintf('\n\t\tNew instance %4.0f\n', random('Uniform', 1, 99999));

    LOAD_NEURON_PARAMETERS = 1;
    LOAD_POTENTIALS = 0;
    LOAD_SYNAPSE_DATA = 0;
    LOAD_SPIKES = 0;
    LOAD_SYNAPSES = 0;
    LOAD_CURRENTS = 0;
    
    CALCULATE_AVERAGE_POTENTIAL = LOAD_POTENTIALS * 1;              %Plots an average potential vs time
    CALCULATE_AVERAGE_CURRENT = LOAD_CURRENTS * 1;                  %Plots average synaptic current (with SD) vs time
    CALCULATE_BG_CURRENT_HISTOGRAM = LOAD_NEURON_PARAMETERS * 1;    %Plots a histogram of background currents of neurons
    
    CALCULATE_ACTIVITY_HISTOGRAM = LOAD_SPIKES * 1;                 %Plots an activity histogram (activity of network vs time)
    CALCULATE_AVERAGE_ISI = LOAD_SPIKES * 0;                        %Plots average InterSpike Interval vs time
    CALCULATE_ISI_HISTOGRAM = LOAD_SPIKES * 0;                      %Plots a histogram of ISI
    
    FIND_BURSTS = CALCULATE_ACTIVITY_HISTOGRAM * 1;                 %Finds bursts and marks them on other plots
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

    %% Analysing neuron_parameters.txt
    if LOAD_NEURON_PARAMETERS == 1
        fid = fopen('data/neuron_parameters.txt');
        fscanf(fid,'\n________Neurons potential threshold:\n');
        Vth = fscanf(fid,'%f ');
        Vth = Vth';
        fscanf(fid,'\n________Neurons Vrest:\n');
        Vrest = fscanf(fid,'%f ');
        Vrest = Vrest';
        fscanf(fid,'\n________Neurons Ibg:\n');
        Ibg = fscanf(fid,'%f ');
        Ibg = Ibg';
        fclose(fid);
        fprintf('Neuron parameters loaded.\n');
    end
  
    %% Analyzing output.txt
    if LOAD_POTENTIALS == 1
        fid = fopen('data/output.txt');
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

    %% Analyzing synapse.txt
    if LOAD_SYNAPSE_DATA == 1
        fid = fopen('data/synapse.txt');
        X = zeros(1, number_of_exports);
        Y = zeros(1, number_of_exports);
        Z = zeros(1, number_of_exports);
        U = zeros(1, number_of_exports);
        time2 = zeros(1, number_of_exports);
        for i=1 : 1 : number_of_exports
            time2(i) = fscanf(fid, '_________Time = %f\n');
            X(i) = fscanf(fid, '%f ', 1);
            Y(i) = fscanf(fid, '%f ', 1);
            Z(i) = fscanf(fid, '%f ', 1);
            U(i) = fscanf(fid, '%f ', 1);
            fscanf(fid, '\n');
        end
        fclose(fid);
        fprintf('Synapse loaded.\n');
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
    end

    %% Analyzing export/synapses.txt
    if LOAD_SYNAPSES == 1
        fid = fopen('export/synapses.txt', 'r');
        Nsyn = fscanf(fid, 'Number of synapses = %d');
        synnum = zeros(1, N);

        for i=1 : 1 : Nsyn
            buf = fscanf(fid, '%s', 1);
            while(~strcmp(buf,'from'))
                buf = fscanf(fid, '%s', 1);
            end
            buf00 = fscanf(fid, '%d', 1);
            synnum(buf00+1) = synnum(buf00+1) + 1;
        end
        fclose(fid);
        fprintf('Connections loaded.\n');
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
                V_average(i) =  V_average(i) + V(i, j);
            end
            V_average(i) = V_average(i) / N;
        end
        fprintf('Average potential calculated.\n');
    end
    
    %% Calculate average current
    if CALCULATE_AVERAGE_CURRENT == 1
        I_average = zeros(1, number_of_exports);
        I_sd = zeros(1, number_of_exports);
        for i=1 : 1 : number_of_exports
            for j=1 : 1 : N
                I_average(i) =  I_average(i) + I(i, j);
            end
            I_average(i) = I_average(i) / N;
            for j=1 : 1 : number_of_exports
                I_sd(i) = I_sd(i) + (I_average(i) - I(i, j))^2;
            end
            I_sd(i) = sqrt(I_sd(i) / (N - 1));
        end
        fprintf('Average current calculated.\n');
    end
    
    %% BG current histogram
    if CALCULATE_BG_CURRENT_HISTOGRAM == 1
        I_bg_min = min(min(Ibg));
        I_bg_max = max(max(Ibg));
        I_bg_N = 20;    %SET AMOUNT OF HISTOGRAM PARTS HERE
        I_bg_d = (I_bg_max - I_bg_min) / (I_bg_N - 1);
        I_bg_hist = zeros(1, I_bg_N);
        for i=1 : 1 : N
            buf00 = floor((Ibg(i) - I_bg_min) / I_bg_d) + 1;
            I_bg_hist(buf00) = I_bg_hist(buf00) + 1;
        end
        I_bg_ind = I_bg_min : I_bg_d : I_bg_max;
        I_bg_thrs = 15;%pA
        I_bg_N1 = 0;
        I_bg_N2 = 0;
        for i=1 : 1 : N
            if(Ibg(i)>I_bg_thrs)
                I_bg_N1 = I_bg_N1 + 1;
            else
                I_bg_N2 = I_bg_N2 + 1;
            end
        end
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

    %% Average InterSpike Interval vs time
    if CALCULATE_AVERAGE_ISI == 1
        ISI = 1./(integrated_spikes * N);

        %Overall ISI
        spikes_in_neuron = ones(1, N);
        for i=1 : 1 : num
            spikes_in_neuron(rastr(i)+1) = spikes_in_neuron(rastr(i)+1) + 1;
        end
        overall_ISI = max_rastr_time ./ spikes_in_neuron;
        min_ISI_hist = 0;
        max_ISI_hist = max(max(overall_ISI));
        N_ISI_hist = N / 5; %SET AMOUNT OF HISTOGRAM PARTS HERE
        dISI_hist = (max_ISI_hist - min_ISI_hist) / (N_ISI_hist - 1);
        average_hist_ISI = zeros(1, N_ISI_hist);
        for i=1 : 1 : N
            buf00 = floor(overall_ISI(i) / dISI_hist) + 1;
            average_hist_ISI(buf00) = average_hist_ISI(buf00) + 1;
        end
        average_hist_ISI_ind = min_ISI_hist : dISI_hist : max_ISI_hist;
        fprintf('Average InterSpike Interval calculated.\n');
    end
  
    %% True ISI calculation - ISI histogram
    if CALCULATE_ISI_HISTOGRAM == 1
        spike_per_neuron = zeros(N, 1);
        for i=1 : 1 : length(rastr)
            neuronum = rastr(i) + 1;
            neurotime = rastr_time(i);
            spike_per_neuron(neuronum, 1) = spike_per_neuron(neuronum, 1) + 1;
            spike_per_neuron(neuronum, spike_per_neuron(neuronum, 1)+1) = neurotime;
        end
        diff_spike_per_neuron = zeros(size(spike_per_neuron));
        min_true_ISI = 0;
        max_true_ISI = 0;
        for i=1 : 1 : N
            diff_spike_per_neuron(i, 1) = spike_per_neuron(i, 1);
            diff_spike_per_neuron(i, 2) = spike_per_neuron(i, 2);
            if diff_spike_per_neuron(i, 2) > max_true_ISI
                max_true_ISI = diff_spike_per_neuron(i, 2);
            end
            if spike_per_neuron(i,1) > 1
                for j=2 : 1 : spike_per_neuron(i, 1)
                    diff_spike_per_neuron(i, j+1) = ...
                        spike_per_neuron(i, j+1) - spike_per_neuron(i, j);
                    if diff_spike_per_neuron(i, j+1) > max_true_ISI
                        max_true_ISI = diff_spike_per_neuron(i, j+1);
                    end
                end
            end
        end
        N_true_ISI = max_true_ISI / 10; %SET AMOUNT OF HISTOGRAM PARTS HERE
        true_ISI = zeros(1, N_true_ISI+1);
        d_true_ISI = (max_true_ISI - min_true_ISI) / N_true_ISI;
        for i=1 : 1 : N
            if diff_spike_per_neuron(i,1) > 0
                for j=1 : 1 : diff_spike_per_neuron(i, 1)
                    buf00 = floor(diff_spike_per_neuron(i, j+1) / d_true_ISI) + 1;
                    true_ISI(buf00) = true_ISI(buf00) + 1;
                end
            end
        end
        true_ISI_ind = min_true_ISI : d_true_ISI : max_true_ISI;
        fprintf('InterSpike Interval histogram calculated.\n');
    end
  
    %% Finding busrsts
    if FIND_BURSTS == 1
        burst_threshold = 0.2;          % % of network
        burst_width_threshold = 0.02;   % % of network
        window_size = 100;              %ms
        find_bursts = zeros(1, integrated_spikes_amount);
        for i=1 : 1 : integrated_spikes_amount
            if integrated_spikes(i) > burst_threshold
                find_bursts(i) = integrated_spikes(i);
            end
        end
        i = 1;
        while i < integrated_spikes_amount
            j = i+1;
            while integrated_spikes_time(j) - integrated_spikes_time(i) ...
                    < window_size && j+1 <= integrated_spikes_amount
                if find_bursts(i) > find_bursts(j)
                    find_bursts(j) = 0;
                else
                    find_bursts(i) = 0;
                    break;
                end
                j = j+1;
            end
            i = j;
        end
        amount_of_bursts = 0;
        burst_amps = zeros(1);
        burst_time = zeros(1);
        burst_spikes = zeros(1);
        burst_start = zeros(1);
        burst_end = zeros(1);
        for i=1 : 1 : integrated_spikes_amount
            if find_bursts(i) > 0
                amount_of_bursts = amount_of_bursts + 1;
                burst_amps(amount_of_bursts) = find_bursts(i) * N;
                burst_time(amount_of_bursts) = integrated_spikes_time(i);

                j = i-1;
                burst_spikes(amount_of_bursts) = integrated_spikes(i);
                while integrated_spikes(j) > burst_width_threshold
                    burst_spikes(amount_of_bursts) = ...
                        burst_spikes(amount_of_bursts) + integrated_spikes(j);
                    j = j - 1;
                    if j < 1
                        j = 1;
                        break;
                    end
                end
                burst_start(amount_of_bursts) = integrated_spikes_time(j);

                j = i+1;
                while integrated_spikes(j) > burst_width_threshold
                    burst_spikes(amount_of_bursts) = ...
                        burst_spikes(amount_of_bursts) + integrated_spikes(j);
                    j = j + 1;
                    if j > integrated_spikes_amount
                        j = integrated_spikes_amount;
                        break;
                    end
                end
                burst_end(amount_of_bursts) = integrated_spikes_time(j);
            end
        end
        burst_spikes = burst_spikes * N;
        burst_length = burst_end - burst_start;
        fprintf('Bursts found.\n');
    end
  
    %% Bursts statistics: average size, length, etc.
    if BURSTS_STATISTICS == 1
        % Average
        burst_length_av = 0;
        burst_activity_amp_av = 0;
        burst_spikes_amp_av = 0;
        for i=1 : 1 : amount_of_bursts
            burst_spikes_amp_av = burst_spikes_amp_av + burst_spikes(i);
            burst_length_av = burst_length_av + burst_length(i);
            burst_activity_amp_av = burst_activity_amp_av + burst_amps(i);
        end
        burst_spikes_amp_av = burst_spikes_amp_av / amount_of_bursts;
        burst_length_av = burst_length_av / amount_of_bursts;
        burst_activity_amp_av = burst_activity_amp_av / amount_of_bursts;
        
        % SD
        burst_length_sd = 0;
        burst_activity_amp_sd = 0;
        burst_spikes_amp_sd = 0;
        for i=1 : 1 : amount_of_bursts
            burst_length_sd = burst_length_sd + (burst_length(i) - burst_length_av)^2;
            burst_activity_amp_sd = burst_activity_amp_sd + (burst_activity_amp_av - burst_amps(i))^2;
            burst_spikes_amp_sd = burst_spikes_amp_sd + (burst_spikes_amp_av - burst_spikes(i))^2;
        end
        burst_spikes_amp_sd = sqrt(burst_spikes_amp_sd / (amount_of_bursts - 1));
        burst_length_sd = sqrt(burst_length_sd / (amount_of_bursts - 1));
        burst_activity_amp_sd = sqrt(burst_activity_amp_sd / (amount_of_bursts - 1));

        fprintf('Normal burst statistics calculated.\n');
    end
  
    %% InterBurst Interval calculation
    if CALCULATE_IBI == 1
        IBI = zeros(size(burst_amps));
        IBI(1) = burst_time(1);
        for i=2 : 1 : amount_of_bursts
            IBI(i) = burst_time(i) - burst_time(i-1);
        end
        fprintf('InterBurst Interval calculated.\n');
    end
    
    %% InterBurst Interval histogram
    if CALCULATE_IBI_HISTOGRAM == 1
        IBI_min = 0;
        IBI_max = max(max(IBI));
        IBI_N = 1000;   %SET AMOUNT OF HISTOGRAM PARTS HERE
        IBI_d = (IBI_max - IBI_min) / (IBI_N -1);
        IBI_hist = zeros(1, IBI_N);
        for i=1 : 1 : amount_of_bursts
            IBI_hist(floor(IBI(i) / IBI_d) + 1) = ...
                IBI_hist(floor(IBI(i) / IBI_d) + 1) + 1;
        end
        IBI_hist_ind = IBI_min : IBI_d : IBI_max;
        fprintf('InterBurst Interval histogram calculated.\n');
    end

    %% Burst amplitudes as spikes in synapses
    if FULL_BURST_AMPLITUDE == 1
        burst_full_amps = zeros(size(burst_time));
        for i=1 : 1 : length(rastr)
            for j=1 : 1 : length(burst_time)
                if burst_time(j) - rastr_time(i) > -time_int_spikes && ...
                      burst_time(j) - rastr_time(i) <= 0
                    burst_full_amps(j) = burst_full_amps(j) + synnum(rastr(i)+1);
                end
            end
        end
        average_full_burst_amp = 0;
        for i=1 : 1 : length(burst_time)
            average_full_burst_amp = average_full_burst_amp + burst_full_amps(i);
        end
        average_full_burst_amp = average_full_burst_amp / length(burst_time);
        sd_full_burst_amp = 0;
        for i=1 : 1 : amount_of_bursts
            sd_full_burst_amp = sd_full_burst_amp + (average_full_burst_amp - burst_full_amps(i))^2;
        end
        sd_full_burst_amp = sqrt(sd_full_burst_amp/length(burst_time));
        fprintf('Synapse bursts amplitudes calculated.\n');
    end

    fprintf('\t\tData analysed.\n');
%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotnum = 1;    %number of the free figure
    %% Average potential
    if CALCULATE_AVERAGE_POTENTIAL == 1
        figure(plotnum);
        plot(time13, V_average / 15);
        axis([0 max(max(time13)) 0 1]);
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
    
    %% Synapse data
    if LOAD_SYNAPSE_DATA == 1
        figure(plotnum);
        hold on;
        title('Average x, y, z, u.');
        plot(time2, X, 'b');
        plot(time2, Y, 'k');
        plot(time2, Z, 'y');
        plot(time2, U, 'r');
        axis([0 time_length 0 1]);
        legend('X', 'Y', 'Z', 'U');
        xlabel('Time, ms');
        ylabel('Fraction');
        hold off;
        plotnum = plotnum + 1;
    end

    %% IBI
    if PLOT_IBI == 1
        figure(plotnum);
        hold on;
        plot(burst_time, IBI, 'k');
        title('InterBurst Interval');
        xlabel('Time, ms');
        ylabel('InterBurst Interval, ms');
        hold off;
        plotnum = plotnum + 1;
    end

    %% Average ISI
    if CALCULATE_AVERAGE_ISI == 1
        figure(plotnum);
        hold on;
        plot(integrated_spikes_time, ISI, 'k');
        title('Average InterSpikes Interval');
        axis([0 time_length 0 1]);
        xlabel('Time, ms');
        ylabel('InterSpikes Interval, ms');
        hold off;
        plotnum = plotnum + 1;
    end

    %% Average ISI histogram
    if CALCULATE_AVERAGE_ISI == 1
        figure(plotnum);
        hold on;
        plot(average_hist_ISI_ind, average_hist_ISI, 'k');
        title('InterSpikes Interval Histogram');
        xlabel('Time, ms');
        ylabel(' ');
        hold off;
        plotnum = plotnum + 1;
    end

    %% ISI histogram
    if CALCULATE_ISI_HISTOGRAM == 1
        figure(plotnum);
        hold on;
        plot(true_ISI_ind, true_ISI, 'k*');
        title('InterSpikes Interval Histogram')
        xlabel('Time, ms');
        ylabel(' ');
        hold off;
        plotnum = plotnum + 1;
    end

    %% IBI histogram
    if CALCULATE_IBI_HISTOGRAM == 1
        figure(plotnum);
        hold on;
        plot(IBI_hist_ind, IBI_hist, 'k*');
        title('InterBurst Interval Histogram')
        xlabel('Time, ms');
        ylabel(' ');
        hold off;
        plotnum = plotnum + 1;
    end

    %% Currents
    if LOAD_CURRENTS == 1
        figure(plotnum);
        hold on;
        errorbar(I_time, I_average, I_sd, 'k');
        title('Average currents with sd')
        xlabel('Time, ms');
        ylabel('Synaptic current, pA');
        axis([0 time_length min(min(I))*1.1 max(max(I))*1.1]);
        hold off;
        plotnum = plotnum + 1;
    end
    
    %% Background currents
    if CALCULATE_BG_CURRENT_HISTOGRAM == 1
        figure(plotnum);
        hold on;
        plot(I_bg_ind, I_bg_hist / N, 'k*');
        arr = min_noise : (max_noise - min_noise) / 1000 : max_noise;
        plot(arr, normalDistrTheor(mean_stim, sd_stim, arr), 'r');%.*max(max(I_bg_hist))*10
        plot([I_bg_thrs I_bg_thrs], [0 normalDistrTheor(mean_stim, sd_stim, I_bg_thrs)],'g');
        title('Background current histogram')
        xlabel('I_{bg}, pA');
        ylabel(' ');
        legend('Experimental distribution', 'Theoretical', 'I_crit');
        hold off;
        plotnum = plotnum + 1;
    end

%%%%%%%%%%%%%%%%%%%%%%%%PRINTING DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Burst data:
    if FIND_BURSTS == 1
        fprintf('Burst info: Total %d bursts:', amount_of_bursts);
        for i=1 : 1 : amount_of_bursts
            fprintf(' %d', burst_amps(i));
        end
        fprintf('\n');
    end
    if BURSTS_STATISTICS == 1
        fprintf('\nBursts'' neuron amps: Mean: %f, sd: %f.', burst_activity_amp_av, burst_activity_amp_sd);
        fprintf('\nBursts'' lengths: Mean: %f, sd: %f.', burst_length_av, burst_length_sd);
        fprintf('\nBursts'' spikes amounts: Mean: %f, sd: %f.', burst_spikes_amp_av, burst_spikes_amp_sd);
    end
    if FULL_BURST_AMPLITUDE == 1
        fprintf('\nBurst unnormal info. Mean: %f, sd: %f.', average_full_burst_amp, sd_full_burst_amp);
    end
    
    %% BG current data:
    if CALCULATE_BG_CURRENT_HISTOGRAM == 1
        fprintf('\n\nBackground current info. N1 (active): %d (%.1f%%), N2(not-active): %d (%.1f%%).', I_bg_N1, I_bg_N1 / N*100, I_bg_N2, I_bg_N2/N*100);
    end
%     fid = fopen('current_hist_data.txt', 'w');
%     fprintf(fid, 'I_bg | Amount of neurons with this Ibg.\n\n');
%     for i=1 : 1 : I_bg_N
%         fprintf(fid, '%.5f %d\n', I_bg_ind(i), I_bg_hist(i));
%     end
%     fclose(fid);

    fprintf('\n');
    
end

function answ = normalDistrTheor(average, sd, index)
global pi
    answ = 1 ./ (sqrt(2.*pi).*sd) .* exp(- (index - average).^2 ./ (2 .* sd.^2));
end
