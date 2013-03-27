%Print Analysed Raw Data
function pard
  fprintf('\n\t\tNew instance %4.0f\n', random('Uniform', 1, 99999));
  
  %Analyzing parameters.txt
  fid = fopen('data/parameters.txt');
  fscanf(fid,'________Simulation parameters:\n');
  ton = fscanf(fid,'Type of neurons _=_ %s');
  tos = fscanf(fid,'\nType of synapses _=_ %s');
  N = fscanf(fid,'Neurons in simulation _=_ %f;');
  p = fscanf(fid,'\nProbability of connection _=_ %f;');
  time_length = fscanf(fid,'\nLength of simulation (msec) _=_ %f;');
  tbe = fscanf(fid,'\nTime between exports (msec) _=_ %f');
  dt = fscanf(fid,';\nTime-step (msec) _=_ %f;');
  fclose(fid);

  number_of_exports = floor(time_length / dt)+1;
  time = zeros(1, number_of_exports);
  for i=1 : 1 : number_of_exports
      time(i) = (i-1) * dt;
  end
  fprintf('Parameters loaded.\n');

%   %Analysing neuron_parameters.txt
%   fid = fopen('data/parameters.txt');
%   fscanf(fid,'\n________Neurons potential threshold:\n');
%   Vth = fscanf(fid,'%f ');
%   Vth = Vth';
%   fscanf(fid,'\n________Neurons Vrest:\n');
%   Vrest = fscanf(fid,'%f ');
%   Vrest = Vrest';
%   fclose(fid);
%   fprintf('Neuron parameters loaded.\n');
  
%   %Analyzing output.txt
%   fid = fopen('data/output.txt');
%   V = zeros(number_of_exports, N);
%   for i=1 : 1 : number_of_exports
%       time(i) = fscanf(fid, '_________Time = %f');
%       fscanf(fid, '\n____Potentials:\n');
%       asd = fscanf(fid, '%f ', [1 inf]);
%       for j=1 : 1 : N
%           V(i, j) = asd(j);
%       end
%       fscanf(fid, '\n');
%   end
%   fclose(fid);
%   fprintf('Output loaded.\n');

%   %Analyzing synapse.txt
%   fid = fopen('data/synapse.txt');
%   X = zeros(1, number_of_exports);
%   Y = zeros(1, number_of_exports);
%   Z = zeros(1, number_of_exports);
%   U = zeros(1, number_of_exports);
%   time2 = zeros(1, number_of_exports);
%   for i=1 : 1 : number_of_exports
%       time2(i) = fscanf(fid, '_________Time = %f\n');
%       X(i) = fscanf(fid, '%f ', 1);
%       Y(i) = fscanf(fid, '%f ', 1);
%       Z(i) = fscanf(fid, '%f ', 1);
%       U(i) = fscanf(fid, '%f ', 1);
%       fscanf(fid, '\n');
%   end
%   fclose(fid);
%   fprintf('Synapse loaded.\n');
  
  %Analyzing spikes.txt
  fid = fopen('data/spikes.txt');
  rastr1 = fscanf(fid, '%f %f', [2 inf]);
  fclose(fid);
  rastr = rastr1(2,:);
  rastr_time = rastr1(1, :);
  rastr = rastr.';
  rastr_time = rastr_time.';
  num = size(rastr, 1);
  fprintf('Spikes loaded.\n');
  max_rastr_time = max(rastr_time);
  
%   %Analyzing export/synapses.txt
%   fid = fopen('export/synapses.txt', 'r');
%   Nsyn = fscanf(fid, 'Number of synapses = %d');
%   synnum = zeros(1, N);
% 
%   for i=1 : 1 : Nsyn
%       buf = fscanf(fid, '%s', 1);
%       while(~strcmp(buf,'from'))
%           buf = fscanf(fid, '%s', 1);
%       end
%       buf00 = fscanf(fid, '%d', 1);
%       synnum(buf00+1) = synnum(buf00+1) + 1;
%   end
%   fclose(fid);
%   fprintf('Connections loaded.\n');
%   
%   fprintf('\tSimulation loaded.\n\n');
  %Calculating stuff
%   V_average = zeros(1, number_of_exports);
%   for i=1 : 1 : number_of_exports
%       for j=1 : 1 : N
%         V_average(i) =  V_average(i) + V(i, j);
%       end
%       V_average(i) = V_average(i) / N;
%   end
%   
%   fprintf('Average calculated.\n');
  
  %Activity
  time_int_spikes = 1;
  integrated_spikes_amount = floor(time_length / time_int_spikes)+1;
  integrated_spikes = zeros(1, integrated_spikes_amount);
  integrated_spikes_time = zeros(1, integrated_spikes_amount);
  for i=1 : 1 : num
      integrated_spikes(floor(rastr_time(i) / time_int_spikes)+1) = ...
          integrated_spikes(floor(rastr_time(i) / time_int_spikes)+1) + 1;
  end
  for i=1 : 1 : integrated_spikes_amount
      integrated_spikes_time(i) = time_int_spikes * (i-1);
  end
  integrated_spikes = integrated_spikes / N;
  
  fprintf('Spikes calculated.\n');
  
%   %Average InterSpike Interval vs time
%   ISI = 1./(integrated_spikes * N);
%   
%   %Overall ISI
%   spikes_in_neuron = ones(1, N);
%   for i=1 : 1 : num
%       spikes_in_neuron(rastr(i)+1) = spikes_in_neuron(rastr(i)+1) + 1;
%   end
%   overall_ISI = max_rastr_time ./ spikes_in_neuron;
%   max_ISI_hist = max(max(overall_ISI));
%   min_ISI_hist = 0;
%   N_ISI_hist = N / 5;
%   dISI_hist = (max_ISI_hist - min_ISI_hist) / N_ISI_hist;
%   hist_ISI = zeros(1, N_ISI_hist);
%   for i=1 : 1 : N
%       buf00 = floor(overall_ISI(i) / dISI_hist);
%       hist_ISI(buf00) = hist_ISI(buf00) + 1;
%   end
  
  %True ISI calculation
  spike_per_neuron = zeros(N, 1);
  for i=1 : 1 : length(rastr)
      neuronum = rastr(i)+1;
      neurotime = rastr_time(i);
      spike_per_neuron(neuronum, 1) = spike_per_neuron(neuronum, 1) + 1;
      spike_per_neuron(neuronum, spike_per_neuron(neuronum, 1)+1) = neurotime;
  end
  diff_spike_per_neuron = zeros(size(spike_per_neuron));
  max_true_ISI = 0;
  min_true_ISI = 0;
  for i=1 : 1 : N
      diff_spike_per_neuron(i, 1) = spike_per_neuron(i, 1);
      diff_spike_per_neuron(i, 2) = spike_per_neuron(i, 2);
      if diff_spike_per_neuron(i, 2) > max_true_ISI
          max_true_ISI = diff_spike_per_neuron(i, 2);
      end
      if spike_per_neuron(i,1) > 1
          for j=2 : 1 : spike_per_neuron(i, 1)
              diff_spike_per_neuron(i, j+1) = spike_per_neuron(i, j+1) - spike_per_neuron(i, j);
              if diff_spike_per_neuron(i, j+1) > max_true_ISI
                  max_true_ISI = diff_spike_per_neuron(i, j+1);
              end
          end
      end
  end
  N_true_ISI = max_true_ISI;
  true_ISI = zeros(1, N_true_ISI+1);
  d_true_ISI = (max_true_ISI - min_true_ISI) / N_true_ISI;
  for i=1 : 1 : N
      if diff_spike_per_neuron(i,1) > 0
          for j=1 : 1 : diff_spike_per_neuron(i, 1)
              buf00 = floor(diff_spike_per_neuron(i, j+1) / d_true_ISI)+1;
              true_ISI(buf00) = true_ISI(buf00) + 1;
          end
      end
  end
  true_ISI_ind = min_true_ISI : d_true_ISI : max_true_ISI;
  
  %Burst calculating
  burst_threshold = 0.2;    % % of network
  window_size = 100;    %ms
  find_bursts = zeros(1, integrated_spikes_amount);
  for i=1 : 1 : integrated_spikes_amount
      if integrated_spikes(i) > burst_threshold
          find_bursts(i) = integrated_spikes(i);
      end
  end
  i = 1;
  while i < integrated_spikes_amount
      j = i+1;
      while integrated_spikes_time(j) - integrated_spikes_time(i) < window_size && ...
              j+1 <= integrated_spikes_amount
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
  for i=1 : 1 : integrated_spikes_amount
      if find_bursts(i)>0
          amount_of_bursts = amount_of_bursts + 1;
          burst_amps(amount_of_bursts) = find_bursts(i) * N;
          burst_time(amount_of_bursts) = integrated_spikes_time(i);
      end
  end
  fprintf('Activity analyzed.\n');
  
  average_burst_amp = 0;
  for i=1 : 1 : amount_of_bursts
      average_burst_amp = average_burst_amp + burst_amps(i);
  end
  average_burst_amp = average_burst_amp / amount_of_bursts;
  sd_burst_amp = 0;
  for i=1 : 1 : amount_of_bursts
      sd_burst_amp = sd_burst_amp + (average_burst_amp - burst_amps(i))^2;
  end
  sd_burst_amp = sqrt(sd_burst_amp/amount_of_bursts);
  fprintf('Normal bursts calculated.\n');
  
  %InterBurst Interval calculation
  IBI = zeros(size(burst_amps));
  IBI(1) = burst_time(1);
  for i=2 : 1 : amount_of_bursts
      IBI(i) = burst_time(i) - burst_time(i-1);
  end
  fprintf('InterBurst Interval calculated.\n');
  
%   burst_full_amps = zeros(size(burst_time));
%   for i=1 : 1 : length(rastr)
%       for j=1 : 1 : length(burst_time)
%           if burst_time(j) - rastr_time(i) > -time_int_spikes && ...
%                   burst_time(j) - rastr_time(i) <= 0
%               burst_full_amps(j) = burst_full_amps(j) + synnum(rastr(i)+1);
%           end
%       end
%   end
%   average_full_burst_amp = 0;
%   for i=1 : 1 : length(burst_time)
%       average_full_burst_amp = average_full_burst_amp + burst_full_amps(i);
%   end
%   average_full_burst_amp = average_full_burst_amp / length(burst_time);
%   sd_full_burst_amp = 0;
%   for i=1 : 1 : amount_of_bursts
%       sd_full_burst_amp = sd_full_burst_amp + (average_full_burst_amp - burst_full_amps(i))^2;
%   end
%   sd_full_burst_amp = sqrt(sd_full_burst_amp/length(burst_time));
%   fprintf('Synapse bursts calculated.\n');
  
%   %Plotting data
%   figure(1);
%   plot(time, V_average);
%   title('Average potential');
%   xlabel('Time, ms');
%   ylabel('Potential, mV');
  
  figure(2);
  hold on;
  plot(rastr_time, rastr, 'k.','MarkerSize',1);
  for i=1 : 1 : amount_of_bursts
      plot([burst_time(i) burst_time(i)], [0 N], 'r');
  end
  plot([max_rastr_time max_rastr_time], [0 N], 'g.-');
  title('Raster plot of spiking neurons');
  axis([0 time_length 0 N]);
  xlabel('Time, ms');
  ylabel('Number of neuron');
  hold off;
  
  figure(3);
  hold on;
  plot(integrated_spikes_time, integrated_spikes);
  plot([max_rastr_time max_rastr_time], [0 N], 'g.-');
  for i=1 : 1 : amount_of_bursts
      plot([burst_time(i) burst_time(i)], [0 burst_amps(i)/N], 'r');
  end
  title('Fraction of neurons spiking per time');
  axis([0 time_length 0 1]);
  xlabel('Time, ms');
  ylabel('Fraction of network');
  hold off;
  
%   figure(4);
%   hold on;
%   title('Average x, y, z, u.');
%   plot(time2, X, 'b');
%   plot(time2, Y, 'k');
%   plot(time2, Z, 'y');
%   plot(time2, U, 'r');
%   axis([0 time_length 0 1]);
%   legend('X', 'Y', 'Z', 'U');
%   xlabel('Time, ms');
%   ylabel('Fraction');
%   hold on;

%   figure(5);
%   hold on;
%   plot(burst_time, IBI, 'k');
%   title('InterBurst Interval');
% %   axis([0 time_length 0 2000]);
%   xlabel('Time, ms');
%   ylabel('InterBurst Interval, ms');
%   hold off;

%   figure(6);
%   hold on;
%   plot(integrated_spikes_time, ISI, 'k');
%   title('Average InterSpikes Interval');
%   axis([0 time_length 0 1]);
%   xlabel('Time, ms');
%   ylabel('InterSpikes Interval, ms');
%   hold off;

%   figure(7);
%   hold on;
%   plot(1 : 1 : length(hist_ISI), hist_ISI, 'k');
%   title('InterSpikes Interval Histogram');
% %   axis([0 time_length 0 50]);
%   xlabel('Time, ms');
%   ylabel(' ');
%   hold off;

  figure(8);
  hold on;
  plot(true_ISI_ind, true_ISI, 'k');
  title('InterSpikes Interval Histogram')
  axis([min_true_ISI max_true_ISI 0 1000]);
  xlabel('Time, ms');
  ylabel(' ');
  hold off;

  fprintf('Burst info: Bursts:');
  for i=1 : 1 : amount_of_bursts
      fprintf(' %d', burst_amps(i));
  end
  fprintf('\n');
  fprintf('Mean: %f, sd: %f.', average_burst_amp, sd_burst_amp);
%   fprintf('\nBurst unnormal info. Mean: %f, sd: %f.', average_full_burst_amp, sd_full_burst_amp);
  
  fprintf('\n');
end

