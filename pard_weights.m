function pard_weights
global pi
pi = 3.14159265359;
    fprintf('\t\tNew instance %5.0f\n', random('u',1,99999));
    %DO NOT FORGET to set right amount of neurons in data
    N = 1000;    
    
    %This part is about w
    weightmin = 0;
    weightmax = 1;
    histparts = 101;    %Intervals in histogram
    weight_histograms = [0];    %Where to plot histograms (time)
    %You can set something like this: [0 500 1500 3000] and script will plot
    %them in increasing order
    weight_index = weightmin : (weightmax - weightmin) / (histparts - 1) : weightmax;
    
    %This part is about A*w
    Aweightmin = -200;  
    Aweightmax = 200;
    Ahistparts = 101;   %Intervals in histogram
    Aweight_histograms = [0];   %where to plot histograms (time)
    Aweight_index = Aweightmin : (Aweightmax - Aweightmin) / (Ahistparts - 1) : Aweightmax;
    
    %Finding A from export/synapse.txt
    A = zeros(1, N*N);
    fid = fopen('export/synapses.txt', 'r');
    Nsyn = fscanf(fid, 'Number of synapses = %d');
    for i=1 : 1 : Nsyn
        num = fscanf(fid,'%d');
        from = fscanf(fid,': from %d');
        to = fscanf(fid,' to %d');
        fscanf(fid,'%s', 11);
        buf00 = fscanf(fid,'%f', 1);
        A(1, num+1) = buf00;
        fscanf(fid, '%s', 1);
    end
    fclose(fid);
    
    i = 1;
    j = 1;
    Aj = 1;
    fid = fopen('data/weights.txt', 'r');
    while(1)
        buf1 = fscanf(fid, '%c', 1);
        if buf1 == '_'
        else
            break;
        end
        time(i) = fscanf(fid, '________Time = %f', 1);
        arr = fscanf(fid, '%f');
        
        [mean(i) sd(i)] = findMeanDropNull(arr, length(arr));
        
        if isItHere(time(i), weight_histograms)
            line = makeHistogramData(arr, length(arr), histparts, weightmin, weightmax);
            weight_times(j) = time(i);
            for k = 1 : 1 : length(line)
                weight_hist_data(j, k) = line(k);
            end
            j = j + 1;
        end
        
        if isItHere(time(i), Aweight_histograms)
            Aarr = arr.*(A');
            Aline = makeHistogramData(Aarr, length(Aarr), Ahistparts, Aweightmin, Aweightmax);
            Aweight_times(Aj) = time(i);
            for k = 1 : 1 : length(Aline)
                Aweight_hist_data(Aj, k) = Aline(k);
            end
            Aj = Aj + 1;
        end
        
%         R(i) = 0;
%         if(i>1)
%             for k=1 : 1 : length(arr)
%                 R(i) = R(i) + (arr2(k) - arr(k))^2;
%             end
%         end
%         arr2 = arr;
        
        disp(time(i));
        i = i + 1;
        if(time(i-1)>max(weight_histograms)+50 && time(i-1)>max(Aweight_histograms)+50)
            break;
        end
    end
    
    figure(1);
    errorbar(time, mean, sd);
%     axis([0 max(time) -0.1 1.1]);
    xlabel('Time, ms');
    ylabel('Average weight');
    weight_hist_data = weight_hist_data / length(arr);
    max_hist = 1.2 * max(max(weight_hist_data));
    min_hist = 0;
    
    for i=1 : 1 : length(weight_times)
        figure(i+1);
        hold on
        plot(weight_index, weight_hist_data(i, :), 'k*');
        title('Histogram of weights'); 
        xlabel('Weight');
        ylabel('Part of synapses with this weight');
        axis([weightmin weightmax min_hist max_hist]);
        hold off;
    end
    
    Aweight_hist_data = Aweight_hist_data / length(arr);
    Amax_hist = 1.2 * max(max(Aweight_hist_data));
    Amin_hist = 0;
    for i=1 : 1 : length(Aweight_times)
        figure(i+1+length(weight_times));
        hold on;
        plot(Aweight_index, Aweight_hist_data(i, :), 'k*');
        title('Histogram of weights'); 
        xlabel('Weight');
        ylabel('Part of synapses with this weight');
        arr = Aweightmin : (Aweightmax - Aweightmin) / 1000 : Aweightmax;
        plot(arr,(0.64 .* normalDistrTheor(38, 38/2, arr) + ...
                  0.16 .* normalDistrTheor(54, 54/2, arr) + ...
                  0.20 .* normalDistrTheor(-72, 72/2, arr) ) * Amax_hist * 50 / 1000, 'r');
        legend('Experimental distribution', 'Theoretical');
        axis([Aweightmin Aweightmax Amin_hist Amax_hist]);
        hold off;
    end
    
%     figure(length(weight_times)+length(Aweight_times) + 1);
%     plot(time, R);

%     figure(length(weight_times)+length(Aweight_times) + 2);
%     A = 8;
%     B = 9;
%     min_diff_hist = (min(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     max_diff_hist = (max(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     plot(weight_index, weight_hist_data(B, :) - weight_hist_data(A, :), 'k*');
%     title('Histogram of difference of weights'); 
%     xlabel('Weight');
%     ylabel('Part of synapses with this weight');
%     axis([weightmin weightmax min_diff_hist max_diff_hist]);
    
    
    fclose(fid);
end

function [mean sd] = findMeanDropNull(A, N)
    mean = 0;
    num = 0;
    sd = 0;
    
    for i=1 : 1 : N
        mean = mean + A(i);
        if A(i) > 0
            num = num + 1;
        end
    end
    mean = mean / num;
    if num == 0
        mean = 0;
        return;
    end
    
    for i=1 : 1 : N
        if A(i) > 0
            sd = sd + (mean - A(i))^2;
        end
    end
    if num>1
        sd = sqrt(sd / (num - 1));
    end
end
function line = makeHistogramData(A, N, parts, mini, maxi)
    if max(A) > maxi || min(A) < mini
        fprintf('\n\t\tRANGE ERROR\n');
        return;
    end
    
    line = zeros(1, parts);
    dx = (maxi - mini) / (parts - 1);
    
    for i=1 : 1 : N
        if A(i)~=0
            line(floor((A(i) - mini)/dx) + 1) = ...
                line(floor((A(i) - mini)/dx) + 1) + 1;
        end
    end
end
function answ = isItHere(elem, line)
    answ = 0;
    for i=1 : 1 : length(line)
        if line(i) == elem
            answ = 1;
        end
    end
end
function answ = normalDistrTheor(average, sd, index)
global pi
    answ = 1 ./ (sqrt(2*pi)*sd) .* exp(- (index - average).^2 ./ (2 * sd^2));
end