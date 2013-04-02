function pard_weights
    fprintf('\t\tNew instance %5.0f\n', random('u',1,99999));
    N = 10;
    
    weightmin = 0;
    weightmax = 1;
    histparts = 101;
    weight_histograms = [0];
    for i=1 : 1 : histparts
        weight_index(i) = (i - 1) * (weightmax - weightmin) / (histparts - 1);
    end
    
    Aweightmin = -150;
    Aweightmax = 150;
    Ahistparts = 10;
    Aweight_histograms = [0];
    for i=1 : 1 : Ahistparts
        Aweight_index(i) = (i - 1) * (Aweightmax - Aweightmin) / (Ahistparts - 1);
    end
    
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
            line = makeHistogramData(Aarr, length(Aarr), Ahistparts, Aweightmin, Aweightmax);
            Aweight_times(Aj) = time(i);
            disp(line);
            for k = 1 : 1 : length(line)
                Aweight_hist_data(Aj, k) = line(k);
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
%         if(time(i-1)>max(weight_histograms)+50)
%             break;
%         end
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
        plot(weight_index, weight_hist_data(i, :), 'k*');
        title('Histogram of weights'); 
        xlabel('Weight');
        ylabel('Part of synapses with this weight');
        axis([weightmin weightmax min_hist max_hist]);
    end
    
    Aweight_hist_data = Aweight_hist_data / length(arr);
    Amax_hist = 1.2 * max(max(Aweight_hist_data));
    Amin_hist = 0;
    for i=1 : 1 : length(Aweight_times)
        figure(i+30+length(weight_times));
        plot(Aweight_index, Aweight_hist_data(i, :), 'k*');
        title('Histogram of weights'); 
        xlabel('Weight');
        ylabel('Part of synapses with this weight');
        axis([Aweightmin Aweightmax Amin_hist Amax_hist]);
    end
    
%     figure(length(weight_times)+20);
%     A = 8;
%     B = 9;
%     min_diff_hist = (min(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     max_diff_hist = (max(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     plot(weight_index, weight_hist_data(B, :) - weight_hist_data(A, :), 'k*');
%     title('Histogram of difference of weights'); 
%     xlabel('Weight');
%     ylabel('Part of synapses with this weight');
%     axis([weightmin weightmax min_diff_hist max_diff_hist]);
    
%     figure(length(weight_times)+2);
%     plot(time, R);
    
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
        if A(i)>0
            line(floor(A(i)/dx) + 1) = line(floor(A(i)/dx) + 1) + 1;
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