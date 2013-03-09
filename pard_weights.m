function pard_weights
    fprintf('\t\tNew instance %f\n', rand());
    
    weightmin = 0;
    weightmax = 1;
    histparts = 101;
    weight_histograms = [1200 1250 1750 1800 2730 2780 2950 3000 4870 4920];
    for i=1 : 1 : histparts
        weight_index(i) = (i - 1) * (weightmax - weightmin) / (histparts - 1);
    end
    
    i = 1;
    j = 1;
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
        
        R(i) = 0;
        if(i>1)
            for k=1 : 1 : length(arr)
                R(i) = R(i) + (arr2(k) - arr(k))^2;
            end
        end
        arr2 = arr;
        
        disp(time(i));
        i = i + 1;
%         if(time(i-1)>max(weight_histograms)+50)
%             break;
%         end
    end
    
    figure(1);
    errorbar(time, mean, sd);
    axis([0 max(time) -0.1 1.1]);
    xlabel('Time, ms');
    ylabel('Average weight');
    
    weight_hist_data = weight_hist_data / length(arr);
    max_hist = 1.2 * max(max(weight_hist_data));
    min_hist = 0;
    
%     for i=1 : 1 : length(weight_times)
%         figure(i+1);
%         plot(weight_index, weight_hist_data(i, :), 'k*');
%         title('Histogram of weights'); 
%         xlabel('Weight');
%         ylabel('Part of synapses with this weight');
%         axis([weightmin weightmax min_hist max_hist]);
%     end
%     
%     figure(length(weight_times)+2);
%     A = 8;
%     B = 9;
%     min_diff_hist = (min(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     max_diff_hist = (max(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     plot(weight_index, weight_hist_data(B, :) - weight_hist_data(A, :), 'k*');
%     title('Histogram of difference of weights'); 
%     xlabel('Weight');
%     ylabel('Part of synapses with this weight');
%     axis([weightmin weightmax min_diff_hist max_diff_hist]);
%     
%     figure(length(weight_times)+3);
%     A = 9;
%     B = 10;
%     min_diff_hist = (min(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     max_diff_hist = (max(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     plot(weight_index, weight_hist_data(B, :) - weight_hist_data(A, :), 'k*');
%     title('Histogram of difference of weights'); 
%     xlabel('Weight');
%     ylabel('Part of synapses with this weight');
%     axis([weightmin weightmax min_diff_hist max_diff_hist]);
%     
%     figure(length(weight_times)+4);
%     A = 8;
%     B = 9;
%     max_diff_hist = abs(max(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     min_diff_hist = 0;
%     plot(weight_index, abs(weight_hist_data(B, :) - weight_hist_data(A, :)), 'k*');
%     title('Histogram of difference of weights'); 
%     xlabel('Weight');
%     ylabel('Part of synapses with this weight');
%     axis([weightmin weightmax min_diff_hist max_diff_hist]);
%     
%     figure(length(weight_times)+5);
%     A = 9;
%     B = 10;
%     max_diff_hist = abs(max(weight_hist_data(B, :) - weight_hist_data(A, :)))*1.2;
%     min_diff_hist = 0;
%     plot(weight_index, abs(weight_hist_data(B, :) - weight_hist_data(A, :)), 'k*');
%     title('Histogram of difference of weights'); 
%     xlabel('Weight');
%     ylabel('Part of synapses with this weight');
%     axis([weightmin weightmax min_diff_hist max_diff_hist]);
    
    figure(length(weight_times)+5);
    plot(time, R);
    
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