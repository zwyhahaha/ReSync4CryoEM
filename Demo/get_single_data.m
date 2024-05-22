function [Dist,time]=get_single_data(use_real,n,p,SNR,round,i,is_random)
    if use_real
        exp_name = ['SNR_' num2str(1/SNR) 'K_' num2str(K) '_' num2str(round)];
    else
        exp_name = ['P_' num2str(p) 'K_' num2str(n) 'filter_0.1_' num2str(round)];
    end
    if is_random
        dist_filename = strrep('Data/plot_BSGD/plot_data/dist.mat','dist',exp_name);
        csv_filename = strrep('Data/plot_BSGD/csv/demo.csv', 'demo', exp_name);
    else
        dist_filename = strrep('Data/plot_benchmarks/plot_data/dist.mat','dist',exp_name);
        csv_filename = strrep('Data/plot_benchmarks/csv/demo.csv', 'demo', exp_name);
    end
    if exist(dist_filename, 'file') == 2
        load(dist_filename,'all_dist');
        resultsTable = readtable(csv_filename);
        Dist = all_dist{i};
        time = resultsTable{i+1,'Time'};
    else
        disp('check setting!');
    end
end

