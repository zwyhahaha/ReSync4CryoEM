function compare_real(n,SNR,filter_ratio,algorithms)
%% load data
K=n;
fprintf('Start testing... SNR=%0.2f, N=%d \n \n',SNR,n); 
exp_name = ['SNR_' num2str(1/SNR) 'K_' num2str(K) 'filter_' num2str(filter_ratio)];
C_save_dir = ['Data/cl_matrix/SNR_' num2str(1/SNR) 'K_' num2str(K) '_noisy_cl.mat'];
ref_rot_save_dir = ['Data/cl_matrix/SNR_' num2str(1/SNR) 'K_' num2str(K) '_ref_rot.mat'];
if exist(C_save_dir, 'file') == 2
    load(C_save_dir, 'C');
    load(ref_rot_save_dir, 'ref_rot');
    disp('succesfully load data!');
else
    [C, ref_rot, ~] = gen_and_save_cl(K, SNR);
    disp('file not exists, creating...');
end
Cij=C; R_orig=ref_rot; p=get_p_from_snr(SNR);

%% running algorithms
% initialized rotations
[R_init, ~] = Eigenvector_Relaxation(Cij); 
tic;
mse_eig = cmpt_mse(R_init, R_orig);
time_eig=toc;
fprintf('estimation error: %d \n', mse_eig); 

% set ReSync defult parameters
ReSync_parameters.max_iter = 201;
ReSync_parameters.stepsize = 0.01/ (n*p);
ReSync_parameters.decay = 0.98;
ReSync_parameters.stop_threshold = 1e-5;
ReSync_parameters.check_freq = 10;
ReSync_parameters.filter_ratio = filter_ratio;
ReSync_parameters.p = p;

resultsTable = table();
resultsTable(1,:)={"Eig",mse_eig, 1, time_eig};
all_dist = cell(length(algorithms),1);

for i = 2:1+length(algorithms)
    [Dist,time,iter]=algo_wrapper(Cij, R_init, R_orig, ReSync_parameters,algorithms(i-1));
    all_dist{i-1} = Dist;
    resultsTable(i,:)={algorithms(i-1),Dist(end), iter, time};
end
%% save results
resultsTable.Properties.VariableNames = {'Method', 'MSE', 'Iter', 'Time'};
disp(resultsTable);
csv_filename = strrep('Data/real_result/demo.csv', 'demo', exp_name);
writetable(resultsTable, csv_filename);

end

