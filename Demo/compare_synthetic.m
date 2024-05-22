function compare_synthetic(n,p,filter_ratio,round)
%% generate data
n_theta=360;
fprintf('Start testing... p=%0.2f, N=%d \n \n',p,n); 
exp_name = ['P_' num2str(p) 'K_' num2str(n) 'filter_' num2str(filter_ratio) '_' num2str(round)];
model_out = Uniform_Topology(n,p,n_theta);
R_orig = model_out.R_orig; % ground truth rotations (3 by 3 by n)
Cij = model_out.Cij; % given common lines (3 by n by n)
ReSync_parameters.target=get_target(n,p);

%% running algorithms
[R_init, ~] = Eigenvector_Relaxation(Cij); % initialized rotations (3 by 3 by n)
tic;
mse_eig = cmpt_mse(R_init, R_orig);
time_eig=toc;
fprintf('estimation error: %d \n', mse_eig); 

% set ReSync defult parameters
ReSync_parameters.max_iter = 201;
ReSync_parameters.stepsize = 0.01/ (n*p);
ReSync_parameters.decay = 0.98;
ReSync_parameters.stop_threshold = 1e-5;
ReSync_parameters.check_freq = 3;
ReSync_parameters.filter_ratio = filter_ratio;
ReSync_parameters.p = p;

resultsTable = table();
resultsTable(1,:)={"Eig",mse_eig, 1, time_eig};
[Dist,time,iter]=algo_wrapper(Cij, R_init, R_orig, ReSync_parameters,"LUD-PGD");
mse_PGD = Dist(end);
resultsTable(2,:)={"LUD-PGD",mse_PGD, iter, time};
ReSync_parameters.target = mse_PGD;
algorithms = ["ReSync","shuffle-BCD","shuffle-SGD","shuffle-block(homo)"];
ReSync_parameters.check_freq = 3;

for i = 3:2+length(algorithms)
    [Dist,time,iter]=algo_wrapper(Cij, R_init, R_orig, ReSync_parameters,algorithms(i-2));
    resultsTable(i,:)={algorithms(i-2),Dist(end), iter, time};
end
%% save results
resultsTable.Properties.VariableNames = {'Method', 'MSE', 'Iter', 'Time'};
disp(resultsTable);
csv_filename = strrep('Data/synthetic_result/demo.csv', 'demo', exp_name);
writetable(resultsTable, csv_filename);

end

