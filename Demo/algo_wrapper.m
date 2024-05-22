function [Dist,time,iter] = algo_wrapper(Cij, R_init, R_orig, ReSync_parameters,solver)
n = size(Cij, 3);
filter_ratio = ReSync_parameters.filter_ratio;
p = ReSync_parameters.p;
ReSync_parameters.filter_size = round(filter_ratio * n);
if strcmpi(solver, 'ReSync')
    if isfield(ReSync_parameters,'target'), ReSync_parameters.target = ReSync_parameters.target; end
    ReSync_parameters.stepsize = 0.1/ (n*p);
    ReSync_parameters.decay = 0.9;
    tic;
    [~, Dist, mse_time, iter] = ReSync(Cij, R_init, R_orig, ReSync_parameters);
    time = toc;
    time = time - mse_time;
    fprintf('full time: %0.2f \n', time); 
elseif strcmpi(solver, 'ReSync_norm')
    if isfield(ReSync_parameters,'target'), ReSync_parameters.target = ReSync_parameters.target; end
    ReSync_parameters.stepsize = 0.1/ (n*p);
    ReSync_parameters.decay = 0.9;
    tic;
    [~, Dist, mse_time, iter] = ReSync_norm(Cij, R_init, R_orig, ReSync_parameters);
    time = toc;
    time = time - mse_time;
    fprintf('full time: %0.2f \n', time); 
elseif strcmpi(solver, 'SDR-ADMM')
    pars.alpha = 2/3;
    pars.solver = 'ADMM';
    [R_ADMM,time] = est_orientations_LUD_C(Cij,pars);
    mse_admm = cmpt_mse(R_ADMM, R_orig);
    Dist = [mse_admm];
    iter = 1;
    fprintf('SDP_ADMM estimation error: %d \n', mse_admm);
    fprintf('SDR-ADMM time: %0.2f \n', time); 
elseif strcmpi(solver, 'LS-PGD')
    tic;
    [~,  Dist, mse_time, iter]= my_ProjGradRotLS(Cij(1:2,:,:), R_init, R_orig, ReSync_parameters); % LS-PGD
    time = toc;
    time = time - mse_time;
    fprintf('LS_PGD time: %0.2f \n', time); 
elseif strcmpi(solver, 'LUD-PGD')
    tic;
    [~,  Dist, mse_time, iter]= my_ProjGradRotLUD(Cij(1:2,:,:), R_init, R_orig, ReSync_parameters); % LUD-PGD
    time = toc;
    time = time - mse_time;
    fprintf('LUD_PGD time: %0.2f \n', time); 
elseif strcmpi(solver, 'LUD-IRLS-PGD')
    tic;
    [~,  Dist, mse_time, iter]= my_ProjGradRotIterWLS(Cij(1:2,:,:), R_init, R_orig, ReSync_parameters); % LUD-IRLS-PGD
    time = toc;
    time = time - mse_time;
    fprintf('LUD_PGD time: %0.2f \n', time); 
elseif strcmpi(solver, 'shuffle-BCD')
    if isfield(ReSync_parameters,'target'), ReSync_parameters.target = 5*ReSync_parameters.target; end
    ReSync_parameters.stepsize = 0.1/ (n*p);
    ReSync_parameters.decay=0.98^(1/filter_ratio);
    tic;
    [~, Dist, mse_time, iter] = shuffle_BCD_ReSync(Cij, R_init, R_orig, ReSync_parameters);
    time = toc;
    time = time - mse_time;
    fprintf('shuffle_BCD time: %0.2f \n', time); 
elseif strcmpi(solver, 'shuffle-SGD')
    if isfield(ReSync_parameters,'target'), ReSync_parameters.target = 5*ReSync_parameters.target; end
    ReSync_parameters.stepsize = 0.01/ (n*p*filter_ratio);
    ReSync_parameters.decay = (1-0.2*filter_ratio)^(1/filter_ratio);
    ReSync_parameters.check_freq = 3;
    tic;
    [~, Dist, mse_time, iter] = shuffle_SGD_ReSync(Cij, R_init, R_orig, ReSync_parameters);
    time = toc;
    time = time - mse_time;
    fprintf('shuffle_SGD time: %0.2f \n', time);
elseif strcmpi(solver, 'shuffle-block(homo)')
    if isfield(ReSync_parameters,'target'), ReSync_parameters.target = 5*ReSync_parameters.target; end
    ReSync_parameters.stepsize = 0.1/ (n*p*filter_ratio);
    ReSync_parameters.decay = (1-0.2*filter_ratio)^(1/filter_ratio);
    ReSync_parameters.use_homofilter=1;
    ReSync_parameters.check_freq = 1;
    tic;
    [~, Dist, mse_time, iter] = shuffle_block_ReSync(Cij, R_init, R_orig, ReSync_parameters);
    time = toc;
    time = time - mse_time;
    iter=ceil(iter*filter_ratio);
    fprintf('shuffle_block time: %0.2f \n', time); 
else
    disp('unsupported algorithm');
end
end

