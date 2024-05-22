clear all; close all; clc;

%% for formal experiments, run this setting
%% Be careful that generating data needs a quite long time...
N =  [1000, 3000, 5000];
SNR = [1/32, 1/64, 1/128];

%% for quick execution, we recommend the following setting
N =  [500];
SNR = [1/32];

%% main
%% select benchmark algorithms
algorithms = ["SDR-ADMM", "ReSync-norm"]; % they require long running time, ignore them for quick execution
algorithms = ["LS-PGD","LUD-PGD","LUD-IRLS-PGD","ReSync","shuffle-BCD","shuffle-SGD","shuffle-block(homo)"];
filter_ratio = 0.1;

for snr = SNR
    for n = N
        compare_real(n,snr,filter_ratio,algorithms);
    end
end
