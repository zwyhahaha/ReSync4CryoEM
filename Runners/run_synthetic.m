clear all; close all; clc;

%% for formal experiments, run this setting
rounds = 10;
filter_ratio = 0.1;
N =  [1000, 3000, 5000];
P = [0.5, 0.3, 0.1, 0.05];

%% for quick execution, run this setting
rounds = 1;
filter_ratio = 0.1;
N =  [500];
P = [0.5, 0.3, 0.1, 0.05];
%% main
parpool(length(P))

parfor i=1:length(P)
    p = P(i);
    for n=N
        for round = 1:rounds
           compare_synthetic(n,p,filter_ratio,round)
        end
    end
end
    
delete(gcp('nocreate'))