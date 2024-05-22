clear all; close all; clc;


algorithms = ["LS-PGD","LUD-PGD","LUD-IRLS-PGD","ReSync"];
rounds = 10;
filter_ratio = 0.1;
N =  [3000];
P = [0.1, 0.3, 0.5];

parpool(length(P))

parfor i=1:length(P)
    p = P(i);
    for n=N
        for round = 1:rounds
           plot_benchmark(n,p,filter_ratio,round,algorithms)
        end
    end
end

delete(gcp('nocreate'))


