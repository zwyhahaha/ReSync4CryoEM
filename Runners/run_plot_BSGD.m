
rounds = 10;
filter_ratios = [1.0,0.5,0.2,0.1,0.05];
N =  [3000];

P = [0.1, 0.3, 0.5];
parpool(length(P))

parfor i=1:length(P)
    p = P(i);
    for n=N
        for round = 1:rounds
           plot_block_ReSync(n,p,filter_ratios,round)
        end
    end
end

delete(gcp('nocreate'))
