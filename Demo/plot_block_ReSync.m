function plot_block_ReSync(n,p,filter_ratios,round)
%%
n_theta=360;

fprintf('Start testing... p=%0.2f, N=%d \n \n',p,n); 
exp_name = ['P_' num2str(p) 'K_' num2str(n) '_' num2str(round)];
model_out = Uniform_Topology(n,p,n_theta);
R_orig = model_out.R_orig; % ground truth rotations (3 by 3 by n)
Cij = model_out.Cij; % given common lines (3 by n by n)
%%
[R_init, ~] = Eigenvector_Relaxation(Cij); % initialized rotations (3 by 3 by n)
tic;
mse_eig = cmpt_mse(R_init, R_orig);
time_eig=toc;
fprintf('estimation error: %d \n', mse_eig); 

% set ReSync defult parameters
ReSync_parameters.max_iter = 101;
ReSync_parameters.stepsize = 0.01/ (n*p);
ReSync_parameters.decay = 0.98;
ReSync_parameters.stop_threshold = -0.1;
ReSync_parameters.check_freq = 2;
fig = figure;
markerList = {'o', '+', '*', '^', 'x', '_', '|', '.', 'v', '<', '>', 'p', 'h'};

ReSync_parameters.p = p;

resultsTable = table();
resultsTable(1,:)={"Eig",mse_eig, 1, time_eig};

all_dist = cell(length(filter_ratios),1);
algorithms = strings(1,length(filter_ratios));

for i = 2:1+length(filter_ratios)
    ReSync_parameters.filter_ratio = filter_ratios(i-1);
    [Dist,time,iter]=algo_wrapper(Cij, R_init, R_orig, ReSync_parameters,"shuffle-block(homo)");
    all_dist{i-1} = Dist;
    algo_name = ['$\rho=$' num2str(filter_ratios(i-1))];
    algorithms(i-1) = algo_name;
    resultsTable(i,:)={algo_name,Dist(end), iter, time};
    n = length(Dist);
    x = linspace(0, time, n);
    loglog(x,Dist,['-', markerList{i-1}],'LineWidth',2,'MarkerIndices', 1:10:n,'MarkerSize',8);
    hold on
end
%% save results
hold off
set(gcf, 'Color', 'white');
set(gca, 'LineWidth' , 1.7, 'FontName', 'Times New Roman','FontSize',18);
legend(algorithms,'Interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','NorthEast')
xlabel('Time (sec)','Interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('MSE','Interpreter','latex','FontName','Times New Roman','FontSize',20)

fig_filename = strrep('Data/plot_BSGD/fig/my_fig.fig','my_fig',exp_name);
savefig(fig_filename);
dist_filename = strrep('Data/plot_BSGD/plot_data/dist.mat','dist',exp_name);
save(dist_filename,'all_dist');

resultsTable.Properties.VariableNames = {'Method', 'MSE', 'Iter', 'Time'};
disp(resultsTable);
csv_filename = strrep('Data/plot_BSGD/csv/demo.csv', 'demo', exp_name);
writetable(resultsTable, csv_filename);

end