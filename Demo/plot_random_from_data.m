function plot_random_from_data(use_real,n,p,SNR,rounds)
if use_real
    fprintf('Start loading... SNR=%0.2f, N=%d \n \n',SNR,n); 
    exp_name = ['SNR_' num2str(1/SNR) 'K_' num2str(K)];
else
    fprintf('Start testing... p=%0.2f, N=%d \n \n',p,n); 
    exp_name = ['P_' num2str(p) 'K_' num2str(n)];
end

markerList = {'o', '+', '*', '^', 'x', '_', '|', '.', 'v', '<', '>', 'p', 'h'};
len = 5;
for i=1:len
    Dists = cell(rounds,1);
    times = 0;
    for round = 1:rounds
        [Dists{round},time] = get_single_data(use_real,n,p,SNR,round,i,1);
        times = times + time;
    end
    Dist = array_avg(Dists);
    time = times/rounds;
    k = length(Dist);
    x = linspace(0, time, k);
%     semilogy(log(1000*x),Dist,['-', markerList{i}],'LineWidth',2,'MarkerIndices', 1:10:k,'MarkerSize',8);
    loglog(x,Dist,['-', markerList{i}],'LineWidth',2,'MarkerIndices', 1:10:k,'MarkerSize',8);
    hold on
end

algorithms = ["$\rho=$1","$\rho=$0.5","$\rho=$0.2","$\rho=$0.1","$\rho=$0.05"];
hold off
set(gcf, 'Color', 'white');
set(gca, 'LineWidth' , 1.7, 'FontName', 'Times New Roman','FontSize',18);
legend(algorithms,'Interpreter','latex','FontName','Times New Roman','FontSize',20,'Location','NorthEast')
xlabel('Time (sec)','Interpreter','latex','FontName','Times New Roman','FontSize',20)
ylabel('MSE','Interpreter','latex','FontName','Times New Roman','FontSize',20)

fig_filename = strrep('plot/plot_random/my_fig.fig','my_fig',exp_name);
savefig(fig_filename);

end

