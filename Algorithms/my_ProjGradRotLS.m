function [R_est,Dist,stepTime, iter]= my_ProjGradRotLS(C, R_init, R_orig, Param)
% 2019-10-14
%% Object:
%  \sum_{ij} \min 1/2*|| RiCij - RjCji||_2^2
%  s.t  Ri'Ri = I
%  Ri^{k+1} = Ri - t*(RiCij - RjCji)*Cij^{T}
%% Proximal Gradient with Backtracking Line-Search
fprintf('Starting LS_PGD... \n');

n  = size(C,3);
R_old = R_init;
MaxIter   = Param.max_iter;
stop_threshold = Param.stop_threshold;
SolRE = stop_threshold;
check_freq = Param.check_freq;

if ~exist('R_old','var')    % initial the rotation matrix value
    R_old   = rand(3,3,n);
    for i = 1:n
        [U,~,V]   = svd(R_old(:,:,i));
        R_old(:,:,i) = U * V';
    end
end


TOL   = 1.0e-14;
Dist = [];

tk    = 1; %[0.99,0.7,0.618,0.5,0.2,0.08]
alpha = 0.618; %0.618;%0.99; % SNR = 1/16 alpha =0.618
R_new = R_old;

% totalStartTime = tic;
stepTime = 0;
for iter = 1:MaxIter
    flag = 1;
    kk   = 0;
    [GradRi,JRK] = FunCostAndGradp2q2(R_old(:,1:2,:),C);
    Grad_Ri      = zeros(3,3,n);
    Grad_Ri(:,1:2,:) = GradRi;
    
    tk = min(0.618,tk*2.5);
    while flag
        kk = kk+1;
        GradObj = 0;
        for i =1:n
            temp    = R_old(:,:,i) - tk * Grad_Ri(:,:,i);
            [U,~,V] = svd(temp,0);
            R_new(:,:,i) = U*V';
            
            tmp     = Grad_Ri(:,:,i).* (R_new(:,:,i)-R_old(:,:,i));
            GradObj = GradObj + sum(tmp(:)); % Grad_R*(Rnew-Rold)
        end
        
        [~,JR] = FunCostAndGradp2q2(R_new(:,1:2,:),C);
        JRKnew = JRK + GradObj + (norm(R_new(:)-R_old(:))^2)/(2*tk);
        
        if JR > JRKnew && kk<20
            tk = alpha*tk;
        else
            flag = 0;
        end
    end
    %fprintf('step=%2e\n',tk);
    
    

    %% Make sure that we got true rotations.
    if exist('R_orig','var')        
        R_est = zeros(3,3,n);
        for k=1:n
            R_est(:,:,k) = [R_new(:,1:2,k),cross(R_new(:,1,k),R_new(:,2,k))];
            R = R_est(:,:,k);
            erro = norm(R*R.'-eye(3));
            if erro > TOL || abs(det(R)-1)> TOL
                [U,~,V] = svd(R);
                R_est(:,:,k) = U*V.';
            end
        end
        
    end
    if mod(iter, check_freq) == 0
       stepStartTime = tic;
       mse = cmpt_mse(R_est, R_orig);
       stepTime = stepTime + toc(stepStartTime);
       fprintf('estimation error: %d in iteration %d \n',mse, iter); 
       Dist = [Dist, mse];
       REiter = norm(R_new(:) - R_old(:))/norm(R_new(:));
       if isfield(Param,'target')
            if mse < Param.target, break; end
        else
            if REiter <SolRE, break; end
       end
    end
    R_old = R_new;
        
%     ObjValue(iter) = JR; 
end

R_est = zeros(3,3,n);
for k=1:n
    R_est(:,:,k) = [R_new(:,1:2,k),cross(R_new(:,1,k),R_new(:,2,k))];
    R = R_est(:,:,k);
    erro = norm(R*R.'-eye(3));
    if erro > TOL || abs(det(R)-1)> TOL
        [U,~,V] = svd(R);
        R_est(:,:,k) = U*V.';
    end
end

% totalTime = toc(totalStartTime);
% finalTime=totalTime-stepTime;
% figure(9); subplot(2,1,1); semilogy(REiter); subplot(2,1,2); semilogy(MSEiter);