function [R_est,Dist,stepTime, Witer]= my_ProjGradRotIterWLS(C, R_init, R_orig, Param)
% 2019-10-14
%% Object:
%  \sum_{ij} \min 1/2*wij*|| RiCij - RjCji||_2^2
%  s.t  Ri'Ri = I
%  Ri^{k+1} = Ri - t*wij*(RiCij - RjCji)*Cij^{T}
%% Proximal Gradient with Backtracking Line-Search
fprintf('Starting LUD_IRLS... \n');

Rflag = 2;
K  = size(C,3);
R_old = R_init;
MaxIter   = Param.max_iter;
stop_threshold = Param.stop_threshold;
SolRE = stop_threshold;
check_freq = Param.check_freq;
if ~exist('R_old','var')    % initial the rotation matrix value
    if Rflag == 1
        R_old   = rand(3,3,K);
        for i = 1:K
            [U,~,V]   = svd(R_old(:,:,i));
            R_old(:,:,i) = U * V';
        end
        MSEiter1 = []; REiter1 = [];
    else
        SolRE = 1e-3;
        [R_old,  MSEiter1, REiter1]= ProjGradRotLS(C, Param);
    end
end

W     = ones(K,K);

TOL   = 1.0e-14;
Dist = [];


tk    = 1; %[0.99,0.7,0.618,0.5,0.2,0.08]
alpha = 0.618; %0.618;%0.99; % SNR = 1/16 alpha =0.618
R_new = R_old;

stepTime = 0;
Witer = 0;
for iter = 1:MaxIter
    %     flag = 1;
    %     kk   = 0;
    %     Param1  = Param;
    %     Param1.InitRot = R_old;  Param1.MaxIter = 10; Param1.SolRE = 1e-4;
    %
    %
    %
    %     [R_new,  MSEiter2, REiter2, W]= ProjGradRotWLS(C,W, Param1);
    %     MSEiter = [MSEiter MSEiter2];
    %     REiter  = [REiter REiter2];
    
    
    R_old2 = R_old;
    clear MSEiter2 REiter2 ObjValue
    
    for ii = 1:100
        Witer = Witer + 1;
        if Witer > MaxIter, break; end
        flag = 1;
        kk   = 0;
        [GradRi,JRK,B] = FunCostAndGradp2q2w(R_old(:,1:2,:),C,W);
        Grad_Ri      = zeros(3,3,K);
        Grad_Ri(:,1:2,:) = GradRi;
        
        tk = min(0.618,tk*2.5);
        while flag
            kk = kk+1;
            GradObj = 0;
            for i =1:K
                temp    = R_old(:,:,i) - tk * Grad_Ri(:,:,i);
                [U,~,V] = svd(temp,0);
                R_new(:,:,i) = U*V';
                
                tmp     = Grad_Ri(:,:,i).* (R_new(:,:,i)-R_old(:,:,i));
                GradObj = GradObj + sum(tmp(:)); % Grad_R*(Rnew-Rold)
            end
            
            [~,JR,B] = FunCostAndGradp2q2w(R_new(:,1:2,:),C,W);
            JRKnew = JRK + GradObj + (norm(R_new(:)-R_old(:))^2)/(2*tk);
            
            if JR > JRKnew && kk<20
                tk = alpha*tk;
            else
                flag = 0;
            end
        end
        %fprintf('step=%2e\n',tk);
        %W = B;
        REiter2(ii) = norm(R_new(:) - R_old(:))/norm(R_new(:));
        
        R_old = R_new;
        
        ObjValue(ii) = JR;
        
        %% Make sure that we got true rotations.
        if exist('R_orig','var')
            R_est = zeros(3,3,K);
            for k=1:K
                R_est(:,:,k) = [R_new(:,1:2,k),cross(R_new(:,1,k),R_new(:,2,k))];
                R = R_est(:,:,k);
                erro = norm(R*R.'-eye(3));
                if erro > TOL || abs(det(R)-1)> TOL
                    [U,~,V] = svd(R);
                    R_est(:,:,k) = U*V.';
                end
            end
%             [MSEiter2(ii),~,~] = check_MSE(est_rots, truerot);
%             MSEiter2(ii) = cmpt_mse(R_est, R_orig);
        end
        if REiter2(ii) <1e-4,  break; end
        if mod(Witer, check_freq) == 0
           stepStartTime = tic;
           mse = cmpt_mse(R_est, R_orig);
           stepTime = stepTime + toc(stepStartTime);
           fprintf('estimation error: %d in iteration %d \n',mse, Witer);
           Dist = [Dist, mse];
    %        if REiter <SolRE, fprintf('iter = %d\n',iter); break; end
        end
        %if iter > 2, break; end
    end
%     Witer = Witer + ii;
    W   = B;
%     MSEiter = [MSEiter MSEiter2];
%     REiter  = [REiter REiter2];
%     tmp   = R_new(:,1:2,:) - R_old2(:,1:2,:);
    
%     REiterb  = norm(tmp(:))/norm(R_new(:));
    
    R_old = R_new;
    
%     if REiterb  <SolRE, fprintf('iter = %d\n',iter); break; end
    mse = cmpt_mse(R_est, R_orig);    
    tmp   = R_new(:,1:2,:) - R_old2(:,1:2,:);
    REiterb  = norm(tmp(:))/norm(R_new(:));   
    if isfield(Param,'target')
        if mse < Param.target, break; end
    else
        if  REiterb<SolRE, break; end
    end
    
end

R_est = zeros(3,3,K);
for k=1:K
    R_est(:,:,k) = [R_new(:,1:2,k),cross(R_new(:,1,k),R_new(:,2,k))];
    R = R_est(:,:,k);
    erro = norm(R*R.'-eye(3));
    if erro > TOL || abs(det(R)-1)> TOL
        [U,~,V] = svd(R);
        R_est(:,:,k) = U*V.';
    end
   
end

% figure(9); subplot(2,1,1); semilogy(REiter); subplot(2,1,2); semilogy(MSEiter);
%
% MSEiter = [MSEiter1 MSEiter];
% REiter  = [REiter1 REiter];