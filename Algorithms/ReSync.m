function [R_est, Dist, stepTime, iter] = ReSync(Cij, R_init, R_orig, ReSync_parameters)
% totalStartTime = tic;
% stop_threshold = 1e-5;
fprintf('Starting ReSync... \n');
    R_est = R_init;
    [~, ~, n] = size(R_init);
    Dist = [];
    stepTime = 0;
    for iter = 1 : ReSync_parameters.max_iter
        % linear decay of stepsize
        step = ReSync_parameters.stepsize * ReSync_parameters.decay^(iter);
        
        % compute the subgradient 
        [~,~, ~, x_egrad] = FunCostAndGradp2q1_ReSync(R_est(:,1:3,:),Cij(1:2,:,:));
        x_egrad = x_egrad*2;
        
        R_new = R_est - step * x_egrad;
        for i = 1 : n
            y = R_new(:,:,i);
            [Q_y, R_y] = qr(y);
            xi = Q_y * diag(sign(sign(diag(R_y))+.5));  % qr-based retraction
            R_new(:,:,i) = xi;
        end
        
        % compute the distance to groundtruch R_orig
        stepStartTime = tic;
        if mod(iter, ReSync_parameters.check_freq) == 0
            mse = cmpt_mse(R_est, R_orig);
%             mse = cmpt_mse(permute(R_est,[2,1,3]), R_orig);
            fprintf('estimation error: %d in iteration %d \n',mse, iter); 
            Dist = [Dist, mse];
            
%             REiter = norm(R_new(:) - R_est(:))/norm(R_new(:));
%             if REiter < ReSync_parameters.stop_threshold
%                 break;
%             end
            if isfield(ReSync_parameters,'target')
                if mse < ReSync_parameters.target, break; end
            else
                if step < ReSync_parameters.stop_threshold, break; end
            end
%             if step < ReSync_parameters.stop_threshold
%                 break;
%             end
        end
        stepTime = stepTime + toc(stepStartTime);
        
        R_est = R_new;

    end
% totalTime = toc(totalStartTime);
% finalTime=totalTime-stepTime;
end