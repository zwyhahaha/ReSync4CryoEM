function [R_est, Dist, stepTime, iter] = ReSync_norm(Cij, R_init, R_orig, ReSync_parameters)
fprintf('Starting ReSync_norm... \n');
    R_est = R_init;
    [~, ~, n] = size(R_init);
    scalar = sqrt(2*n/3);
    Dist = [];
    stepTime = 0;
    for iter = 1 : ReSync_parameters.max_iter
        % linear decay of stepsize
        step = ReSync_parameters.stepsize * ReSync_parameters.decay^(iter);
        
        % compute the subgradient
        [~,~, ~, x_egrad] = FunCostAndGradp2q1_ReSync(R_est(:,1:3,:),Cij(1:2,:,:));
        x_egrad = x_egrad*2;
        
        R_new = R_est - step * x_egrad;
        
        flag = 1;
        R_neww = zeros(2*n, 3);
        for i = 1 : n
            R_neww(2*i-1:2*i,:) = R_new(1:2,:,i);
        end
        R_tmp = R_neww;
        
        h = 0;
        while flag
            h = h + 1;
            for i = 1 : n
                y = R_neww(2*i-1:2*i,:)';
                [Q_y, R_y] = qr(y);
                R_y(1:2, 1:2) = diag(sign(sign(diag(R_y))+.5));
                xi = Q_y * R_y;
                R_tmp(2*i-1:2*i,:) = xi';
            end
            
            re = norm(R_neww-R_tmp, 'fro');
            if re < 1e-5
                flag=0;
            end
            
            [Q_y, R_y] = qr(R_tmp);
            R_y(1:3, 1:3) = diag(sign(sign(diag(R_y))+.5));
            R_tmp = Q_y * R_y;
            R_tmp = R_tmp * scalar;
            
            R_neww = R_tmp;
        end
%         disp(h);
        
        for i = 1 : n
            y = R_neww(2*i-1:2*i,:)';
            [Q_y, R_y] = qr(y);
            R_y(1:2, 1:2) = diag(sign(sign(diag(R_y))+.5));
            xi = Q_y * R_y;
            R_new(1:2,:,i) = xi';
            R_new(3,:,i) = cross(xi(:,1)',xi(:,2)');
        end
        
        % compute the distance to groundtruch R_orig
        stepStartTime = tic;
        if mod(iter, ReSync_parameters.check_freq) == 0
            mse = cmpt_mse(R_est, R_orig);
            fprintf('estimation error: %d in iteration %d \n',mse, iter);
            Dist = [Dist, mse];
            if isfield(ReSync_parameters,'target')
                if mse < ReSync_parameters.target, break; end
            else
                if step < ReSync_parameters.stop_threshold, break; end
            end
            
%             REiter = norm(R_new(:) - R_est(:))/norm(R_new(:));
%             if REiter < stop_threshold
%                 break;
%             end
        end
        stepTime = stepTime + toc(stepStartTime);
        R_est = R_new;

    end
end

