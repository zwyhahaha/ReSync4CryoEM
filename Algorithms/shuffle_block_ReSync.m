function [R_est, Dist, stepTime, iter] = shuffle_block_ReSync(Cij, R_init, R_orig, ReSync_parameters)
%     totalStartTime = tic;
%     stop_threshold = 1e-5;
    R_est = R_init;
    [~, ~, n] = size(R_init);
    Dist = [];
    fprintf('Starting random_ReSync...filter ratio %0.2f \n', ReSync_parameters.filter_size/n);
    num_blocks = ceil(n/ReSync_parameters.filter_size);
%     ReSync_parameters.check_freq = ReSync_parameters.check_freq/ReSync_parameters.filter_ratio;
    stepTime = 0;
    for iter = 1 : ReSync_parameters.max_iter
        indices = randperm(n);
        for idx = 1:num_blocks
%             filter = indices(1:ReSync_parameters.filter_size);
            start_index = (idx - 1) * ReSync_parameters.filter_size + 1;
            end_index = min(idx * ReSync_parameters.filter_size, n);
            filter = indices(start_index:end_index);
            if ReSync_parameters.use_homofilter
                filter2 = filter;
            else
                indices = randperm(n);
                filter2 = indices(1:ReSync_parameters.filter_size);
            end

            % linear decay of stepsize
            step = ReSync_parameters.stepsize * ReSync_parameters.decay^(iter);

            % compute the subgradient 
            [~,~, x_egrad] = FunCostAndGradp2q1_random_ReSync2(R_est(:,:,filter),R_est(:,:,filter2),Cij(1:2,filter,filter2),Cij(1:2,filter2,filter));
            x_egrad = x_egrad*2;

            R_new = R_est;
            R_new(:,:,filter) = R_new(:,:,filter) - step * x_egrad;
            for i = filter
                y = R_new(:,:,i);
                [Q_y, R_y] = qr(y);
                xi = Q_y * diag(sign(sign(diag(R_y))+.5));  % qr-based retraction
                R_new(:,:,i) = xi;
            end

           
            R_est = R_new;

        end
        % compute the distance to groundtruch R_orig
        if mod(iter, ReSync_parameters.check_freq) == 0
            stepStartTime = tic;
            mse = cmpt_mse(R_est, R_orig);
            stepTime = stepTime + toc(stepStartTime);

            fprintf('estimation error: %d in iteration %d \n',mse, iter); 
            Dist = [Dist, mse];

            if isfield(ReSync_parameters,'target')
                if mse < ReSync_parameters.target, break; end
            else
                if step < ReSync_parameters.stop_threshold, break; end
            end
%             if step < ReSync_parameters.stop_threshold
%                 break;
%             end
        end
    end
end

%% The following code is intended for parallelization
% function [R_est, Dist, stepTime, iter] = shuffle_block_ReSync(Cij, R_init, R_orig, ReSync_parameters)
% %     totalStartTime = tic;
% %     stop_threshold = 1e-5;
%     R_est = R_init;
%     [~, ~, n] = size(R_init);
%     Dist = [];
%     fprintf('Starting random_ReSync...filter ratio %0.2f \n', ReSync_parameters.filter_size/n);
%     num_blocks = ceil(n/ReSync_parameters.filter_size);
%     stepTime = 0;
%     for iter = 1 : ReSync_parameters.max_iter
%         indices = randperm(n);
% %         R_new = R_est(indices,indices);
% %         Cij_copy = cell(num_blocks);
%         R_new = cell(num_blocks);
%         filters = cell(num_blocks);
%         for idx = 1:num_blocks
%             start_index = (idx - 1) * ReSync_parameters.filter_size + 1;
%             end_index = min(idx * ReSync_parameters.filter_size, n);
%             filter = indices(start_index:end_index);
% %             Cij_copy{idx} = Cij(1:2,filter,filter);
%             R_new{idx} = R_est(:,:,filter);
%             filters{idx} = filter;
%         end
%         % linear decay of stepsize
%         step = ReSync_parameters.stepsize * ReSync_parameters.decay^(iter);
%         for idx = 1:num_blocks
%             filter=filters{idx};
%             % compute the subgradient 
%             [~,~, x_egrad] = FunCostAndGradp2q1_random_ReSync2(R_est(:,:,filter),R_est(:,:,filter),Cij(1:2,filter,filter),Cij(1:2,filter,filter));
%             x_egrad = x_egrad*2;
%             
%             R_new{idx} = R_est(:,:,filter) - step * x_egrad;
%             for i = 1:length(filter)
%                 y = R_new{idx}(:,:,i);
%                 [Q_y, R_y] = qr(y);
%                 xi = Q_y * diag(sign(sign(diag(R_y))+.5));  % qr-based retraction
%                 R_new{idx}(:,:,i) = xi;
%             end
%         end
% %             start_index = (idx - 1) * ReSync_parameters.filter_size + 1;
% %             end_index = min(idx * ReSync_parameters.filter_size, n);
% %             filter = indices(start_index:end_index);
%         for idx = 1:num_blocks
%             R_est(:,:,filters{idx})=R_new{idx};
%         end
% %         R_est = R_new;
%         % compute the distance to groundtruch R_orig
%         if mod(iter, ReSync_parameters.check_freq/ReSync_parameters.filter_ratio) == 0
%             stepStartTime = tic;
%             mse = cmpt_mse(R_est, R_orig);
%             stepTime = stepTime + toc(stepStartTime);
% 
%             fprintf('estimation error: %d in iteration %d \n',mse, iter); 
%             Dist = [Dist, mse];
% 
%             if step < ReSync_parameters.stop_threshold
%                 break;
%             end
%         end
%     end
% end