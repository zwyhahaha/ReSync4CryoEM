function[model_out]=Uniform_Topology(n,p,n_theta)
    % generate rotation matrices
    R_orig = zeros(3,3,n);
    quat = randn(n, 4);
    quat = normalize(quat, 2, 'norm');
    for k=1:n
        if (quat(k,1) < 0)
            quat(k,:) = -quat(k,:);
        end
    end

    for i = 1:n
        q = quat(i,:);
        R_orig(:,:,i)=my_quat2rotm(q);
    end
    
    % generate common lines
%     Cij_orig = zeros(3,n,n);
%     for i = 1:n
%         for j = 1:i-1
%             crs = cross(R_orig(:,3,i),R_orig(:,3,j))/norm(cross(R_orig(:,3,i),R_orig(:,3,j)));
%             Cij_orig(:,i,j) = R_orig(:,:,i)'*crs;
%             Cij_orig(:,j,i) = R_orig(:,:,j)'*crs;
%         end
%     end
    
    % generate corrupted common lines
%     G = rand(n,n) >= p;
%     G = tril(G,-1);
%     [Ind_j, Ind_i] = find(G==1);
%     m = length(Ind_i);
%     
%     C_corr = rand(2,m);
%     C_corr = normalize(C_corr,1,'norm');
%     C_corr = vertcat(C_corr, zeros(1,m));
%     Cij = Cij_orig;
%     for k = 1:m
%         i=Ind_i(k); j=Ind_j(k); 
%         Cij(:,i,j)=C_corr(:,k);
%     end
    
    clmatrix=gen_clmatrix(R_orig,n_theta);  % generate common line index
    [noisy_cl,~]=corr_clmatrix(clmatrix,n_theta,p);  % perturb common line index
    Cij_orig = gen_Cij(clmatrix, n_theta);  % get common line coordinate
    Cij = gen_Cij(noisy_cl, n_theta);   % get perturbed common line coordinate
%     S=syncmatrix_vote(noisy_cl,n_theta);
    
    model_out.R_orig = R_orig;
    model_out.Cij_orig = Cij_orig;
    model_out.Cij = Cij;
    model_out.noisy_cl = noisy_cl;
%     model_out.S = S;
end

