function [R_est,time] = Eigenvector_Relaxation(Cij)
tic;
    [~,n,~] = size(Cij);
    d = 3;
    Xij = squeeze(Cij(1,:,:));
    Yij = squeeze(Cij(2,:,:));
    S = [Xij.* Xij', Xij.* Yij'; Yij.* Xij', Yij.* Yij'];
    for i = 1 : n
        S(i,i) = 0;
        S(i, i+n) = 0;
        S(i+n, i) = 0;
        S(i+n, i+n) = 0;
    end
    [V,~] = eigs(S, d);
%     V = V * diag([1,1,-1]);
    R_est = zeros(d,d,n);
    normR = 0;
    Q_est = zeros(d,d,n);
    normQ = 0;
    
    for i = 1:n
        A1 = V(i,:)' * sqrt(2* n / 3);
        A2 = V(n+i, :)' * sqrt(2* n / 3);
        A3 = cross(A1, A2);
        Ri = [A1, A2, A3];
        [Ur,Lr,Vr] = svd(Ri);
        S0 = diag([ones(1,d-1),1]);
        normR = norm(S0 - Lr,"fro")^2; % the distance dist(V, SO(d)^n)
        R_est(:,:,i) = Ur*S0*Vr';  % the projection on SO(d)
        
%         A1 = diag([1,1,-1])*V(i,:)' * sqrt(2* n / 3);
%         A2 = diag([1,1,-1])*V(n+i, :)' * sqrt(2* n / 3);
%         A3 = cross(A1, A2);
%         Qi = [A1, A2, A3];
%         Qi = Ri * diag([1,1,-1]);
        Qi = diag([1,1,-1]) * Ri;
%         Qi = [Ri(1,:); Ri(2,:); -Ri(3,:)];
        [Ur,Lr,Vr] = svd(Qi);
        S0 = diag([ones(1,d-1),det(Ur*Vr')]);
        normQ = norm(S0 - Lr,"fro")^2;  % the distance dist(U, SO(d)^n)
        Q_est(:,:,i) = Ur*S0*Vr';  % the projection on SO(d)
    end
    
    if normR > normQ 
        R_est = Q_est;
    end
    time = toc;
    fprintf('Initialization completed! \n \n'); 
end

