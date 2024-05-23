function Cij = gen_Cij(clmatrix, n_theta)
% generate common line coordinate from common line index
N = size(clmatrix,2);
Cij = zeros(3, N, N);
for k1=1:N-1    
    for k2=k1+1:N
        l_k1k2 = clmatrix(k1,k2)-1;
        l_k2k1 = clmatrix(k2,k1)-1;
        PI=4*atan(1.0);
        c_k1k2 = [cos(2*PI*l_k1k2/n_theta); sin(2*PI*l_k1k2/n_theta); 0];
        c_k2k1 = [cos(2*PI*l_k2k1/n_theta); sin(2*PI*l_k2k1/n_theta); 0];
        Cij(:,k1,k2) = c_k1k2;
        Cij(:,k2,k1) = c_k2k1;
    end
end
end

