function  C  = clstack2C( clstack,L )
% Store the common line vectors c_ij , i,j=1,...K from the common line
% stack. C(:,i,j)=(c_ij(1), c_ij(2))';
% 
% Lanhui Wang
% Jul 18, 2012
K=size(clstack,1);
% L=max(clstack(:));
fprintf('K=%d, L=%d\n',K,L);

C=zeros(2,K,K);
for i=1:K
    j=(i+1):K;
    l1 = clstack(i,j)-1;
    l2 = clstack(j,i)-1;
    l1=l1(:);
    l2=l2(:);
    x12 = cos(2*pi*l1/L);
    y12 = sin(2*pi*l1/L);
    x21 = cos(2*pi*l2/L);
    y21 = sin(2*pi*l2/L);
    C(1,i,j)=x12;
    C(2,i,j)=y12;
    C(1,j,i)=x21;
    C(2,j,i)=y21;
end

end

