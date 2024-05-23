function mse = cmpt_mse(rotations,rots_ref)
    K = size(rotations, 3);
    TOL=1.0e-12;
    J=[1 0 0; 0 1 0; 0 0 -1]; % Reflection matrix
    
    % Make sure that we got rotations.
    for k=1:K
        R=rotations(:,:,k);
        erro=norm(R*R.'-eye(3));
        if erro>TOL
            fprintf('Trnaformation %d is not orthogonal, err=%e  tol=%e\n',k,erro,TOL);
        end

        errd=abs(det(R)-1);
        if errd>TOL
            fprintf('Determinant of %d diffrs from 1, err=%e  tol=%e\n',k,errd,TOL);    
        end

        % Enforce R to be a rotation (in case the error is large)
        [U,~,V]=svd(R);
        rotations(:,:,k)=U*V.';
    end
    
    for k1=1:K-1
        for k2=k1+1:K
            R1=rotations(:,:,k1);
            R2=rotations(:,:,k2);
            R=R1.'*R2;
            
            R1ref=rots_ref(:,:,k1);
            R2ref=rots_ref(:,:,k2);
            inv_R1ref=R1ref.';
            inv_R2ref=R2ref.';
            
            % The resulting rotations should satisfy the same ratio
            % equaions as the true orientations. Specifically, the ration
            % of each pair of estimated rotation should equal one of the
            % following two rotations:
            Rref1=inv_R1ref.'*inv_R2ref;
            Rref2=J*Rref1*J;
            
            err1=norm(Rref1-R,'fro')/norm(Rref1,'fro');
            err2=norm(Rref2-R,'fro')/norm(Rref2,'fro');
            
        end
    end
    
    % Register estimated rotations to true ones, and compute the difference
    % between the two.
    rot=zeros(3*K,3);  % The K estimated rotation matrices stacked as a matrix with 3K rows.
    rot1=zeros(3*K,3); % True true K rotation matrices stacked as a matrix with 3K rows.
    rot2=zeros(3*K,3); % Reflected matrices of rot1, which are also a valid solution for the rotations.
    
    for k=1:K
        R=rotations(:,:,k);
        rot(3*(k-1)+1:3*k,:)=R.';
        Rref=rots_ref(:,:,k);
        rot1(3*(k-1)+1:3*k,:)=Rref;
        rot2(3*(k-1)+1:3*k,:)=J*Rref*J;
    end
    
    % Compute the two possible orthogonal matrices which register the
    % estimated rotations to the true ones.
    O1=rot.'*rot1./K;       
    O2=rot.'*rot2./K;
    
    % In cany case, enforce the registering matrix O to be a rotation.
    if err1<err2
        [U,~,V]=svd(O1); % Use O1 as the registering matrix
        flag=1;
    else
        [U,~,V]=svd(O2); % Use O2 as the registering matrix
        flag=2;
    end
    O=U*V.';    
    
    % Plot estimation errors
    diff=zeros(K,1);
    mse=0;
    for k=1:K
        R=rotations(:,:,k);
        Rref=rots_ref(:,:,k);
        if flag==2
            Rref=J*Rref*J;
        end
        diff(k)=norm(R.'*O-Rref,'fro');
        mse=mse+diff(k).^2;
    end
    mse=mse/K;
%    hist(diff)
 
end

