function vol = my_NewSizeVol(volume,N)

nv=size(volume,1);
% If the volume is too small for the given size, upsample it as needed.
if N>nv+1 % For compaibility with gen_projections, allow one pixel aliasing.
    % More precisely, it should be N>nv, however, by using nv+1 the
    % results match those of gen_projections.
    if mod(N-nv,2)==1
        error('Upsampling from odd to even sizes or vice versa is currently not supported');
    end
    dN=floor((N-nv)/2);
    fv=cfftn(volume);
    padded_volume=zeros(N,N,N);
    padded_volume(dN+1:dN+nv,dN+1:dN+nv,dN+1:dN+nv)=fv;
    volume=icfftn(padded_volume);
    assert(norm(imag(volume(:)))/norm(volume(:))<1.0e-5);
    nv=N; % The new volume size
end
vol = volume;
end

