function [C, ref_rot, common_lines_matrix] = gen_and_save_cl(K, SNR)
    load ('Data/cleanrib.mat'); 
    
    k = size(volref,1);  % data set size
    n = 129; % 65;129
    V = my_NewSizeVol(volref,n);   % upsampling
    volref = V;

    ref_rot = rand_rots(K);
    A = OpNufft3D(ref_rot,n);

    % ground-truth projection
    projs = A * volref;
    % figure; viewstack(projs,5,5);

    % add noise to projections according to SNR
    [noisy_projs, ~] = ProjAddNoise(projs, SNR);
    % figure;viewstack(noisy_projs,5,5);

    % Mask the projections with a radius mask of size 'rad'
    masked_r = 45; % mask radius
    masked_projs=mask_fuzzy(noisy_projs,masked_r);
    % figure;viewstack(masked_projs,5,5);

    %% 3. compute polar Fourier Transform of projections
    n_theta = 360; %360;%72
    n_r = 100;     %100;%33
    [npf,~]=cryo_pft(masked_projs,n_r,n_theta,'single');

    %% 4. find common lines from Fourier projections in 3D
    max_shift=0;
    shift_step=1;
    common_lines_matrix = commonlines_gaussian(npf,max_shift,shift_step);
    C = clstack2C( common_lines_matrix,n_theta ); % common lines matrix
    % Find reference common lines and compare
    [ref_clstack,~]=gen_clmatrix(ref_rot,n_theta);
    p = comparecl( common_lines_matrix, ref_clstack, n_theta, 10 );
    %% 5. save common line matrix
    noisy_cl_save_dir = ['Data/cl_matrix/SNR_' num2str(1/SNR) 'K_' num2str(K) '_noisy_cl.mat'];
    noisy_cl_index_save_dir = ['Data/cl_matrix/SNR_' num2str(1/SNR) 'K_' num2str(K) '_noisy_cl_index.mat'];
    ref_rot_save_dir = ['Data/cl_matrix/SNR_' num2str(1/SNR) 'K_' num2str(K) '_ref_rot.mat'];
    save(noisy_cl_save_dir, 'C');
    save(noisy_cl_index_save_dir, 'common_lines_matrix');
    save(ref_rot_save_dir, 'ref_rot');
    
end

