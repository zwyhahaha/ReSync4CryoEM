function p = get_p_from_snr(SNR)
if SNR==1
    p=0.966;
elseif SNR==1/2
    p=0.919;
elseif SNR==1/4
    p=0.813;
elseif SNR==1/8
    p=0.638;
elseif SNR==1/16
    p=0.437;
elseif SNR==1/32
    p=0.252;
elseif SNR==1/64
    p=0.115;
elseif SNR==1/128
    p=0.047;
elseif SNR==1/256
    p=0.023;
else
    p=0.01;
end
end

