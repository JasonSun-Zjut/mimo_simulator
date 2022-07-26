function [ nRx_tss_matrix_with_noise ] = CORE_add_noise( nRx_tss_matrix_wout_noise, snr_db )

sigma_noise_matrix  = nRx_tss_matrix_wout_noise / (10^(snr_db/20)) ;
noise               = normrnd(0, sigma_noise_matrix);

nRx_tss_matrix_with_noise =  nRx_tss_matrix_wout_noise + noise;

% CORRECT NEGATIVES
negative_mask = (nRx_tss_matrix_with_noise < 0);
nRx_tss_matrix_with_noise(negative_mask) = 0;

end

