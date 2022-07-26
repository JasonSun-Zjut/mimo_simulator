function [ res ] = CORE_shift_removeHT_and_sample( nRx_singleRun_wout_noise, H_cnt, T_cnt, rx_node, sim_params )
% Dimensions:
% nRx_singleRun_wout_noise (nRx_row_cnt_molType x nsym*ts_step)

ts                      = sim_params.ts_inSeconds;
tss                     = sim_params.tss_inSeconds;
delta_t                 = sim_params.delta_t;

ts_step                 = round(ts / delta_t);
tss_step                = round(tss / delta_t);
sampling_cnt_for_symbol = round(ts / tss);


synch_offset            = rx_node.synch_offset;

% Synchronization and Removing First H_cnt and Last T_cnt Syms
t_start = 1 + synch_offset + (H_cnt * ts_step); % remove first symbol and synch
t_end   = size(nRx_singleRun_wout_noise,2) + synch_offset - (T_cnt*ts_step); % remove last symbol and synch

nsym = round( (t_end - t_start + 1) / ts_step );

if ((rx_node.demod == 3) || (rx_node.demod == 4))
    if(tss_step == 1)
       res = nRx_singleRun_wout_noise(:,t_start:t_end);
    else
       nRx1_wout_noise_tss = sum( reshape(nRx_singleRun_wout_noise(1,t_start:t_end), tss_step, nsym*sampling_cnt_for_symbol) );

       nRx2_wout_noise_tss = sum( reshape(nRx_singleRun_wout_noise(2,t_start:t_end), tss_step, nsym*sampling_cnt_for_symbol) );
       
       res = [nRx1_wout_noise_tss; nRx2_wout_noise_tss];
    end
else
    if(tss_step == 1)
       res = nRx_singleRun_wout_noise(t_start:t_end);
    else
       res = sum( reshape(nRx_singleRun_wout_noise(t_start:t_end), tss_step, nsym*sampling_cnt_for_symbol) );
    end
end

end

