function [ res ] = EXAMPLE_runner_02_Ntx_vs_BER(  )

% Set Parameters
fprintf(1, '\n################# Ntx vs BER (Mod:BCSK;BMoSK) Example Parameters  ############');

delta_t              = 0.015;
fprintf(1,'\n ## delta_t = %f s', delta_t);

dist_inMicroMeters   = 2;
fprintf(1,'\n ## Distance = %d micrometer', dist_inMicroMeters);


nsym                 = 2001;
fprintf(1,'\n ## Num. Symbols = %d ', nsym);

replication          = 100;
fprintf(1,'\n ## Replication = %d ', replication);

ts                   = 0.030;
tss                  = delta_t;
fprintf(1,'\n ## (ts, tss) = (%f, %f) ', ts, tss);

fprintf(1, '\n################# Ntx vs BER (Mod:BCSK;BMoSK) Example Parameters  ############');
fprintf(1, '\n############### If MFSK is used arrange tss and delta_t accordingly  ############');

num_molecules_to_emit = 200;
ntx_list              = 200:100:600; 

% Prepare Variables
[tx_node, rx_node, env_params, sim_params] = prepare_vars4runners_PointSrc(dist_inMicroMeters, [0 0 0], 4, 100, delta_t, num_molecules_to_emit, ts, tss, [0.5 0.5], nsym, replication);

%%%%%% GENERATE SYMBOL SEQUENCES
tx_sym_matrix       = CORE_prepare_symbol_sqces(sim_params.symbol_probs, sim_params.replication, sim_params.nsym);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Run Sims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BCSK (pulse)
fprintf(1, '\n--------------    Running  BCSK     ------------');
tx_node.mod                     = 0; %% BCSK (pulse)
rx_node.demod                   = tx_node.mod;

res.BCSK = runner_Ntx_vs_BER(tx_sym_matrix, tx_node, env_params, sim_params, rx_node, ntx_list);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BMoSK (pulse 2typeMol)
fprintf(1, '\n--------------    Running  BMoSK     ------------');
tx_node.mod                     = 3; %% BMoSK (pulse 2typeMol)
rx_node.demod                   = tx_node.mod;

res.BMoSK = runner_Ntx_vs_BER(tx_sym_matrix, tx_node, env_params, sim_params, rx_node, ntx_list);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BMFSK (cosine)
% fprintf(1, '\n--------------    Running  BMFSK     ------------');
% tx_node.mod                     = 2; %% BMFSK (cosine)
% rx_node.demod                   = tx_node.mod;
% 
% res.BMFSK = runner_Ntx_vs_BER(tx_sym_matrix, tx_node, env_params, sim_params, rx_node, ntx_list);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ber_bcsk_withN    = res.BCSK.NTX_vs_BER_row4ts_withN;
ber_bmosk_withN   = res.BMoSK.NTX_vs_BER_row4ts_withN;


figure;
semilogy(res.BCSK.ntx_list, ber_bcsk_withN(1,:), '--',  'LineWidth', 3, 'Color', [0.9 0.1 0.2]);
hold on
semilogy(res.BMoSK.ntx_list, ber_bmosk_withN(1,:), '-.', 'LineWidth', 3, 'Color', [0.4 0.4 1]);
ylabel('BER')
xlabel('SNR (dB)')
grid on
legend('BCSK', 'BMoSK', 3)
hold off
end



function [res] = runner_Ntx_vs_BER(tx_sym_matrix, tx_node, env_params, sim_params, rx_node, ntx_list)

ntx_list_size     = size(ntx_list, 2); 
replication       = sim_params.replication;

res.ntx_list                        = ntx_list;
res.ntx_list_size                   = ntx_list_size;

res.NTX_vs_BER_row4ts_withN         = ones(1, ntx_list_size);
res.t_per_replication               = zeros(1, ntx_list_size);
res.tx_sym_matrix                   = tx_sym_matrix;

for jj=1:ntx_list_size
   sim_params.molecules_perTs      = ntx_list(jj);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% ADD Header and Tail (can be outside the for loop, it is invariant)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   tx_sym_matrix_HT = CORE_add_HT(ones(replication, 1), tx_sym_matrix, tx_sym_matrix(:, end));
   % Added "1"        to the HEAD for warm up and time synchronization offset
   % Added "1" or "0" to the TAIL, last symbol twice for synchronization shift
   H_cnt = 1;
   T_cnt = 1;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% Modulate + SimulateDiffusion for each row
   %%%
   %%% Dimensions of nRx_raw_matrix_wout_noise: (nRx_row_cnt_molType x nsym*ts_step x replication)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [ nRx_raw_matrix_wout_noise, stats ] = CORE_sim_replicator( ...
      'CORE_sim_diffusion_3d_P2S_FAST', ...% Simulator Name ##You Can Change This##
      tx_sym_matrix_HT, ... % Tx Symbol Sequences + HT for each replication
      tx_node,    ...% Tx node properties
      rx_node,    ...% Rx node properties
      env_params, ...% Environment properties
      sim_params );
   
   stats.tElapsed_all
   res.t_per_replication(jj) = stats.tElapsed_all / replication;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%  * PREPARE Ntx Variations for WithNoise and DEMOD
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if ((rx_node.demod == 0) || (rx_node.demod == 1))
      % Find effective threshold
      fprintf('\nEstimating Best THR with partial data');
      rx_node.threshold = find_best_thr(nRx_raw_matrix_wout_noise(:,:,1:5), tx_sym_matrix_HT(1:5,:), H_cnt, T_cnt, rx_node, sim_params);
      fprintf('\nEstimating Best THR with partial data [DONE.best_thr = %d]', rx_node.threshold);
      
      res.NTX_vs_BER_row4ts_withN(jj) = prepare_NTX_vs_BER(nRx_raw_matrix_wout_noise, tx_sym_matrix, H_cnt, T_cnt, rx_node, sim_params);
   elseif ((rx_node.demod == 2) || (rx_node.demod == 2.1))
      % Find effective offset
      fprintf('\nEstimating Best OFFSET with partial data');
      rx_node.synch_offset = find_best_synch_offset(nRx_raw_matrix_wout_noise(:,:,1:5), tx_sym_matrix_HT(1:5,:), H_cnt, T_cnt, rx_node, sim_params);
      fprintf('\nEstimating Best OFFSET with partial data [DONE.best_synch_offset = %d]', rx_node.synch_offset);
      
      res.NTX_vs_BER_row4ts_withN(jj) = prepare_NTX_vs_BER(nRx_raw_matrix_wout_noise, tx_sym_matrix, H_cnt, T_cnt, rx_node, sim_params);
   elseif ((rx_node.demod == 3) || (rx_node.demod == 4))
      % No THR no OFFSET
      res.NTX_vs_BER_row4ts_withN(jj) = prepare_NTX_vs_BER(nRx_raw_matrix_wout_noise, tx_sym_matrix, H_cnt, T_cnt, rx_node, sim_params);
   else
      % TODO checks
   end
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%     UTILITY FUNCTIONS     %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res] = prepare_NTX_vs_BER(nRx_raw_matrix_wout_noise, tx_sym_matrix, H_cnt, T_cnt, rx_node, sim_params)

sym_cnt             = size(tx_sym_matrix, 2);
replication         = size(tx_sym_matrix, 1);

demod_sym_matrix    = zeros(replication, sym_cnt);

% SHIFT, REMOVE HT, and SAMPLE  + Add NOISE + DEMOD each row
    for jj=1:replication
       nRx_singleRun_wout_noise_tss = CORE_shift_removeHT_and_sample( nRx_raw_matrix_wout_noise(:,:,jj), H_cnt, T_cnt, rx_node, sim_params );
       demod_sym_matrix(jj, :) = CORE_block_demodulate(nRx_singleRun_wout_noise_tss, rx_node, sim_params);
    end

res = nnz(demod_sym_matrix-tx_sym_matrix) / numel(tx_sym_matrix);

end


function [best_thr] = find_best_thr(nRx_raw_matrix_wout_noise_partial, tx_sym_matrix_HT_partial, H_cnt, T_cnt, rx_node, sim_params)
%%% Dimensions of                   nRx_raw_matrix_wout_noise: (nRx_row_cnt_molType x nsym*ts_step x replication)

ts                   = sim_params.ts_inSeconds;
tss                  = sim_params.tss_inSeconds;
rcv_in_ts            = round(ts/tss);

row_sz               = size(tx_sym_matrix_HT_partial, 1); % This is replication count for partial data

for ii=1:row_sz
    % Shift, crop HT, and sample
    nRx_matrix_wout_noise_Shifted_RemHT_Sampled(1,:,ii) = CORE_shift_removeHT_and_sample(nRx_raw_matrix_wout_noise_partial(1,:,ii), H_cnt, T_cnt, rx_node, sim_params);
end
tx_sym_matrix_RemHT_part                    = tx_sym_matrix_HT_partial(:, 1+H_cnt:end-T_cnt);

sym_cnt              = size(nRx_matrix_wout_noise_Shifted_RemHT_Sampled, 2)/rcv_in_ts;
replication          = size(nRx_matrix_wout_noise_Shifted_RemHT_Sampled, 3);

tstats_in_ts         = sum(reshape(nRx_matrix_wout_noise_Shifted_RemHT_Sampled(1,:,1), rcv_in_ts, sym_cnt));
threshold_min        = round(min(tstats_in_ts));
threshold_max        = round(max(tstats_in_ts));

thr_cnt              = 30;
delta_increment      = (threshold_max - threshold_min) / thr_cnt;

thresholds = threshold_min:delta_increment:threshold_max;

thresholds_col_cnt = size(thresholds,2);

ber_vs_thr = ones(1, thresholds_col_cnt);
for ii= 1:thresholds_col_cnt
    
    rx_node.threshold = thresholds(ii);
    
    % DEMOD each row and count errors
    err_cnt = 0;
    for jj=1:replication
        demod_syms = CORE_block_demodulate(nRx_matrix_wout_noise_Shifted_RemHT_Sampled(:,:,jj), rx_node, sim_params);
        err_cnt = err_cnt + nnz(demod_syms - tx_sym_matrix_RemHT_part(jj,:));
    end
    
    ber_vs_thr(ii) = err_cnt/numel(tx_sym_matrix_RemHT_part);
end

left_idx    = 1;
right_idx   = thresholds_col_cnt;
left_best   = 1;
right_best  = 1;
for ii=1:thresholds_col_cnt
    curr_left_ber = ber_vs_thr(ii);
    curr_right_ber = ber_vs_thr(thresholds_col_cnt-ii+1);
    if (curr_left_ber < left_best)
        left_idx = ii;
        left_best = curr_left_ber;
    end
    if (curr_right_ber < right_best)
        right_idx = thresholds_col_cnt-ii+1;
        right_best = curr_right_ber;
    end
end
avg_index = (2*left_idx+right_idx)/3;

best_thr = threshold_min + (avg_index-1)*delta_increment;

        
end



function [res_best_offset] = find_best_synch_offset(nRx_raw_matrix_wout_noise_partial, tx_sym_matrix_HT_partial, H_cnt, T_cnt, rx_node, sim_params)
%%% Dimensions of                                   nRx_raw_matrix_wout_noise: (nRx_row_cnt_molType x nsym*ts_step x replication)

ts                      = sim_params.ts_inSeconds;
delta_t                 = sim_params.delta_t;

half_ts_cnt             = round(0.5 * (ts / delta_t));

row_sz                  = size(tx_sym_matrix_HT_partial, 1);

best_offset_list        = zeros(row_sz, 1);
for ii=1:row_sz
    % Crop TxSymMatrix
    tx_sym_rowii_RemHT          = tx_sym_matrix_HT_partial(ii, 1+H_cnt:end-T_cnt);
    
    best_berr = 1;
    best_offset = 0;
    for jj=0:half_ts_cnt
        rx_node.synch_offset    = jj;
        nRx_matrix_wout_noise_Shifted_RemHT_Sampled = CORE_shift_removeHT_and_sample(nRx_raw_matrix_wout_noise_partial(:,:,ii), H_cnt, T_cnt, rx_node, sim_params);
        
        demod_syms  = CORE_block_demodulate(nRx_matrix_wout_noise_Shifted_RemHT_Sampled, rx_node, sim_params);
        
        curr_berr   = nnz(tx_sym_rowii_RemHT-demod_syms)/numel(tx_sym_rowii_RemHT);
        
        if (curr_berr < best_berr)
            best_berr = curr_berr;
            best_offset = jj;
        end
    end
    % Now we have best offset for row ii
    best_offset_list(ii) = best_offset;
end
res_best_offset = round(mean(best_offset_list)+0.2);

end
