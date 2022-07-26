function [ res ] = EXAMPLE_runner_01_Thr_vs_BER(  )

fprintf(1, '\n################# Threshold vs BER (Mod:BCSK) Example Parameters  ############');

delta_t              = 0.04;
fprintf(1,'\n ## delta_t = %f s', delta_t);

dist_inMicroMeters   = 2;
fprintf(1,'\n ## Distance = %d micrometer', dist_inMicroMeters);

num_molecules_to_emit = 500;
fprintf(1,'\n ## Num. Emitted Molecules = %d ', num_molecules_to_emit);

nsym                 = 500;
fprintf(1,'\n ## Num. Symbols = %d ', nsym);

replication          = 111;
fprintf(1,'\n ## Replication = %d ', replication);
%symbol duration%
ts                   = 0.08;
%sampling duration%
tss                  = delta_t;
fprintf(1,'\n ## (ts, tss) = (%f, %f) ', ts, tss);
%Binary  Concentration Shift Keying%
fprintf(1, '\n################# Threshold vs BER (Mod:BCSK) Example Parameters  ############');


% Prepare Variables
[tx_node, rx_node, env_params, sim_params] = prepare_vars4runners_PointSrc(dist_inMicroMeters, [0 0 0], 4, 100, delta_t, num_molecules_to_emit, ts, tss, [0.5 0.5], nsym, replication);

%%%%%% GENERATE SYMBOL SEQUENCES
%%%%%%111*2001 symbols
tx_sym_matrix       = CORE_prepare_symbol_sqces(sim_params.symbol_probs, sim_params.replication, sim_params.nsym);

% Run 
res = runner_01_Thr_vs_BER( tx_sym_matrix, tx_node, rx_node, env_params, sim_params);


% PLOT
figure;
semilogy(res.thr_vec(1:500), res.ber_vec(1:500), '--r', 'LineWidth', 2)
hold on
grid on

set(gca, 'fontsize',11, 'fontweight','d');
xlabel('Threshold')
ylabel('BER')

end


function [res] = runner_01_Thr_vs_BER(tx_sym_matrix, tx_node, rx_node, env_params, sim_params)
tx_node.mod                = 0; %% BCSK (pulse)
rx_node.demod              = tx_node.mod;

% Run Sim
[ nRx_raw_matrix_wout_noise, stats ] = CORE_sim_replicator( ...
           'CORE_sim_diffusion_3d_P2S_FAST', ...% Simulator Name 
           tx_sym_matrix, ... % Tx Symbol Sequences
           tx_node,    ...% Tx node properties
           rx_node,    ...% Rx node properties
           env_params, ...% Environment properties
           sim_params );

res.nRx_raw_matrix_wout_noise = nRx_raw_matrix_wout_noise;
res.stats = stats;

% DeMod
fprintf(1,'\n\n Demodulation Process Starts');
max_thr = sim_params.molecules_perTs;

ber_vec = ones(1, max_thr+1); % one more threshold due to thr=0
% Try ALL Thresholds @ Rx Side
for thr=0:max_thr
   rx_node.threshold = thr;
   % Demod & eval BER
   ber_vec(thr+1) = eval_BER(nRx_raw_matrix_wout_noise, tx_sym_matrix, rx_node, sim_params);
end
fprintf(1,'\n\n Demodulation Process Ends');
res.thr_vec = 0:max_thr;
res.ber_vec = ber_vec;
end


function [res] = eval_BER(nRx_raw_matrix_wout_noise, tx_sym_matrix, rx_node, sim_params)

sym_cnt             = size(tx_sym_matrix, 2);
replication         = size(tx_sym_matrix, 1);

demod_sym_matrix    = zeros(replication, sym_cnt);

% For Each Line (i.e., Replication) demodulate bits
    for jj=1:replication
       nRx_singleRun_wout_noise = nRx_raw_matrix_wout_noise(:,:,jj);
       demod_sym_matrix(jj, :) = CORE_block_demodulate(nRx_singleRun_wout_noise, rx_node, sim_params);
    end

res = nnz(demod_sym_matrix-tx_sym_matrix) / numel(tx_sym_matrix);

end