function [res_Signal] = Example_runner_diffusion_passive_rx(n, NA,transmitter_coordinate,receiver_coordinate,r,diffusion_coefficient)
fprintf(1, '\n################# 开始模拟Tx在 %dth时隙释放分子的过程  ############', n);
delta_t              = 5*10^(-6);
fprintf(1,'\n ## delta_t = %f s', delta_t);
num_molecules_to_emit = NA;
fprintf(1,'\n ## Num. Emitted Molecules = %d ', num_molecules_to_emit);
nsym                 = n;
ts                   = 0.5*10^(-3);
dist_inMicroMeters = norm(receiver_coordinate - transmitter_coordinate);
[tx_node, rx_node, env_params, sim_params] = prepare_vars4_diffusion_runners_PointSrc(dist_inMicroMeters, transmitter_coordinate, receiver_coordinate, r, diffusion_coefficient(1), diffusion_coefficient(2), diffusion_coefficient(3), delta_t, num_molecules_to_emit, ts, delta_t, [0.5 0.5], nsym, 200);
%mobile perfect passive
rx_node.p_react  = 4;
res_Signal = runner_diffusion_passive_rx(n, tx_node, rx_node, env_params, sim_params);
end

function [res] = runner_diffusion_passive_rx(n, tx_node, rx_node, env_params, sim_params)
tx_sym_matrix              = repmat([ones(1,n-1),1], sim_params.replication, 1);

tx_node.mod                = 0; %% BCSK (pulse)
rx_node.demod              = tx_node.mod;

[ nRx_raw_matrix_wout_noise, stats ] = CORE_sim_replicator( ...
           'CORE_sim_diffusion_3d_P2S_wAbsorption', ...% Simulator Name ##You Can Change This##
           tx_sym_matrix, ... % Tx Symbol Sequences
           tx_node,    ...% Tx node properties
           rx_node,    ...% Rx node properties
           env_params, ...% Environment properties
           sim_params );

res.nRx_raw_matrix_wout_noise = nRx_raw_matrix_wout_noise;
res.stats = stats;
res.nRx_avg =  sum(nRx_raw_matrix_wout_noise(:,:,:), 3) / size(nRx_raw_matrix_wout_noise(:,:,:), 3);
end
