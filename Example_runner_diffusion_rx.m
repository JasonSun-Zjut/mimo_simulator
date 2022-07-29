function [ res_Signal] = Example_runner_diffusion_rx()
fprintf(1, '\n################# Single Emission  Scenario Parameters  ############');
delta_t              = 0.001;
fprintf(1,'\n ## delta_t = %f s', delta_t);
num_molecules_to_emit = 5000;
fprintf(1,'\n ## Num. Emitted Molecules = %d ', num_molecules_to_emit);
nsym                 = 1;
ts                   = 1;
fprintf(1, '\n################# Single Emission  Scenario Parameters  ############');
fprintf(1, '\n################# Tx2 Parameters  ############');
% Prepare Variables
emission_pt = [0 0 0]; 
receiver_pt = [[8 0 0];[8 -9 0]];
dist_inMicroMeters = 4;
%dist_inMicroMeters, emission_pt, receiver_pt, r_r, D, D_tx, D_rx, delta_t, molecules_perTs, ts_inSeconds, tss_inSeconds, symbol_probs, nsym, replication%
[tx_node, rx_node, env_params, sim_params] = prepare_vars4_diffusion_runners_PointSrc(dist_inMicroMeters, emission_pt, receiver_pt, 4, 50, 0, 0, delta_t, num_molecules_to_emit, ts, 0.001, [0.5 0.5], nsym, 500);
% Run 
res_Signal = runner_diffusion_rx( tx_node, rx_node, env_params, sim_params);
% Plot Arrival Rate
signal_resolution_merge = 1;
x_time = (delta_t*signal_resolution_merge):(delta_t*signal_resolution_merge):2;
%y_1_signal_avg = sum(reshape( res_Signal.nRx_1_avg, signal_resolution_merge, size(res_Signal.nRx_1_avg,2)/signal_resolution_merge));
%y_2_signal_avg = sum(reshape( res_Signal.nRx_2_avg, signal_resolution_merge, size(res_Signal.nRx_2_avg,2)/signal_resolution_merge));
figure(1);
plot(x_time(1:end), res_Signal.nRx_1_avg(1:end), '--r', 'LineWidth', 2)
hold on
plot(x_time(1:end), res_Signal.nRx_2_avg(1:end), '--g', 'LineWidth', 2)

%% Analytical
%{
r_r = 1;
Vobs = 4*pi*r_r^3/3;
D = 100;
D_rx = 30;
r_1_0 = norm(receiver_pt(1,:) - emission_pt);
r_2_0 = norm(receiver_pt(2,:) - emission_pt);
%relative time of observation for a fixed value of t
tao = (delta_t*signal_resolution_merge):delta_t:2;
m_1 = num_molecules_to_emit*Vobs.*exp(-r_1_0^2./(4*(D+D_rx).*tao))./(4*pi*(D+D_rx).*tao).^(3/2);
m_2 = num_molecules_to_emit*Vobs.*exp(-r_2_0^2./(4*(D+D_rx).*tao))./(4*pi*(D+D_rx).*tao).^(3/2);
%tao = 1e-10:(delta_t*signal_resolution_merge):2;
%m_1_avg = sum(reshape( m_1, signal_resolution_merge, size(m_1,2)/signal_resolution_merge));
%m_2_avg = sum(reshape( m_2, signal_resolution_merge, size(m_2,2)/signal_resolution_merge));
figure(1)
plot(tao(1:end), m_1(1:end), '-b', 'LineWidth', 1)
plot(tao(1:end), m_2(1:end), '-k', 'LineWidth', 1)
%}
end

function [res] = runner_diffusion_rx(tx_node, rx_node, env_params, sim_params)
tx_sym_matrix              = repmat([1,0], sim_params.replication, 1);

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

res.nRx_1_avg =  sum(nRx_raw_matrix_wout_noise(:,1:2000,:), 3) / size(nRx_raw_matrix_wout_noise(:,1:2000,:), 3);
res.nRx_2_avg =  sum(nRx_raw_matrix_wout_noise(:,2001:end,:), 3) / size(nRx_raw_matrix_wout_noise(:,2001:end,:), 3);
end