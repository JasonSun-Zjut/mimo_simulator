function [ res_Signal ] = EXAMPLE_runner_00_SingleBurst( )
fprintf(1, '\n################# Single Emission  Scenario Parameters  ############');
delta_t              = 0.001;
fprintf(1,'\n ## delta_t = %f s', delta_t);
%Tx1距离Rx1的距离
dist_inMicroMeters   = 5;
fprintf(1,'\n ## Distance = %d micrometer', dist_inMicroMeters);

num_molecules_to_emit = 2000;
fprintf(1,'\n ## Num. Emitted Molecules = %d ', num_molecules_to_emit);

nsym                 = 1;
ts                   = 1;
fprintf(1, '\n################# Single Emission  Scenario Parameters  ############');
fprintf(1, '\n################# Tx2 Parameters  ############');
% Prepare Variables
emission_pt = [0 0 0];
D = 100;
D_rx = 100;
r_r = 1;
%dist, emission_pt, r_r, D, D_rx, delta_t, molecules_perTs, ts_inSeconds, tss_inSeconds, symbol_probs, nsym, replication%
[tx_node, rx_node, env_params, sim_params] = prepare_vars4runners_PointSrc(dist_inMicroMeters, emission_pt, r_r, D, D_rx, delta_t, num_molecules_to_emit, ts, 0.001, [0.5 0.5], nsym, 1000);

% Run 
res_Signal = runner_00_SingleBurst( tx_node, rx_node, env_params, sim_params);
% Plot Arrival Rate
signal_resolution_merge = 2;
x_time = (delta_t*signal_resolution_merge):(delta_t*signal_resolution_merge):2;
y_signal_avg = sum(reshape( res_Signal.nRx_avg, signal_resolution_merge, size(res_Signal.nRx_avg,2)/signal_resolution_merge));
y_signal_n   = sum(reshape( res_Signal.nRx_raw_matrix_wout_noise(:,:,1), signal_resolution_merge, size(res_Signal.nRx_avg,2)/signal_resolution_merge));
figure(1);
plot(x_time(1:1000), y_signal_avg(1:1000), '--r', 'LineWidth', 2)
hold on


%Analytical
Vobs = 4*pi*r_r^3/3;
r0 = dist_inMicroMeters+r_r;
%relative time of observation for a fixed value of t
tao = 1e-10:delta_t:2;
m = num_molecules_to_emit*Vobs.*exp(-r0^2./(4*(D+D_rx).*tao))./(4*pi*(D+D_rx).*tao).^(3/2);
tao = 1e-10:(delta_t*signal_resolution_merge):2;
m_avg = sum(reshape( m, signal_resolution_merge, size(m,2)/signal_resolution_merge));
figure(1)
plot(tao(1:end), m_avg(1:end), '--g', 'LineWidth', 2)
%{
grid on
plot(x_time(1:100), y_signal_n(1:100), '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.7])
set(gca, 'fontsize',11, 'fontweight','d');
xlabel('Time')
ylabel('Num. Received Molecules')
legend('Avg Signal', '1st realization')
%}
res_Signal.plot1_x_time          = x_time;
res_Signal.plot1_y_signal_avg    = y_signal_avg;
res_Signal.plot1_y_signal_n      = y_signal_n;


% Plot Cumulative Arrivals
%{
signal_resolution_merge = 2;
x_time2 = (delta_t*signal_resolution_merge):(delta_t*signal_resolution_merge):2;
% [0,t]内累计接收到的分子个数
y_cuml_signal_avg = cumsum( sum(reshape( res_Signal.nRx_avg, signal_resolution_merge, size(res_Signal.nRx_avg,2)/signal_resolution_merge)) );
y_cuml_signal_n   = cumsum( sum(reshape( res_Signal.nRx_raw_matrix_wout_noise(:,:,1), signal_resolution_merge, size(res_Signal.nRx_avg,2)/signal_resolution_merge)) );
figure;
plot(x_time2(1:1000), y_cuml_signal_avg(1:1000), '--r', 'LineWidth', 2)
hold on
plot(x_time2(1:1000), y_cuml_signal_n(1:1000), '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.7])
set(gca, 'fontsize',11, 'fontweight','d');
xlabel('Time')
ylabel('Cumulative Num. Received Molecules')
legend('Avg Signal', '1st realization')
%}
end


function [res] = runner_00_SingleBurst(tx_node, rx_node, env_params, sim_params)
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

res.nRx_avg =  sum(nRx_raw_matrix_wout_noise, 3) / size(nRx_raw_matrix_wout_noise, 3);

end
