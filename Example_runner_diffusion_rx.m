function [ b_Rx_1, b_Rx_2] = Example_runner_diffusion_rx(n, NA,coordinate_vector,r,diffusion_coefficient)
% Description:
% NA : number of molecules to emit
% coordinate_vector: [[Tx1];[Rx1];[Rx2]]
% r : radius of Rx
% diffusion_coefficient : [D_mol,D_tx,D_rx]
% n: nth time slot emit molecules
%% experiment
fprintf(1, '\n################# 开始模拟Tx在 %dth时隙释放分子的过程  ############', n);
delta_t              = 0.001;
fprintf(1,'\n ## delta_t = %f s', delta_t);
num_molecules_to_emit = NA;
fprintf(1,'\n ## Num. Emitted Molecules = %d ', num_molecules_to_emit);
nsym                 = 1;
ts                   = 1;
% Prepare Variables
emission_pt = coordinate_vector(1,:); 
receiver_pt = [coordinate_vector(2,:);coordinate_vector(3,:)];
dist_inMicroMeters = norm(receiver_pt(1,:)-emission_pt)-r;
%dist_inMicroMeters, emission_pt, receiver_pt, r_r, D, D_tx, D_rx, delta_t, molecules_perTs, ts_inSeconds, tss_inSeconds, symbol_probs, nsym, replication%
[tx_node, rx_node, env_params, sim_params] = prepare_vars4_diffusion_runners_PointSrc(dist_inMicroMeters, emission_pt, receiver_pt, r, diffusion_coefficient(1), diffusion_coefficient(2), diffusion_coefficient(3), delta_t, num_molecules_to_emit, ts, 0.001, [0.5 0.5], nsym, 50);
% Run 
res_Signal = runner_diffusion_rx(n, tx_node, rx_node, env_params, sim_params);
signal_resolution_merge = 1;
x_time = 1e-10:(delta_t*signal_resolution_merge):ts;
%% non_linear_model
nRx_1 = res_Signal.nRx_1_avg;
nRx_2 = res_Signal.nRx_2_avg;
y_1 = cumsum(nRx_1)./NA;
y_2 = cumsum(nRx_2)./NA;
V_rx = 4*pi*r^3/3;
d_tx_rx_1 = norm(receiver_pt(1,:)-emission_pt)-r;
d_tx_rx_2 = norm(receiver_pt(2,:)-emission_pt)-r;
D1 = diffusion_coefficient(1) + diffusion_coefficient(3);
D2 = diffusion_coefficient(2) + diffusion_coefficient(3);
model_f_Rx_1 = @(b,x)b(1)*V_rx.*erfc(d_tx_rx_1./((4*D2)^b(2).*(D1/D2.*x+(n-1)*ts).^b(3)))./(4*pi*D2*d_tx_rx_1);
model_f_Rx_2 = @(b,x)b(1)*V_rx.*erfc(d_tx_rx_2./((4*D2)^b(2).*(D1/D2.*x+(n-1)*ts).^b(3)))./(4*pi*D2*d_tx_rx_2);
b_init = 0.5 + (1-0.5).*rand(3,1);
Rx_1_mdl = fitnlm(x_time(1:end),y_1(1:end),model_f_Rx_1,b_init);
Rx_2_mdl = fitnlm(x_time(1:end),y_2(1:end),model_f_Rx_2,b_init);
b_Rx_1 = Rx_1_mdl.Coefficients.Estimate;
b_Rx_2 = Rx_2_mdl.Coefficients.Estimate;
fprintf(1, '\n############ Tx 在 %dth time slot: mobiled TX_1----RX_1 fitted parameter ############\n', n);
disp(b_Rx_1(:,1));
fprintf(1, '\n############ Tx 在 %dth time slot: mobiled TX_1----RX_2 fitted parameter ############\n', n);
disp(b_Rx_2(:,1));
%% plot
%y_1_signal_avg = sum(reshape( res_Signal.nRx_1_avg, signal_resolution_merge, size(res_Signal.nRx_1_avg,2)/signal_resolution_merge));
%y_2_signal_avg = sum(reshape( res_Signal.nRx_2_avg, signal_resolution_merge, size(res_Signal.nRx_2_avg,2)/signal_resolution_merge));
%figure;
%plot(x_time(1:end), res_Signal.nRx_1_avg(1:end), '--r', 'LineWidth', 2)
%hold on
%plot(x_time(1:end), res_Signal.nRx_2_avg(1:end), '--g', 'LineWidth', 2)


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

function [res] = runner_diffusion_rx(n, tx_node, rx_node, env_params, sim_params)
tx_sym_matrix              = repmat([zeros(1,n-1),1], sim_params.replication, 1);

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
ts_step       =  round( sim_params.ts_inSeconds / sim_params.delta_t );
res.nRx_1_avg =  sum(nRx_raw_matrix_wout_noise(:,(n-1)*ts_step+1:n*ts_step,:), 3) / size(nRx_raw_matrix_wout_noise(:,(n-1)*ts_step+1:n*ts_step,:), 3);
res.nRx_2_avg =  sum(nRx_raw_matrix_wout_noise(:,(end-ts_step)+1:end,:), 3) / size(nRx_raw_matrix_wout_noise(:,(end-ts_step)+1:end,:), 3);
end