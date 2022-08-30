function [] = non_linear_model(r_r, D, D_tx, D_rx, emission_pt, receiver_pt)

%{
syms b1 b2 b3 tao t V_rx d_tx_rx_1 d_tx_rx_2 D1 D2
model_f_Rx_1 = b1*V_rx*erfc(d_tx_rx_1/((4*D2)^b2*(D1/D2*tao+t)^b3))/(4*pi*D2*d_tx_rx_1);
model_f_Rx_2 = b1*V_rx*erfc(d_tx_rx_2/((4*D2)^b2*(D1/D2*tao+t)^b3))/(4*pi*D2*d_tx_rx_2);
%}

% diffusion coefficient
diffusion_molecules = D;
diffusion_tx = D_tx;
diffusion_rx = D_rx;
D1 = diffusion_molecules + diffusion_rx;
D2 = diffusion_tx + diffusion_rx;

% initial distance from Tx to Rx1 and Rx2
tx_center = emission_pt;
rx_center = receiver_pt;
d_tx_rx_1 = norm(rx_center(1,:) - tx_center) - r_r;
d_tx_rx_2 = norm(rx_center(2,:) - tx_center) - r_r;

V_rx = 4*pi*r_r^3/3;

%% mobile fitted coefficient [500 10 10] topology [d = 10 r = 4 h = 2]
load 'molecules_data_mobile[10 4 2][500 10 10].mat'
X = 1e-10:0.001:2;
figure;
plot(X(1:end), nRx_1_avg(1:end), '--r', 'LineWidth', 2)
hold on
plot(X(1:end), nRx_2_avg(1:end), '--g', 'LineWidth', 2)
nRx_1 = nRx_1_avg;
nRx_2 = nRx_2_avg;
y_1 = cumsum(nRx_1);
y_2 = cumsum(nRx_2);
y_1 = y_1./5000;
y_2 = y_2./5000;
model_f_Rx_1 = @(b,x)b(1)*V_rx.*erfc(d_tx_rx_1./((4*D2)^b(2).*(D1/D2.*x+1).^b(3)))./(4*pi*D2*d_tx_rx_1);
model_f_Rx_2 = @(b,x)b(1)*V_rx.*erfc(d_tx_rx_2./((4*D2)^b(2).*(D1/D2.*x+1).^b(3)))./(4*pi*D2*d_tx_rx_2);
% r = a + (b-a).*rand(N,1)
b_init = 0 + (1-0).*rand(3,1);
Rx_1_mdl = fitnlm(X(1:1000),y_1(1001:end),model_f_Rx_1,b_init);
Rx_2_mdl = fitnlm(X(1:1000),y_2(1001:end),model_f_Rx_2,b_init);
b_Rx_1 = Rx_1_mdl.Coefficients.Estimate;
b_Rx_2 = Rx_2_mdl.Coefficients.Estimate;
fprintf(1, '\n############ mobiled TX_1----RX_1 fitted parameter ############\n');
disp(b_Rx_1(:,1));
fprintf(1, '\n############ mobiled TX_1----RX_2 fitted parameter ############\n');
disp(b_Rx_2(:,1));

fitted_Rx_1 = 5000*b_Rx_1(1,1)*V_rx.*erfc(d_tx_rx_1./((4*D2)^b_Rx_1(2,1).*(D1/D2.*X+1).^b_Rx_1(3,1)))./(4*pi*D2*d_tx_rx_1);
fitted_Rx_2 = 5000*b_Rx_2(1,1)*V_rx.*erfc(d_tx_rx_2./((4*D2)^b_Rx_2(2,1).*(D1/D2.*X+1).^b_Rx_2(3,1)))./(4*pi*D2*d_tx_rx_2);
p = plot(X(1:end), nRx_1(1:end), X(1:end), nRx_2(1:end), X(1:end-1), diff(fitted_Rx_1), X(1:end-1),diff(fitted_Rx_2));
% Tx1-Rx1 simulation
p(1).Color = '#0072BD';
p(1).LineWidth = 2.5;
p(1).LineStyle = '--';
% Tx1-Rx2 simulation
p(2).Color = '#4DBEEE';
p(2).LineWidth = 2.5;
p(2).LineStyle = '--';
% Tx1-Rx1 fitted
p(3).Color = '#D95319';
p(3).LineWidth = 1.5;
p(3).LineStyle = '-.';
% Tx1-Rx2 fitted
p(4).Color = '#EDB120';
p(4).LineWidth = 1.5;
p(4).LineStyle = '-.';
hold on ;

%% no mobile fitted diffusion coefficient [50 0 0]
%{
load molecules_data_no_mobile.mat;
X = time;
nRx_1 = nRx_1_avg;
nRx_2 = nRx_2_avg;
y_1 = cumsum(nRx_1);
y_2 = cumsum(nRx_2);
y_1 = y_1./5000;
y_2 = y_2./5000;
model_f_Rx_1 = @(b,x)erfc(d_tx_rx_1./((4*diffusion_molecules)^b(2).*(x).^b(3))).*b(1).*r_r./(d_tx_rx_1+r_r);
model_f_Rx_2 = @(b,x)erfc(d_tx_rx_2./((4*diffusion_molecules)^b(2).*(x).^b(3))).*b(1).*r_r./(d_tx_rx_2+r_r);
b_init = 0 + (1-0).*rand(3,1);
Rx_1_mdl = fitnlm(X,y_1,model_f_Rx_1,b_init);
Rx_2_mdl = fitnlm(X,y_2,model_f_Rx_2,b_init);
b_Rx_1 = Rx_1_mdl.Coefficients.Estimate;
b_Rx_2 = Rx_2_mdl.Coefficients.Estimate;
fprintf(1, '\n########### fixed TX_1----RX_1 fitted parameter ############\n');
disp(b_Rx_1(:,1));
fprintf(1, '\n########### fixed TX_1----RX_2 fitted parameter ############\n');
disp(b_Rx_2(:,1));
fitted_Rx_1 = 5000*erfc(d_tx_rx_1./((4*diffusion_molecules)^b_Rx_1(2).*(X).^b_Rx_1(3))).*b_Rx_1(1).*r_r./(d_tx_rx_1+r_r);
fitted_Rx_2 = 5000*erfc(d_tx_rx_2./((4*diffusion_molecules)^b_Rx_2(2).*(X).^b_Rx_2(3))).*b_Rx_2(1).*r_r./(d_tx_rx_2+r_r);
p = plot(X(1:end), nRx_1(1:end),  X(1:end), nRx_2(1:end),  X(1:end-1), diff(fitted_Rx_1), X(1:end-1),diff(fitted_Rx_2));
% Tx1-Rx1 simulation
p(1).Color = '#0000FF';
p(1).LineWidth = 3;
p(1).LineStyle = '--';
% Tx1-Rx2 simulation
p(2).Color = '#FF00FF';
p(2).LineWidth = 3;
p(2).LineStyle = '--';
% Tx1-Rx1 fitted
p(3).Color = '#000000';
p(3).LineWidth = 1.5;
p(3).LineStyle = '-.';
% Tx1-Rx2 fitted
p(4).Color = '#FFFFFF';
p(4).LineWidth = 1.5;
p(4).LineStyle = '-.';
%}
end