function [] = non_linear_model(r_r, D, D_tx, D_rx, emission_pt, receiver_pt)
% diffusion coefficient
diffusion_molecules = D;
diffusion_tx = D_tx;
diffusion_rx = D_rx;
D1 = diffusion_molecules + diffusion_rx;
D2 = diffusion_tx + diffusion_rx;

% initial distance from Tx to Rx1 and Rx2
tx_center = emission_pt;
rx_center = receiver_pt;
d_tx_rx_1 = norm(rx_center(1,:) - tx_center);
d_tx_rx_2 = norm(rx_center(2,:) - tx_center);

V_rx = 4*pi*r_r^3/3;

load molecules_data.mat;
X = time;
y_1 = nRx_1_avg;
y_2 = nRx_2_avg;
y_1 = y_1./2000;
y_2 = y_2./2000;
plot(X(1:end), y_1(1:end), '-g', 'LineWidth', 2);
hold on;
plot(X(1:end), y_2(1:end), '-b', 'LineWidth', 2);
model_f_Rx_1 = @(b,x)b(1)*V_rx.*exp(-d_tx_rx_1^2./((4*D1)^b(2).*x.^b(3)))./((4*pi*D1)^(3/2).*x.^b(4));
model_f_Rx_2 = @(b,x)b(1)*V_rx.*exp(-d_tx_rx_2^2./((4*D1)^b(2).*x.^b(3)))./((4*pi*D1)^(3/2).*x.^b(4));
% r = a + (b-a).*rand(N,1)
b_init = 1 + (2-1).*rand(4,1);
Rx_1_mdl = fitnlm(X,y_1,model_f_Rx_1,b_init);
Rx_2_mdl = fitnlm(X,y_2,model_f_Rx_2,b_init);
b_Rx_1 = Rx_1_mdl.Coefficients.Estimate;
b_Rx_2 = Rx_2_mdl.Coefficients.Estimate;
fprintf(1, '\n################# TX_1----RX_1 fitted parameter ############\n');
disp(b_Rx_1(:,1));
fprintf(1, '\n################# TX_1----RX_2 fitted parameter ############\n');
disp(b_Rx_2(:,1));
figure(1)
f_Rx_1 = b_Rx_1(1)*V_rx.*exp(-d_tx_rx_1^2./((4*D1)^b_Rx_1(2).*X.^b_Rx_1(3)))./((4*pi*D1)^(3/2).*X.^b_Rx_1(4));
f_Rx_2 = b_Rx_2(1)*V_rx.*exp(-d_tx_rx_2^2./((4*D1)^b_Rx_2(2).*X.^b_Rx_2(3)))./((4*pi*D1)^(3/2).*X.^b_Rx_2(4));
plot(X(1:end), f_Rx_1(1:end), '-.r', 'LineWidth', 1.5);
hold on;
plot(X(1:end), f_Rx_2(1:end), '-.w', 'LineWidth', 1.5);
end