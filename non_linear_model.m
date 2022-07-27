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
y = nRx_1_avg;
y = y./2000;
model_fun = @(b,x)b(1)*V_rx.*exp(-d_tx_rx_1^2./(4*D1.*x))./((4*pi*D1.*x).^(3/2));
b_init = zeros(1,1);
b_init(1) = 1;
mdl = fitnlm(X,y,model_fun,b_init);
end