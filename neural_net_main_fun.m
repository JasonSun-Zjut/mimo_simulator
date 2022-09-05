%%
clear;
node = 2;
net = feedforwardnet([2 * node + 2, 2 * node, 2 * node]);
net.trainfcn='trainbr';
net.layers{2}.transferFcn='softmax';
net.layers{3}.transferFcn='poslin';
% num, node, r, n, NA, T, tau, diffusion_coefficient, distance_low_bound, distance_upper_bound, molecule_low_bound, molecule_upper_bound
% 初始数据集个数
num = 10;
% 接收节点半径
r = 3;
%分子扩散试验释放个数
NA = 5000;
T = 1;
tau = 0.005;
% 分子扩散系数、传输节点扩散系数、接收节点扩散系数
diffusion_coefficient = [50, 10, 10];
distance_low_bound = 10;
distance_upper_bound = 15;
molecule_low_bound = 2000;
molecule_upper_bound = 5000;

data_set = train_data_set_prepare(num, node, r, n, NA, T, tau, diffusion_coefficient, distance_low_bound, distance_upper_bound, molecule_low_bound, molecule_upper_bound);

net = train(net, data_set(1:end,1:node)', data_set(1:end, node+1:2*node)');