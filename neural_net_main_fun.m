%%
clear;
node = 2;
net = feedforwardnet([2 * node + 2, 2 * node, 2 * node]);
net.trainfcn='trainbr';
net.layers{2}.transferFcn='softmax';
net.layers{3}.transferFcn='poslin';
molecule_allocate_low_bound = 5000;
molecule_allocate_upper_bound = 10000;
distance_low_bound = 10;
distance_upper_bound = 30;
Init_Data_Number = 20;
Rand_Input = sort(randi([distance_low_bound, distance_upper_bound], Init_Data_Number, node), 2);
Rand_Output = randi([molecule_allocate_low_bound, molecule_allocate_upper_bound], 20, node);
Rand_Data_Set = [Rand_Input, Rand_Output];
% init network
net = train(net, Rand_Input', Rand_Output');

% train network
Train_Data_Number = 100;
Train_Input = sort(randi([distance_low_bound, distance_upper_bound], Train_Data_Number, node), 2);
for i = 1 : size(Train_Input, 1)
     net_output = sim(net, Train_Input(i,:)');
     sub_range_length = 100;
     [trail_solution, freq_variable] = Diversification_Generation(node, molecule_allocate_low_bound, molecule_allocate_upper_bound, sub_range_length, Train_Data_Number);
     % error probabiltiy calculate
end










%%
%{
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
diffusion_coefficient = [5000, 1, 1];
molecule_low_bound = 2000;
molecule_upper_bound = 5000;
%data_set = train_data_set_prepare(node, r, n, NA, T, tau, diffusion_coefficient, distance, molecule_low_bound, molecule_upper_bound);
%net = train(net, data_set(1:end,1:node)', data_set(1:end, node+1:2*node)');
for i = 1 : 10
rand_input = round(10 + (30 - 10) * rand(2, 1));
if(rand_input(1) > rand_input(2))
    temp = rand_input(1);
    rand_input(1) = rand_input(2);
    rand_input(2) = temp;
end
point = [0,0];
rx_1_point = rand_coordinate_generate(point,rand_input(1));
rx_2_point = rand_coordinate_generate(point,rand_input(2));
coordinate_vector = [[point,0]; [rx_1_point,0]; [rx_2_point,0]];
[Rx_1_estimate_coefficient, Rx_2_estimate_coefficient] = Example_runner_diffusion_rx(1, NA, coordinate_vector, r, diffusion_coefficient);
estimate_coefficient = [Rx_1_estimate_coefficient, Rx_2_estimate_coefficient];
% 通过ann获得分配方案
molecule_number_by_ann = sim(net, rand_input);
Pe = BER(molecule_number_by_ann, estimate_coefficient, diffusion_coefficient, T, tau, 5, 3, rand_input);
[local_optimum_molecule_allocated, local_optimum] = fminsearch_util(rand_input, 5, molecule_number_by_ann, 3, diffusion_coefficient, estimate_coefficient, T, tau, molecule_low_bound, molecule_upper_bound);
end
function Pe = BER(molecule_allocated, estimate_coefficient, diffusion_coefficient, T, tau, n, r, distance)
        % molecule_allocated: 围绕初始点寻找是否存在局部极小值
        molecule_allocated = round(molecule_allocated);
        mu_0 = 0;
        for k = 1:1:n-1
        %[Rx1_estimate_coefficient, Rx2_estimate_coefficient] = Example_runner_diffusion_rx(k, NA, coordinate_vector, r, diffusion_coefficient);
        probability_function = probability_mobile(estimate_coefficient(1:end, 1), estimate_coefficient(1:end, 2), r, diffusion_coefficient, (k-1)*T,(n-k)*T+tau,distance);
            mu_0 = mu_0 + 0.5 * molecule_allocated(1) * probability_function(1) +  0.5 * molecule_allocated(2) * probability_function(2);
        end
        %[Rx1_estimate_coefficient, Rx2_estimate_coefficient] = Example_runner_diffusion_rx(n, NA, coordinate_vector, r, diffusion_coefficient);
        probability_function = probability_mobile(estimate_coefficient(1:end, 1), estimate_coefficient(1:end, 2), r, diffusion_coefficient, (n-1)*T, tau, distance);
        mu_0 = 0.5 * molecule_allocated(2) * probability_function(2);
        %% Tx_1 transmit bit 0
        mu_1 = mu_0 + molecule_allocated(1) * probability_function(1);
        threshold = ceil((mu_1+mu_0)/log(mu_1/mu_0));
        Pe = 0.5 * ( (1 - poisscdf(threshold, mu_0)) + poisscdf(threshold, mu_1));
end
%}