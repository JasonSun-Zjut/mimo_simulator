function [data_set] = train_data_set_prepare(num, node, r, n, NA, T, tau, diffusion_coefficient, distance_low_bound, distance_upper_bound, molecule_low_bound, molecule_upper_bound)
%   此处显示详细说明
% NA : 发送节点释放的分子个数
% n  : calculate Pe at nth time slot
% molecule_low_bound: 释放个数的下界
% molecule_upper_bound: 释放个数的上界
% num: 数据个数

% 初始化发送节点到接受节点的距离
rand_distance = round(distance_low_bound + (distance_upper_bound - distance_low_bound) * rand(node, num));
data_set = zeros(num,5);
for i = 1 : num
while(rand_distance(1, i) == rand_distance(2, i))
    rand_distance(1:2, i) = round(distance_low_bound + (distance_upper_bound - distance_low_bound) * rand(node, 1));
end
fprintf(1, '\n############ distance:[%d, %d] ############\n',rand_distance(1, i), rand_distance(2, i));
%% Experiment
% 实验参数准备
emit_point = round(rand(1,2)*10);
rx_1_point = rand_coordinate_generate(emit_point, rand_distance(1, i));
rx_2_point = rand_coordinate_generate(emit_point, rand_distance(2, i));
coordinate_vector = [[emit_point,0]; [rx_1_point,0]; [rx_2_point,0]];
fprintf(1, '\n############ point:[%d, %d, %d] ############\n',coordinate_vector(1,1:end), coordinate_vector(2,1:end), coordinate_vector(3,1:end));
[Rx_1_estimate_coefficient, Rx_2_estimate_coefficient] = Example_runner_diffusion_rx(1, NA,coordinate_vector, r, diffusion_coefficient);

molecule_to_allocated = [round((molecule_low_bound + molecule_upper_bound)/2), round((molecule_low_bound + molecule_upper_bound)/2)];
%s = sum(rand_distance(1 : node, i));
distance = rand_distance(1 : node, i);
molecule_allocated = molecule_to_allocated;
estimate_coefficient = [Rx_1_estimate_coefficient, Rx_2_estimate_coefficient];
fprintf(1, '\n############ molecule_allocated:[%d, %d] ############\n',molecule_allocated(1), molecule_allocated(2));
[local_optimum_molecule_allocated, local_optimum] = fminsearch_util(distance, n, molecule_allocated, r, diffusion_coefficient, estimate_coefficient, T, tau, molecule_low_bound, molecule_upper_bound);
data_set(i, 1:end) = [rand_distance(1 : node, i)', local_optimum_molecule_allocated, local_optimum];
end
end

