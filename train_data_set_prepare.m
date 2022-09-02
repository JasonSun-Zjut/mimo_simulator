function [] = train_data_set_prepare(num, node, r, n, NA, T, tau, diffusion_coefficient, distance_low_bound, distance_upper_bound, molecule_low_bound, molecule_upper_bound)
%   此处显示详细说明
% NA : 发送节点释放的分子个数
% n  : calculate Pe at nth time slot
% molecule_low_bound: 释放个数的下界
% molecule_upper_bound: 释放个数的上界
% num: 数据个数

% 初始化发送节点到接受节点的距离
rand_distance = round(distance_low_bound + (distance_upper_bound - distance_low_bound) * rand(node, num));
%% Experiment
% 实验参数准备
emit_point = round(rand(1,2)*10);
rx_1_point = rand_coordinate_generate(emit_point, rand_distance(1));
rx_2_point = rand_coordinate_generate(emit_point, rand_distance(2));
coordinate_vector = [emit_point; rx_1_point; rx_2_point];
[Rx_1_estimate_coefficient, Rx_2_estimate_coefficient] = Example_runner_diffusion_rx(n, NA,coordinate_vector, r, diffusion_coefficient);

molecule_to_allocated = sum(round(molecule_low_bound + (molecule_upper_bound - molecule_low_bound) * rand(node, 1)));
for i = 1 : num
    s = sum(rand_distance(1 : node, i));
    distance = rand_distance(1 : node, i);
    molecule_allocated = rand_distance(1 : node, i).' / s * molecule_to_allocated;
    estimate_coefficient = [Rx_1_estimate_coefficient; Rx_2_estimate_coefficient];
    [local_optimum_molecule_allocated, local_optimum] = fminsearch_util(distance, n, molecule_allocated, r, diffusion_coefficient, estimate_coefficient, T, tau, molecule_low_bound, molecule_upper_bound);
end
end

