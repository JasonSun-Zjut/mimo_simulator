function [x, fval] = fminsearch_util(distance, n, molecule_allocated, r, diffusion_coefficient, estimate_coefficient, T, tau, molecule_low_bound, molecule_upper_bound)
% Input parameter : 
% distance: distance vector between [Tx1,Tx2,...Txn] and Rx
% r : radius of Rx
% diffusion_coefficient: [Dmol, Dtx, Drx]
% estimate_coefficient : model coefficient
%                        from Example_runner_diffusion_rx.m
% NA: number of molecules to emit in experiment
% molecule_allocated: 根据初始距离分配的个数
% Output parameter:
% local_optimum: distance 与 起始点的局部最优点构成数据对，用于neural network训练

%% Numberial analytical
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',molecule_allocated,...
'objective',@BER,'lb',[molecule_low_bound, molecule_low_bound],'ub',[molecule_upper_bound, molecule_upper_bound]);
[x,fval] = run(gs,problem);
fprintf(1, '\n############ Pe : %d ############\n', fval);
% construct data set
function Pe = BER(molecule_allocated)
        % molecule_allocated: 围绕初始点寻找是否存在局部极小值
        molecule_allocated = round(molecule_allocated);
        mu_0 = 0;
        for k = 1:1:n-1
        %[Rx1_estimate_coefficient, Rx2_estimate_coefficient] = Example_runner_diffusion_rx(k, NA, coordinate_vector, r, diffusion_coefficient);
        probability_function = probability_mobile(estimate_coefficient(k, 1:3)', estimate_coefficient(k, 4:end)', r, diffusion_coefficient, (k-1)*T,(n-k)*T+tau,distance);
            mu_0 = mu_0 + 0.5 * molecule_allocated(1) * probability_function(1) +  0.5 * molecule_allocated(2) * probability_function(2);
        end
        %[Rx1_estimate_coefficient, Rx2_estimate_coefficient] = Example_runner_diffusion_rx(n, NA, coordinate_vector, r, diffusion_coefficient);
        probability_function = probability_mobile(estimate_coefficient(n, 1:3)', estimate_coefficient(n, 4:end)', r, diffusion_coefficient, (n-1)*T, tau, distance);
        mu_0 = 0.5 * molecule_allocated(2) * probability_function(2);
        %% Tx_1 transmit bit 0
        mu_1 = mu_0 + molecule_allocated(1) * probability_function(1);
        threshold = ceil((mu_1+mu_0)/log(mu_1/mu_0));
        Pe = 0.5 * ( (1 - poisscdf(threshold, mu_0)) + poisscdf(threshold, mu_1));
end
%{
function stop = output(x, optimvalues, state)
        stop = false;
        if(any(x(:)>molecule_upper_bound)||any(x(:)<molecule_low_bound))
            stop = true;
        end
        if isequal(state,'iter')
          history = [history; x];
        end
end
%}
end