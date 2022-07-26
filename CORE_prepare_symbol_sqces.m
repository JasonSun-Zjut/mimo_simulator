function [ tx_sym_matrix ] = CORE_prepare_symbol_sqces(symbol_probs, replication, nsym)
% Works for any number of symbol types
% Adds 1symbol at the head and the tail
random_mat = random('Uniform', 0, 1, [replication, nsym]);

% define probability regions on (0,1)
lb = 0;
ub = 0;
current_symbol = 0;
for prob=symbol_probs
   ub = ub + prob;
   
   % between lb and ub is the current symbol
   bigger_than_lb = random_mat > lb;
   smaller_than_ub = random_mat < ub;
   random_mat(bigger_than_lb & smaller_than_ub) = current_symbol;
   
   current_symbol = current_symbol + 1;
   lb = ub;
end

tx_sym_matrix = random_mat;

end