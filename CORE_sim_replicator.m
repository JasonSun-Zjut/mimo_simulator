%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2014, 
% Chan-Byoung Chae (Yonsei University)
% H. Birkan YILMAZ (Yonsei University)
% All rights reserved.
%
% Updated and extended versions of this simulator can be found at:
% http://www.cmpe.boun.edu.tr/~yilmaz/
% http://www.cbchae.org/
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the Yonsei University.
% 4. Neither the name of the organization nor the
%    names of its contributors may be used to endorse or promote products
%    derived from this software without specific prior written permission.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ nRx_raw_matrix_wout_noise, stats ] = CORE_sim_replicator( ...
   sim_function_name,   ...% Simulator Name
   tx_sym_matrix, ...% Tx Symbol Sequences for each replication
   tx_node,    ...% Tx node properties
   rx_node,    ...% Rx node properties
   env_params, ...% Environment properties
   sim_params )   % Simulation parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tx_node.emission_point   : Coordinates of emission point
% tx_node.mod              : Modulator index (see block_modulate())
%
% rx_node.center           : Center Coordinates of rx_node center
% rx_node.r_inMicroMeters  : Radius of rx_node
% rx_node.synch_offset     : Time synchronization shift for received signal 
%                            (it must be in the interval [0, ts-delta_t] in terms of delta_t miliseconds)
% rx_node.demod            : DeModulator index (see block_demodulate())
%
% env_params.D_inMicroMeterSqrPerSecond   : Diffusion coefficient
% env_params.destruction_limit            : Destruction boundary
% 
% sim_params.delta_t                      : Simulation step time
% sim_params.molecules_perTs              : Number of molecules per Ts
% sim_params.ts_inSeconds                 : Symbol duration (Ts)
% sim_params.tss_inSeconds                : Sampling duration (Tss)
% sim_params.nsym                         : Symbol sequence length (N_sym)
% sim_params.replication                  : How many different symbol sqces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% OUTPUT Value:
% 
% nRx_raw_matrix_wout_noise               : Number of received molecules (without noise)
%                                           for each symbol sequence of
%                                           tx_sym_matrix(rep x nsym)
%
% Dimensions : (nRx_row_cnt_molType x nsym*ts_step x replication)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart_all = tic;

replication     = sim_params.replication;

ts_step         = round( sim_params.ts_inSeconds / sim_params.delta_t );
nsym            = size(tx_sym_matrix, 2); 
rx_number = size(rx_node.center,1);

if ((tx_node.mod == 3) || (tx_node.mod == 4))
    nRx_row_cnt_molType = 2;
else
    nRx_row_cnt_molType = 1;
end


nRx_raw_matrix_wout_noise = zeros(nRx_row_cnt_molType, rx_number*nsym*ts_step, replication);

simFunc = str2func(sim_function_name);


tElapsed_diffusion = 0;
for i=1:replication
   % Take the ith row for the current random symbol sqces [1,0]
   tx_sym_seq = tx_sym_matrix(i,:);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%   MODULATE  --> tx_timeline[mol_type_cnt x sim_step_cnt]
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%   Tx timeline (each row corresponds different molecule type)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [tx_timeline, mol_type_cnt] = CORE_block_modulate(tx_sym_seq, tx_node, sim_params);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%   SIMULATE DIFFUSION  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   tStart_diffusion = tic;
   nRx_raw_matrix_wout_noise(:, :, i) = simFunc(tx_timeline, mol_type_cnt, tx_node, rx_node, env_params, sim_params );
   tElapsed_diffusion = tElapsed_diffusion + toc(tStart_diffusion);
   
   
   if (mod(i,10)==1)
      fprintf('\nReplication *** (%d / %d)', i, replication);
   end
end
fprintf('\nReplication *** (%d / %d)', replication, replication);
fprintf('\n##########################\nPreparing Timing Stats');
tElapsed_all = toc(tStart_all);

stats.tElapsed_all         = tElapsed_all;
stats.tElapsed_diffusion   = tElapsed_diffusion;

fprintf('\nPreparing Timing Stats  [DONE]\n');
end





%%%%%%%%%%%%%%%%%%%%%%%%%%% UTILITY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

