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

function [ nRx_wout_noise ] = CORE_sim_diffusion_3d_P2S_FAST( ...
   tx_timeline,     ...% Modulated timeline to transmit
   mol_type_cnt,    ...% Molecule type count
   tx_node,         ...% Tx node properties
   rx_node,         ...% Rx node properties
   env_params,      ...% Environment properties
   sim_params )   % Simulation parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates diffusion channel in 3D env
% with a POINT source and an absorbing SPHERICAL receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tx_timeline              : Modulated timeline to transmit with dimensions [mol_type_cnt x sim_step_count] 
%
% tx_node.emission_point   : Coordinates of emission point
% tx_node.mod              : Modulator index (see block_modulate())
%
% rx_node.center           : Center Coordinates of rx_node center
% rx_node.r_inMicroMeters  : Radius of rx_node
% rx_node.demod            : DeModulator index (see block_demodulate())
%
% env_params.D_inMicroMeterSqrPerSecond   : Diffusion coefficient
% env_params.destruction_limit            : Destruction boundary
% 
% sim_params.delta_t                      : Simulation step time
% sim_params.molecules_perTs              : Number of molecules per Ts
% sim_params.ts_inSeconds                 : Symbol duration (Ts)
% sim_params.tss_inSeconds                : Sampling duration (Tss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

nRx_wout_noise = perfect_absorption(tx_timeline, mol_type_cnt, tx_node, rx_node, env_params, sim_params );

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ nRx_wout_noise ] = perfect_absorption( ...
   tx_timeline,     ...% Modulated timeline to transmit
   mol_type_cnt,    ...% Molecule type count
   tx_node,         ...% Tx node properties
   rx_node,         ...% Rx node properties
   env_params,      ...% Environment properties
   sim_params )   % Simulation parameters


rx_r_inMicroMeters      = rx_node.r_inMicroMeters;
tx_emission_point       = tx_node.emission_point;

D                       = env_params.D_inMicroMeterSqrPerSecond;

ts                      = sim_params.ts_inSeconds;
delta_t                 = sim_params.delta_t;


distance_inMicroMeters  = sqrt( sum((tx_emission_point-rx_node.center).^2, 2) ) - rx_r_inMicroMeters;

% First find the number of simulation steps
sim_step_cnt = size(tx_timeline,2);


% Window of ISI
window_sym_cnt    = 25;
p_window          = eval_phit_window(window_sym_cnt, distance_inMicroMeters, D, rx_r_inMicroMeters, ts, delta_t);


% Rx timeline Records the number of molecules at RECEIVER at each time step 
% Added DIFFERENT molecule TYPES by considering each row as another molecule type
nRx_wout_noise = zeros (mol_type_cnt, sim_step_cnt);



for t=1:sim_step_cnt
   % Check for Emission for EACH MOL_TYPE
   num_release = tx_timeline(:, t);
   
   if (sum(num_release) > 0)
      % Add new molecules to environment before moving them
      for ii=1:mol_type_cnt
          if (num_release(ii) > 0)
              rand_arrivals      = mybinornd_wpvec(num_release(ii), p_window);
              sz_rand_arrivals   = size(rand_arrivals, 2);
              
              t_start   = t;
              t_end     = t+sz_rand_arrivals-1;
              
              if (t_end > sim_step_cnt)
                 t_end = sim_step_cnt;
              end
              
              nRx_wout_noise(ii, t_start:t_end) = nRx_wout_noise(ii, t_start:t_end) + rand_arrivals(1:(t_end-t_start+1));
          end
      end
   end
   
end % end.of   <for t=1:sim_step_cnt>


end



function [ res_p_window ] = eval_phit_window( window_sym_cnt, distance_inMicroMeters, D, rx_r_inMicroMeters, ts, delta_t)

t_start  = 0:delta_t:(window_sym_cnt*ts-delta_t);
t_end    = delta_t:delta_t:(window_sym_cnt*ts);

multiplier = rx_r_inMicroMeters / (rx_r_inMicroMeters+distance_inMicroMeters);
res_p_window = multiplier * (erfc(distance_inMicroMeters./(4*D*t_end).^0.5) - erfc(distance_inMicroMeters./(4*D*t_start).^0.5));

end

function [ res ] = mybinornd_wpvec( N, p_vec )
% N should be an integer

col_cnt = size(p_vec,2);

res = zeros(1, col_cnt);

for jj=1:col_cnt
    res(jj) = sum(rand(1,N)<p_vec(jj));
end

end


