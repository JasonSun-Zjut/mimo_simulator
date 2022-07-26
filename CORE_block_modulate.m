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

function [ res_timeline, mol_type_cnt ] = CORE_block_modulate( tx_sym_seq, tx_node, sim_params )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tx_sym_seq         : Symbol sequence to modulate
% tx_node            : Modulator index 
%                       0-> BCSK (pulse)
%                       1-> BCSK (square)
%                       2-> BMFSK (cosine)  2.1-> BMFSK wG
%                       3-> BMoSK (pulse)
%                       4-> BRSK 
%
% sim_params.delta_t                      : Simulation step time
% sim_params.molecules_perTs              : Number of molecules per Ts
% sim_params.ts_inSeconds                 : Symbol duration (Ts)
% sim_params.tss_inSeconds                : Sampling duration (Tss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (tx_node.mod == 0)
   res_timeline = mod_bcsk_pulse(tx_sym_seq, sim_params);
   mol_type_cnt = 1;
elseif (tx_node.mod == 1)
   res_timeline = mod_bcsk_square(tx_sym_seq, sim_params);
   mol_type_cnt = 1;
elseif (tx_node.mod == 2)
   res_timeline = mod_bmfsk_cosine(tx_sym_seq, sim_params);
   mol_type_cnt = 1;
elseif (tx_node.mod == 2.1)
   res_timeline = mod_bmfsk_cosine_withGuard(tx_sym_seq, sim_params, 3/16);
   mol_type_cnt = 1;
elseif (tx_node.mod == 3)
   res_timeline = mod_bmosk_pulse(tx_sym_seq, sim_params);
   mol_type_cnt = 2;
elseif (tx_node.mod == 4)
   res_timeline = mod_brsk_pulse(tx_sym_seq, sim_params);
   mol_type_cnt = 2;
else
   error('\nUnsupported MODULATION index !!!\n');
end

end




function [ res_timeline ] = mod_bcsk_pulse( tx_sym_seq, sim_params )
n_sym                   = size(tx_sym_seq, 2);

molecules_perTs         = sim_params.molecules_perTs;
ts                      = sim_params.ts_inSeconds;
delta_t                 = sim_params.delta_t;
sim_time_inSeconds      = n_sym * ts;

% First find the number of simulation steps
sim_step_cnt = round(sim_time_inSeconds / delta_t);

ts_step = round(ts / delta_t);

% Tx timeline allocation
res_timeline = zeros (1, sim_step_cnt);

for i=1:n_sym
   if (tx_sym_seq(i) == 1)
       %release all the molecules at the start of symbol duration%
      res_timeline( :, (i-1)*ts_step+1 ) = molecules_perTs;
   end
end

end




function [ res_timeline ] = mod_bmosk_pulse( tx_sym_seq, sim_params )
n_sym                   = size(tx_sym_seq, 2);

molecules_perTs         = sim_params.molecules_perTs;
ts                      = sim_params.ts_inSeconds;
delta_t                 = sim_params.delta_t;
sim_time_inSeconds      = n_sym * ts;

% First find the number of simulation steps
sim_step_cnt = round(sim_time_inSeconds / delta_t);

ts_step = round(ts / delta_t);

% Tx timeline allocation
res_timeline = zeros (2, sim_step_cnt);

for i=1:n_sym
   if (tx_sym_seq(i) == 0)     % SYM:0 TYPE A Molecule @ index:1
      res_timeline(1, (i-1)*ts_step+1 ) = molecules_perTs;
   elseif (tx_sym_seq(i) == 1) % SYM:1 TYPE B Molecule @ index:2
      res_timeline(2, (i-1)*ts_step+1 ) = molecules_perTs;
   end
end

end




function [ res_timeline ] = mod_brsk_pulse( tx_sym_seq, sim_params )
n_sym                   = size(tx_sym_seq, 2);

molecules_perTs         = sim_params.molecules_perTs;
ts                      = sim_params.ts_inSeconds;
delta_t                 = sim_params.delta_t;
sim_time_inSeconds      = n_sym * ts;

% First find the number of simulation steps
sim_step_cnt = round(sim_time_inSeconds / delta_t);

ts_step = round(ts / delta_t);

% Tx timeline allocation
res_timeline = zeros (2, sim_step_cnt);

for i=1:n_sym
   if (tx_sym_seq(i) == 0)     % SYM:0 TYPE A Molecule @ index:1
      res_timeline(1, (i-1)*ts_step+1 ) = round(0.85*molecules_perTs);
      res_timeline(2, (i-1)*ts_step+1 ) = round(0.15*molecules_perTs);
   elseif (tx_sym_seq(i) == 1) % SYM:1 TYPE B Molecule @ index:2
      res_timeline(1, (i-1)*ts_step+1 ) = round(0.15*molecules_perTs);
      res_timeline(2, (i-1)*ts_step+1 ) = round(0.85*molecules_perTs);
   end
end

end





function [ res_timeline ] = mod_bcsk_square( tx_sym_seq, sim_params )
n_sym                   = size(tx_sym_seq,2);

molecules_perTs         = sim_params.molecules_perTs;
ts                      = sim_params.ts_inSeconds;
tss                     = sim_params.tss_inSeconds;
delta_t                 = sim_params.delta_t;
sim_time_inSeconds      = n_sym * ts;

% First find the number of simulation steps
sim_step_cnt = round(sim_time_inSeconds / delta_t);

ts_step = round(ts / delta_t);
tss_step = round(tss / delta_t);

emission_cnt = ts_step / tss_step;
molecules_perTss = round(molecules_perTs / emission_cnt);

sym1 = reshape( [molecules_perTss*ones(1,emission_cnt); zeros(tss_step-1, emission_cnt)] , 1,ts_step);

% Tx timeline allocation
res_timeline = zeros (1, sim_step_cnt);
for i=1:n_sym
   if (tx_sym_seq(i) == 1)
      res_timeline( (i-1)*ts_step+1 : i*ts_step ) = sym1;
   end
end
end




function [ res_timeline ] = mod_bmfsk_cosine( tx_sym_seq, sim_params )
n_sym                   = size(tx_sym_seq,2);

molecules_perTs         = sim_params.molecules_perTs;
ts                      = sim_params.ts_inSeconds;
tss                     = sim_params.tss_inSeconds;
delta_t                 = sim_params.delta_t;
sim_time_inSeconds      = n_sym * ts;

% First find the number of simulation steps
sim_step_cnt = round(sim_time_inSeconds / delta_t);

ts_step = round(ts / delta_t);
tss_step = round(tss / delta_t);

emission_cnt = ts_step / tss_step;
avg_molecules_perTss = round(molecules_perTs / emission_cnt);

f0 = 1.0/ts_step;
f1 = 2*f0;

tss_index = 0:emission_cnt-1;
sym0 = reshape( [round(avg_molecules_perTss + avg_molecules_perTss*cos(2*pi*f0*tss_index*tss_step)); zeros(tss_step-1, emission_cnt)] , 1,ts_step);
sym1 = reshape( [round(avg_molecules_perTss + avg_molecules_perTss*cos(2*pi*f1*tss_index*tss_step)); zeros(tss_step-1, emission_cnt)] , 1,ts_step);


% Tx timeline allocation
res_timeline = zeros (1, sim_step_cnt);
for i=1:n_sym
   if (tx_sym_seq(i) == 0)
      res_timeline( (i-1)*ts_step+1 : i*ts_step ) = sym0;
   elseif (tx_sym_seq(i) == 1)
      res_timeline( (i-1)*ts_step+1 : i*ts_step ) = sym1;
   end
end
end




function [ res_timeline ] = mod_bmfsk_cosine_withGuard( tx_sym_seq, sim_params, guardRatio )
n_sym                   = size(tx_sym_seq,2);

molecules_perTs         = sim_params.molecules_perTs;
ts                      = sim_params.ts_inSeconds;
tss                     = sim_params.tss_inSeconds;
delta_t                 = sim_params.delta_t;
sim_time_inSeconds      = n_sym * ts;


% First find the number of simulation steps
sim_step_cnt = round(sim_time_inSeconds / delta_t);

ts_step = round(ts / delta_t);
tss_step = round(tss / delta_t);

sampling_cnt = ts_step / tss_step;
guard_cnt = round(sampling_cnt * guardRatio);

emission_cnt = sampling_cnt - guard_cnt;
avg_molecules_perTss = round(molecules_perTs / emission_cnt);



f0 = 1/(ts_step-guard_cnt*tss_step);
f1 = 2*f0;

tss_index = 0:emission_cnt-1;
sym0 = reshape( [round(avg_molecules_perTss + avg_molecules_perTss*cos(2*pi*f0*tss_index*tss_step)) zeros(1,guard_cnt); zeros(tss_step-1, emission_cnt+guard_cnt)] , 1,ts_step);
sym1 = reshape( [round(avg_molecules_perTss + avg_molecules_perTss*cos(2*pi*f1*tss_index*tss_step)) zeros(1,guard_cnt); zeros(tss_step-1, emission_cnt+guard_cnt)] , 1,ts_step);


% Tx timeline allocation
res_timeline = zeros (1, sim_step_cnt);
for i=1:n_sym
   if (tx_sym_seq(i) == 0)
      res_timeline( (i-1)*ts_step+1 : i*ts_step ) = sym0;
   elseif (tx_sym_seq(i) == 1)
      res_timeline( (i-1)*ts_step+1 : i*ts_step ) = sym1;
   end
end
end




