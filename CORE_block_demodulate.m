%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2014, 
% Chan-Byoung Chae (Yonsei University)
% H. Birkan YILMAZ (Yonsei University)
% All rights reserved.
%
% Updated and extended versions of this simulator can be found at:
% http://www.cbchae.org/
% http://www.cmpe.boun.edu.tr/~yilmaz/
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

function [rx_syms ] = CORE_block_demodulate( nRx_tss_sampled, rx_node, sim_params )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rx_sym_wout_noise           : Demodulated Symbol sequence from wout_noise
% rx_sym_with_noise           : Demodulated Symbol sequence from with_noise
%
% rx_node.synch_offset     : Time synchronization shift for received signal (it must be betwwen 0 and ts in miliseconds)
% rx_node.demod            : DeModulator index (see block_demodulate())
% rx_node.threshold        : Threshold for CSK-based demodulations (MFSK based is hardly coded)
%
% sim_params.delta_t                      : Simulation step time
% sim_params.molecules_perTs              : Number of molecules per Ts
% sim_params.ts_inMiliseconds             : Symbol duration (Ts)
% sim_params.tss_inMiliseconds            : Sampling duration (Tss)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (rx_node.demod == 0)
   rx_syms  = demod_bcsk_pulse(nRx_tss_sampled, rx_node, sim_params);
elseif (rx_node.demod == 1)
   rx_syms = demod_bcsk_square(nRx_tss_sampled, rx_node, sim_params);
elseif (rx_node.demod == 2)
   rx_syms = demod_bmfsk_cosine(nRx_tss_sampled, sim_params);
elseif (rx_node.demod == 2.1)
   rx_syms = demod_bmfsk_cosine_wGuard(nRx_tss_sampled, sim_params, 3/16);
elseif (rx_node.demod == 3)
   rx_syms = demod_bmosk_pulse(nRx_tss_sampled, sim_params);
elseif (rx_node.demod == 4)
   rx_syms = demod_brsk_pulse(nRx_tss_sampled, sim_params);
else
   error('\nUnsupported DeMODULATION index !!!\n');
end

end




function [ rx_syms ] = demod_bcsk_pulse( nRx_tss_sampled, rx_node, sim_params )
% Precondition : Head and Tail are trimmed
% NO ISI Filtering is implemented

threshold         = rx_node.threshold;
ts                = sim_params.ts_inSeconds;
tss               = sim_params.tss_inSeconds;

sampling_cnt_for_symbol = round(ts / tss);
nsym = round( size(nRx_tss_sampled,2) / sampling_cnt_for_symbol );

rx_syms = zeros(1, nsym);
for i=1:nsym
   current_signal = nRx_tss_sampled((i-1)*sampling_cnt_for_symbol+1 : i*sampling_cnt_for_symbol);
   
   tstat = sum(current_signal);
   % Give decision for wout_noise
   if (tstat > threshold)
       rx_syms(1,i) = 1;
   else
       rx_syms(1,i) = 0;
   end 
end % end of <for i=1:nsym>

end


function [ rx_syms ] = demod_bcsk_square( nRx_tss_sampled, rx_node, sim_params )
   rx_syms = demod_bcsk_pulse(nRx_tss_sampled, rx_node, sim_params);
end





function [ rx_syms ] = demod_bmosk_pulse( nRx_tss_sampled, sim_params )
% Precondition : Head and Tail are trimmed
% NO ISI Filtering is implemented

ts                      = sim_params.ts_inSeconds;
tss                     = sim_params.tss_inSeconds;

sampling_cnt_for_symbol = round(ts / tss);
nsym                    = round( size(nRx_tss_sampled,2) / sampling_cnt_for_symbol );

rx_syms = zeros(1, nsym);
for i=1:nsym
   current_signal = nRx_tss_sampled(:, (i-1)*sampling_cnt_for_symbol+1 : i*sampling_cnt_for_symbol);
   
   % Take Sum to Get Detection TestStatistics
   tstat = sum(current_signal, 2);
   
   % Give decision for wout_noise
   if (tstat(1) > tstat(2))
       rx_syms(1,i) = 0;
   elseif (tstat(2) > tstat(1))
       rx_syms(1,i) = 1;
   else
       rx_syms(1,i) = randi([0,1]);
   end
end % end of <for i=1:nsym>

end




function [ rx_syms ] = demod_brsk_pulse( nRx_tss_sampled, sim_params )
% Precondition : Head and Tail are trimmed
% NO ISI Filtering is implemented

ts                      = sim_params.ts_inSeconds;
tss                     = sim_params.tss_inSeconds;

sampling_cnt_for_symbol = round(ts / tss);
nsym                    = round( size(nRx_tss_sampled,2) / sampling_cnt_for_symbol );

rx_syms = zeros(1, nsym);
for i=1:nsym
   current_signal = nRx_tss_sampled(:, (i-1)*sampling_cnt_for_symbol+1 : i*sampling_cnt_for_symbol);
   
   % Take Sum to Get Detection TestStatistics
   tstat = sum(current_signal, 2);
   
   % Give decision for wout_noise
   if (tstat(1) > tstat(2))
       rx_syms(1,i) = 0;
   elseif (tstat(2) > tstat(1))
       rx_syms(1,i) = 1;
   else
       rx_syms(1,i) = randi([0,1]);
   end
end % end of <for i=1:nsym>

end






function [ rx_syms ] = demod_bmfsk_cosine( nRx_tss_sampled, sim_params )
% Precondition : Head and Tail are trimmed

ts                      = sim_params.ts_inSeconds;
tss                     = sim_params.tss_inSeconds;
delta_t                 = sim_params.delta_t;

sampling_cnt_for_symbol = round(ts / tss);

ts_step = round(ts / delta_t);
tss_step = round(tss / delta_t);

nsym = round( size(nRx_tss_sampled,2) / sampling_cnt_for_symbol );

rx_syms = zeros(1, nsym);

% PreFilter for Smoothing
input_signal = imfilter(nRx_tss_sampled, [0.12 0.76 0.12], 'replicate');

f0 = 1.0/ts_step;
f1 = 2*f0;
n  = (0:tss_step:(ts_step-tss_step));
symbol_filters = 1+cos(2*pi*repmat([f0;f1],1,sampling_cnt_for_symbol).* [n;n]);

for i=1:nsym
   current_signal = input_signal((i-1)*sampling_cnt_for_symbol+1 : i*sampling_cnt_for_symbol);
   
   % Evaluate Test Stats
   filtered_output = symbol_filters * current_signal';
   
   % Decide the current symbol
   [~, I] = max(filtered_output);
   rx_syms(1,i)  = I(1)-1;
   
end % end of <for i=1:nsym>

end



function [ rx_syms ] = demod_bmfsk_cosine_wGuard( nRx_tss_sampled, sim_params, guardRatio )
% Precondition : Head and Tail are trimmed


ts                      = sim_params.ts_inSeconds;
tss                     = sim_params.tss_inSeconds;
delta_t                 = sim_params.delta_t;

sampling_cnt_for_symbol = round(ts / tss);
guard_cnt = round(sampling_cnt_for_symbol * guardRatio);

ts_step = round(ts / delta_t);
tss_step = round(tss / delta_t);

emission_cnt = sampling_cnt_for_symbol - guard_cnt;

nsym = round( size(nRx_tss_sampled,2) / sampling_cnt_for_symbol );


rx_syms = zeros(1, nsym);


% PreFilter for Smoothing
input_signal = imfilter(nRx_tss_sampled, [0.12 0.76 0.12], 'replicate');

f0 = 1/(ts_step-guard_cnt*tss_step);
f1 = 2*f0;

tss_index = 0:emission_cnt-1;
sym0_filter = [1 + cos(2*pi*f0*tss_index*tss_step) zeros(1,guard_cnt)] ;
sym1_filter = [1 + cos(2*pi*f1*tss_index*tss_step) zeros(1,guard_cnt)] ;

symbol_filters = [sym0_filter; sym1_filter];

for i=1:nsym
   current_signal = input_signal((i-1)*sampling_cnt_for_symbol+1 : i*sampling_cnt_for_symbol);
   
   % Evaluate Test Stats
   filtered_output = symbol_filters * current_signal';
   
   % Decide the current symbol
   [~, I] = max(filtered_output);
   rx_syms(1,i)  = I(1)-1;
   
end % end of <for i=1:nsym>

end






