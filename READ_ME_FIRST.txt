New Features
- Added multiple types of molecule emission and BMoSK modulation support
- Renamed all CORE functions
- Added new EXAMPLES for guidance


Obsolete Features
- 1D, 2D support (those can be implemented very easily)


Important Points (IMPORTANT IMPORTANT IMPORTANT)
- For the simulation of arriving molecules, there are two different styles: 
one is particle tracking based and the other is characteristic function based.
If you are using particle tracking based simulator, note that it will be slower
and you may need to set delta_t to 0.001. If you are using characteristic 
function based simulator, then you may set delta_t to half of the ts (symbol duration).

- P2S means point source (transmitter) to spherical receiver
- S2S means spherical transmitter to spherical receiver

- tss (sampling in a symbol durationmay be more than 1) is important for MFSK 
in which you emit molecules during a symbol duration not just at the start of 
the symbol duration.

- If MFSK is used you need to have a header and trailer bit for synchronization

- If you report bugs (specific cases) I can handle it very quickly, but if you 
report non specific requests, they may have a journey to null heaven. 

- If you want to add new modulations, you should edit 
CORE_block_modulate
CORE_block_demodulate

- To understand the simulator well, you should walk through 
EXAMPLE_runner_00, EXAMPLE_runner_01, and EXAMPLE_runner_02 in order

Important Dimensions (Standardization)
% OUTPUT of single simulation
% Added DIFFERENT molecule TYPES by considering each row as another molecule type
nRx_wout_noise = zeros (mol_type_cnt, sim_step_cnt);

% OUTPUT of sim_replicator
% Dimensions : (nRx_row_cnt_molType x nsym*ts_step x replication)
nRx_raw_matrix_wout_noise = zeros(nRx_row_cnt_molType, nsym*ts_step, replication);

