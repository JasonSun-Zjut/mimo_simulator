function [ tx_node, rx_node, env_params, sim_params ] = prepare_vars4runners_PointSrc(dist, emission_pt, r_r, D,D_rx, delta_t, molecules_perTs, ts_inSeconds, tss_inSeconds, symbol_probs, nsym, replication)

dist_inMicroMeters = dist;

tx_node.emission_point     = emission_pt;

rx_node.r_inMicroMeters    = r_r;
rx_node.center             = [rx_node.r_inMicroMeters+dist_inMicroMeters 0 0];
%%rx_node.center = [10 0 0];
rx_node.synch_offset       = 0; % Param
rx_node.p_react            = 1;
rx_node.D_inMicroMeterSqrPerSecond = D_rx;
env_params.snr_db                         = 30; % NOISE FROM Environment  
env_params.D_inMicroMeterSqrPerSecond     = D;
env_params.destruction_limit              = 7*dist_inMicroMeters + 45 + rx_node.r_inMicroMeters;

sim_params.delta_t                        = delta_t;
sim_params.molecules_perTs                = molecules_perTs;   % Param changes in Example Runner
sim_params.ts_inSeconds                   = ts_inSeconds; % Param changes in Example Runner
sim_params.tss_inSeconds                  = tss_inSeconds; % Param changes in Example Runner

sim_params.symbol_probs    = symbol_probs;
sim_params.nsym            = nsym; 
sim_params.replication     = replication; % Increase Replication


end

