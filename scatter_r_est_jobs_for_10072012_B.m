% calculate stuff for all vocalizations
verbosity=0;  % how much output or intermediate results the user wants to 
              % see
args.read_from_map_cache=false;  % whether to try to use the map cache 
                                % to save time
args.write_to_map_cache=false;  % whether to write to the map cache after 
                               % calculating a map de novo
args.quantify_confidence=false;  % calculate P-vals, CRs (that's
                                 % what makes it not "raw")
args.return_big_things=false;  % don't return the full map or other large
                               % data structures

% identifying info for each trial
date_str=cell(0,1);
letter_str=cell(0,1);
date_str{end+1}='10072012';
letter_str{end+1}='B';
n_vocs_per_trial_max=20;  % max number of vocs to do per trial
%n_vocs_per_trial_max=inf;  % max number of vocs to do per trial

% directories where to find stuff
base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_vocal_structure';
data_analysis_dir_name='Data_analysis2';

% call r_est_from_voc_indicators_and_ancillary() for each voc, collect
% all the results in r_est_blob_per_voc_per_trial
%[r_est_blob_per_voc_per_trial,per_trial_ancillary]= ... = ...
  scatter_multiple_trials(base_dir_name, ...
                          data_analysis_dir_name, ...
                          date_str, ...
                          letter_str, ...
                          n_vocs_per_trial_max, ...
                          @r_est_from_voc_indicators_and_ancillary, ...
                          args, ...
                          verbosity);

% save everything           
%save('r_est_for_10072012_B.mat');

