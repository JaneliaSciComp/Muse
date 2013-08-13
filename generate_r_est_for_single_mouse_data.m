% calculate stuff for all vocalizations
verbosity=0;  % how much output or intermediate results the user wants to 
              % see
args.read_from_map_cache=true;  % whether to try to use the map cache 
                                % to save time
args.write_to_map_cache=true;  % whether to write to the map cache after 
                               % calculating a map de novo
args.quantify_confidence=true;  % calculate P-vals, CRs (that's
                                 % what makes it not "raw")
args.return_big_things=false;  % don't return the full map or other large
                               % data structures

% identifying info for each trial
date_str=cell(0,1);
letter_str=cell(0,1);
date_str{end+1}='06052012';
letter_str{end+1}='D';
date_str{end+1}='06062012';
letter_str{end+1}='E';
date_str{end+1}='06102012';
letter_str{end+1}='E';
date_str{end+1}='06112012';
letter_str{end+1}='D';
date_str{end+1}='06122012';
letter_str{end+1}='D';
date_str{end+1}='06122012';
letter_str{end+1}='E';
date_str{end+1}='06132012';  % this is the one with the least vocs
letter_str{end+1}='D';
date_str{end+1}='06132012';
letter_str{end+1}='E';
n_vocs_per_trial_max=inf;  % max number of vocs to do per trial

% directories where to find stuff
base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_sys_test';
data_analysis_dir_name='Data_analysis';

% call r_est_from_voc_indicators_and_ancillary() for each voc, collect
% all the results in r_est_blob_per_voc_per_trial
[r_est_blob_per_voc_per_trial,per_trial_ancillary]= ... = ...
  map_multiple_trials(base_dir_name, ...
                      data_analysis_dir_name, ...
                      date_str, ...
                      letter_str, ...
                      n_vocs_per_trial_max, ...
                      @r_est_from_voc_indicators_and_ancillary, ...
                      args, ...
                      verbosity);

% save everything           
save('r_est_for_single_mouse_data.mat');

