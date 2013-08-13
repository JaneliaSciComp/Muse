function feval_blob_from_voc_indicators_and_ancillary_executable(output_file_name,input_file_name)

% make sure the common blob_from_voc_indicators_and_ancillary functions get
% included in the executable:
%#function r_est_from_voc_indicators_and_ancillary
%#function r_est_from_voc_indicators_and_ancillary_for_a_few_vocs

load(input_file_name);
if ~exist(output_file_name,'file')         
  blobs = ...
    feval(blob_from_voc_indicators_and_ancillary, ...
          base_dir_name,date_str_this_trial,letter_str_this_trial, ...
          args_common, ...
          argses, ...
          verbosity);  %#ok
  % save(output_file_name,'r_est','dssen_head');
  save(output_file_name,'blobs','args_common','argses','base_dir_name','date_str_this_trial','letter_str_this_trial');
end

end
