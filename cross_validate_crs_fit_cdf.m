% load pre-calculated stuff
load('r_est_raw_for_single_mouse_data.mat');
% r_est_blob_per_voc_per_trial
r_est_blob_per_voc_per_trial_raw=r_est_blob_per_voc_per_trial;
clear r_est_blob_per_voc_per_trial;

% clear out vocs in which r_est is close to mic mean
% these are usually vocs with very low SNR
r_est_blob_per_voc_per_trial= ...
  delete_vocs_near_mic_mean(r_est_blob_per_voc_per_trial_raw, ...
                            per_trial_ancillary);
                
% extract the fields I want to cell arrays
mse_body=extract_field_from_blob(r_est_blob_per_voc_per_trial,'mse_body');
mse_min=extract_field_from_blob(r_est_blob_per_voc_per_trial,'mse_min');
                          
% calc J at the body
n_trials=length(mse_body);
dJ_body=cell(n_trials,1);
for i=1:n_trials
  dJ_body{i}=(mse_body{i}-mse_min{i})./mse_min{i};
end

% make the grid for the cdf
dJ_min=0;
dJ_max=1;
ddJ_want=0.01;
dJ_line=(dJ_min:ddJ_want/50:dJ_max)';

% Do CV for 68% CRs
conf_level=0.68
P_coverage_cv=cross_validation_fit_cdf(dJ_body,dJ_line,conf_level)
P_coverage_cv_mean=mean(P_coverage_cv)
P_coverage_cv_sd=std(P_coverage_cv)

% Do CV for 95% CRs
conf_level=0.95
P_coverage_cv=cross_validation_fit_cdf(dJ_body,dJ_line,conf_level)
P_coverage_cv_mean=mean(P_coverage_cv)
P_coverage_cv_sd=std(P_coverage_cv)

