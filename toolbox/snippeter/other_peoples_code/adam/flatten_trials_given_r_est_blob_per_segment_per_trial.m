function [date_str_flat, ...
          letter_str_flat, ...
          i_segment_within_trial_flat, ...
          localized_flat, ...
          r_est_flat, ...
          r_head_from_video_flat, ...
          r_tail_from_video_flat] = ...
  flatten_trials_given_r_est_blob_per_segment_per_trial(r_est_blob_per_segment_per_trial)

% marshall the values across trials into double arrays (not cell
% arrays)
%K=4;  % number of mics
date_str_flat=cell(0,1);
letter_str_flat=cell(0,1);
i_segment_within_trial_flat=zeros(0,1);
localized_flat=false(0,1);
r_est_flat=zeros(2,0);
r_head_from_video_flat=zeros(2,0);
r_tail_from_video_flat=zeros(2,0);
n_trials=length(r_est_blob_per_segment_per_trial);
for i_trial=1:n_trials
  %n_segments_this_trial=length(r_est_blob_per_segment_per_trial{i_trial});

  date_str_flat_this={r_est_blob_per_segment_per_trial{i_trial}.date_str}';
  date_str_flat= ...
    [date_str_flat;date_str_flat_this];  %#ok

  letter_str_flat_this={r_est_blob_per_segment_per_trial{i_trial}.letter_str}';
  letter_str_flat= ...
    [letter_str_flat;letter_str_flat_this];  %#ok
  
  i_segment_within_trial_flat= ...
    [i_segment_within_trial_flat;[r_est_blob_per_segment_per_trial{i_trial}.i_segment]'];  %#ok

  localized_flat= ...
    [localized_flat;[r_est_blob_per_segment_per_trial{i_trial}.localized]'];  %#ok
  
  r_est_flat= ...
    [r_est_flat [r_est_blob_per_segment_per_trial{i_trial}.r_est] ];  %#ok

  r_head_from_video_flat= ...
    [r_head_from_video_flat [r_est_blob_per_segment_per_trial{i_trial}.r_head_from_video] ];  %#ok
  
  r_tail_from_video_flat= ...
    [r_tail_from_video_flat [r_est_blob_per_segment_per_trial{i_trial}.r_tail_from_video] ];  %#ok
end

end
