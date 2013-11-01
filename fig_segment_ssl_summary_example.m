% calculate stuff for all vocalizations
options.verbosity=0;  % how much output or intermediate results the user wants to 
                      % see
options.read_from_map_cache=false;  % whether to try to use the map cache 
                                    % to save time
options.write_to_map_cache=false;  % whether to write to the map cache after 
                                   % calculating a map de novo
options.quantify_confidence=false;  % calculate P-vals, CRs (that's
                                 % what makes it not "raw")
options.return_big_things=true;  % return the full map and other large
                                 % data structures

% identifying info for the segment
date_str='06132012';
letter_str='D';
i_segment=51;  % this was voc84 in the old-style

% directories where to find stuff
base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_sys_test';
data_analysis_dir_name='Data_analysis10';

% load the trial overhead
trial_overhead=...
  ssl_trial_overhead_cartesian_heckbertian_packaged(base_dir_name,...
                                                    data_analysis_dir_name, ...
                                                    date_str, ...
                                                    letter_str);

% package up the args                                           
args=trial_overhead;
args.i_segment=i_segment;
%args.verbosity=verbosity;

% call the functions
r_est_blobs = ...
  r_ests_from_segment_indicators_and_trial_overhead(args,options);

% unpack the return blob
%i_snippet_pretty=16;
field_names=fieldnames(r_est_blobs);
for i=1:length(field_names)
  eval(sprintf('%s_all_snippets={r_est_blobs(:).%s}'';',field_names{i},field_names{i}));
end
clear r_est_blobs;

% transform things from cell arrays what can
n_snippets=length(tf_rect_name_all_snippets);  
r_est_all_snippets=cell2mat(reshape(r_est_all_snippets,[1 n_snippets]));  %#ok
r_head_from_video_all_snippets=cell2mat(reshape(r_head_from_video_all_snippets,[1 n_snippets]));  %#ok
r_tail_from_video_all_snippets=cell2mat(reshape(r_tail_from_video_all_snippets,[1 n_snippets]));  %#ok

% unpack the trial overhead
R=trial_overhead.R;
Temp=trial_overhead.Temp;
dx=trial_overhead.dx;
x_grid=trial_overhead.x_grid;
y_grid=trial_overhead.y_grid;
in_cage=trial_overhead.in_cage;
fs=trial_overhead.fs;
r_corners=trial_overhead.r_corners;
clear trial_overhead;

% get dims out
N=size(v_all_snippets{1},1);  %#ok % all same length
n_mics=size(R,2);
dt=1/fs;  % s
n_snippets=length(tf_rect_name_all_snippets)  %#ok
n_pairs=nchoosek(n_mics,2);

% need to do outlier filtering on r_est here
[is_outlier,~,r_est_trans,Covariance_matrix] = kur_rce(r_est_all_snippets',1);
is_outlier=logical(is_outlier);
r_est=r_est_trans';

% Get the video position from the snippets
i_start_all_snippets=cell2mat(i_start_all_snippets);  %#ok
i_end_all_snippets=cell2mat(i_end_all_snippets);  %#ok
[r_head_from_video,r_tail_from_video]= ...
  r_head_for_segment_from_snippets(r_head_from_video_all_snippets, ...
                                   r_tail_from_video_all_snippets, ...
                                   i_start_all_snippets, ...
                                   i_end_all_snippets);                                 

% make up some mouse locations
n_fake_mice=3;
%[r_head_from_video_fake,r_tail_from_video_fake]= ...
%  random_mouse_locations(R,r_head_from_video,r_tail_from_video,n_fake_mice);
% These were a nice-looking sample:
r_head_from_video_fake = ...
   [     0.445553008346868         0.709662308265758         0.522328094818333 ; ...
         0.334817684474456         0.419620179150835         0.582494615220904 ];
r_tail_from_video_fake = ...
   [     0.473832772083876         0.636316477631333         0.556192459384125 ; ...
         0.257019105608912         0.381243714500405         0.658029830339911 ];

% Have to transform these last to convential Cartesian coords
r_head_from_video_fake(2,:)=(0.67925-r_head_from_video_fake(2,:))+0.0265;
r_tail_from_video_fake(2,:)=(0.67925-r_tail_from_video_fake(2,:))+0.0265;

       
% caclulate the density at the real+fake mice, and the posterior
% probability
r_head_from_video_with_fake=[r_head_from_video r_head_from_video_fake];
r_tail_from_video_with_fake=[r_tail_from_video r_tail_from_video_fake];
% r_chest_from_video_real_and_fake= ...
%   (3/4)*r_head_from_video_with_fake + ...
%   (1/4)*r_tail_from_video_with_fake ;
% p_chest_real_and_fake=mvnpdf(r_chest_from_video_real_and_fake',r_est',Covariance_matrix);  % density
% P_posterior_chest_from_video_real_and_fake= ...
%   p_chest_real_and_fake/sum(p_chest_real_and_fake);  % posterior probability


% plot the per-snippet estimates, the density, and the mice
colorbar_max=100;  % 1/m^2
are_mice_beyond_first_fake=true;
[h_fig,h_axes,h_axes_cb]= ...
  fig_segment_ssl_summary(r_est,Covariance_matrix, ...
                          r_est_all_snippets,is_outlier, ...
                          R,r_corners, ...
                          r_head_from_video_with_fake,r_tail_from_video_with_fake, ...
                          colorbar_max, ...
                          are_mice_beyond_first_fake);
                        



