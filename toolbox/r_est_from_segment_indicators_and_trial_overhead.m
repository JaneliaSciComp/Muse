function output= ...
  r_est_from_segment_indicators_and_trial_overhead(args,options)

% Estimates position for all snippets in the given segment, gets rid of
% outliers, and returns the overall position estimate for the segment.

% call the functions
r_est_blobs = ...
  r_ests_from_segment_indicators_and_trial_overhead(args,options);

% unpack the return blob
n_snippets=length(r_est_blobs);
field_names=fieldnames(r_est_blobs);
for i=1:length(field_names)
  eval(sprintf('%s_all_snippets={r_est_blobs.%s}'';',field_names{i},field_names{i}));
end
%clear r_est_blobs;

% transform things from cell arrays what can
if n_snippets<3 ,
  %output=rmfield(args,{'x_grid' 'y_grid' 'in_cage'});
  output=struct();
  output.date_str=args.date_str;
  output.letter_str=args.letter_str;
  output.i_segment=args.i_segment;
  output.localized=false;
  output.r_est=[nan nan]';
  output.Covariance_matrix=[nan nan; nan nan];
  n_mice=size(args.r_head_from_video,3);
  output.p_head=nan(n_mice,1);
  output.P_posterior_head=nan(n_mice,1);
  output.r_head_from_video=nan(2,1);
  output.r_tail_from_video=nan(2,1);
  if options.return_big_things ,
    output.rsrp_grid=nan;
  end
  return
end
if options.return_big_things ,
  rsrp_grid_all_snippets=cell2mat(reshape(rsrp_grid_all_snippets,[1 1 n_snippets]));
end
%rsrp_per_pair_grid_all_snippets=cell2mat(reshape(rsrp_per_pair_grid_all_snippets,[1 1 1 n_snippets]));
r_est_all_snippets=cell2mat(reshape(r_est_all_snippets,[1 n_snippets]));
r_head_from_video_all_snippets=cell2mat(reshape(r_head_from_video_all_snippets,[1 n_snippets]));
r_tail_from_video_all_snippets=cell2mat(reshape(r_tail_from_video_all_snippets,[1 n_snippets]));

% % unpack the trial overhead
% R=trial_overhead.R;
% Temp=trial_overhead.Temp;
% dx=trial_overhead.dx;
% x_grid=trial_overhead.x_grid;
% y_grid=trial_overhead.y_grid;
% in_cage=trial_overhead.in_cage;
% fs=trial_overhead.fs;

% need to do outlier filtering on r_est here
[is_outlier,~,r_est_trans,Covariance_matrix] = kur_rce(r_est_all_snippets',1);
is_outlier=logical(is_outlier);
n_outliers=sum(is_outlier);
n_keepers=n_snippets-n_outliers;
if n_keepers<3 ,
  %output=rmfield(args,{'x_grid' 'y_grid' 'in_cage'});
  output=struct();
  output.date_str=args.date_str;
  output.letter_str=args.letter_str;
  output.i_segment=args.i_segment;
  output.localized=false;
  output.r_est=[nan nan]';
  output.Covariance_matrix=[nan nan; nan nan];
  n_mice=size(args.r_head_from_video,3);
  output.p_head=nan(n_mice,1);
  output.P_posterior_head=nan(n_mice,1);
  output.r_head_from_video=nan(2,1);
  output.r_tail_from_video=nan(2,1);
  if options.return_big_things ,
    output.rsrp_grid=nan;
  end
  return
end

%indices_of_outliers=find(is_outlier);
r_est=r_est_trans';  % overall position estimate for the segment
%n_outliers=sum(is_outlier);
%n_keepers=n_snippets-n_outliers;

% filter out the outliers
is_keeper=~is_outlier;
%n_keepers=sum(is_keeper);
if options.return_big_things ,
  rsrp_grid_all_keepers=rsrp_grid_all_snippets(:,:,is_keeper);
end
%r_est_all_keepers=r_est_all_snippets(:,is_keeper);
%r_est_all_outliers=r_est_all_snippets(:,is_outlier);

% take the mean of the maps for all the non-outliers
if options.return_big_things ,
  rsrp_grid=mean(rsrp_grid_all_keepers,3);
end
r_head_from_video=mean(r_head_from_video_all_snippets,2);
r_tail_from_video=mean(r_tail_from_video_all_snippets,2);

% calculate the density at the mice, and the posterior
% probabilities
p_head=mvnpdf(r_head_from_video',r_est',Covariance_matrix);  % density
P_posterior_head=p_head/sum(p_head);  % posterior probability

% put stuff in the return blob
%output=rmfield(args,{'x_grid' 'y_grid' 'in_cage'});
output=struct();
output.date_str=args.date_str;
output.letter_str=args.letter_str;
output.i_segment=args.i_segment;
output.localized=true;
output.r_est=r_est;
output.Covariance_matrix=Covariance_matrix;
output.p_head=p_head;
output.P_posterior_head=P_posterior_head;
output.r_head_from_video=r_head_from_video;
output.r_tail_from_video=r_tail_from_video;
if options.return_big_things ,
  output.rsrp_grid=nan;
end

end
