function [r_est_blobs,overhead]= ...
  r_est_from_segment_indicators(base_dir_name,...
                                data_analysis_dir_name, ...
                                date_str, ...
                                letter_str, ...
                                i_segment, ...
                                are_positions_in_old_style_coords, ...
                                frame_height_in_pels, ...
                                options, ...
                                trial_overhead, ...
                                verbosity)

% do all the stuff we only have to do once per trial, unless it's been
% pre-computed
if isempty(trial_overhead) ,
  [syl_name_all,i_start_all,i_end_all,f_lo_all,f_hi_all, ...
   r_head_all,r_tail_all,R,Temp, ...
   dx,x_grid,y_grid,in_cage]= ...
    ssl_trial_overhead(base_dir_name, ...
                       data_analysis_dir_name, ...
                       date_str, ...
                       letter_str, ...
                       are_positions_in_old_style_coords, ...
                       frame_height_in_pels);
else
  % unpack the contents of trial overhead
  field_names=fieldnames(trial_overhead);
  for i=1:length(field_names) ,
    eval(sprintf('%s=trial_overhead.(field_names{i});',field_names{i}));
  end
end  

% filter out the snippets (i.e. syls) that don't match i_segment
if ~isempty(syl_name_all) ,
  syl_name_all_first=syl_name_all{1};
  if length(syl_name_all_first)~=17 ,
    error('The "syllable" "names" appear to be in the wrong format.');
  end
  i_segment_all_as_string=cellfun(@(s)s(4:9),syl_name_all,'UniformOutput',false);
  i_segment_all=str2double(i_segment_all_as_string);
  keep=(i_segment_all==i_segment);
  syl_name_keep=syl_name_all(keep);
  i_start_keep=i_start_all(keep);
  i_end_keep=i_end_all(keep);
  f_lo_keep=f_lo_all(keep);
  f_hi_keep=f_hi_all(keep);
  r_head_keep=r_head_all(:,keep,:);
  r_tail_keep=r_tail_all(:,keep,:);    
end

% package up all the trial overhead for return
overhead.R=R;
overhead.Temp=Temp;
overhead.dx=dx;
overhead.x_grid=x_grid;
overhead.y_grid=y_grid;
overhead.in_cage=in_cage;

% % pack up the arguments that don't change across snippets
% args_template_and_overhead=merge_scalar_structs(options,overhead);

% do the estimation for each snippet
n_snippets=length(syl_name_keep);
r_est_blobs=struct();
for i_snippet=1:n_snippets ,
  % pack up all the arguments
  %snippet_args=args_template_and_overhead;
  syl_name=syl_name_keep{i_snippet};
  i_start=i_start_keep(i_snippet);
  i_end=i_end_keep(i_snippet);  
  f_lo=f_lo_keep(i_snippet);  
  f_hi=f_hi_keep(i_snippet);  
  r_head=r_head_keep(:,i_snippet);  
  r_tail=r_tail_keep(:,i_snippet);

  % estimate r
  r_est_blob_this = ... 
    r_est_from_tf_rect_indicators_and_ancillary(base_dir_name, ...
                                                date_str, ...
                                                letter_str, ...
                                                R,Temp,dx,x_grid,y_grid,in_cage, ...
                                                syl_name, ...
                                                i_start,i_end, ...
                                                f_lo,f_hi, ...
                                                r_head,r_tail, ...
                                                options, ...
                                                verbosity);

  % throw some args into the returned blob
  %args_field_name=fieldnames(args);
  args_field_names_to_add_to_blob={'syl_name' 'i_start' 'i_end' 'f_lo' 'f_hi' 'r_head' 'r_tail'}';
  for i=1:length(args_field_names_to_add_to_blob)
    r_est_blob_this.(args_field_names_to_add_to_blob{i})=eval(sprintf('%s',args_field_names_to_add_to_blob{i}));
  end
  
  if i_snippet==1
    r_est_blobs=r_est_blob_this;
    fs=r_est_blob_this.fs;
    overhead.fs=fs;
  else
    r_est_blobs(i_snippet)=r_est_blob_this;
  end
end

end

