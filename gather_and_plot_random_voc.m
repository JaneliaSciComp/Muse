% where the files are stored, etc.
base_dir_name='~/egnor_stuff/ssl_vocal_structure_bizarro';
%base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_vocal_structure';
data_analysis_dir_name='Data_analysis10';
m_file_path=mfilename('fullpath');
project_dir_path=fileparts(m_file_path);
job_outputs_dir_name=fullfile(project_dir_path,'job_outputs');

% pick a random output file, load it
dir_struct=dir(fullfile(job_outputs_dir_name,'*.mat'));
output_file_names_rel={dir_struct.name}';
n_output_file_names=length(output_file_names_rel);
i_output_file_name=floor(n_output_file_names*rand(1))+1;
output_file_name_rel=output_file_names_rel{i_output_file_name};
output_file_name_abs=fullfile(job_outputs_dir_name,output_file_name_rel);
output=load(output_file_name_abs);

% unpack the output
field_name=fieldnames(output);
for i=1:length(field_name)
  eval(sprintf('%s=output.%s;',field_name{i},field_name{i}));
end

% pick a random "voc" from the output file, load it
n_vocs=length(blobs);
i_voc=floor(n_vocs*rand(1))+1;
blob=blobs(i_voc);
args=argses(i_voc);

% unpack the blob
field_name=fieldnames(blob);
for i=1:length(field_name)
  eval(sprintf('%s=blob.%s;',field_name{i},field_name{i}));
end

% unpack the common args
field_name=fieldnames(args_common);
for i=1:length(field_name)
  eval(sprintf('%s=args_common.%s;',field_name{i},field_name{i}));
end

% unpack the args
field_name=fieldnames(args);
for i=1:length(field_name)
  eval(sprintf('%s=args.%s;',field_name{i},field_name{i}));
end

% regenerate the grids
[x_grid,y_grid]=grids(R,dx);

% RMSEs are easier to interpret
%rmse_grid=sqrt(mse_grid);
%rmse_grid=zeros(size(x_grid));
rmse_min=sqrt(mse_min);
rmse_crit=sqrt(mse_crit);
rmse_body=sqrt(mse_body);
rms_total=sqrt(ms_total);

%dsiplay some of those
rmse_min_disp=1e3*rmse_min  %#ok
rms_total_disp=1e3*rms_total  %#ok

% extrapolate a mouse ellipse
r_center=(r_head+r_tail)/2;
a_vec=r_head-r_center;  % vector
b=normcols(a_vec)/3;  % scalar, and a guess at the half-width of the mouse

% Compute the CR borders from the objective function and critical value
%in_cr_grid=(mse_grid<=mse_crit);
%r_cr_bounds_ij=bwboundaries(in_cr_grid);
%r_cr_bounds=path_xy_from_path_ij(r_cr_bounds_ij,x_grid,y_grid);

% make a figure of the objective funtion, estimate, and CR
% title_str=sprintf('%s %s %s', ...
%                   date_str_this_trial,letter_str_this_trial,syl_name);
t0_voc=i_start/fs;  % s  (maybe off by dt)
title_str=sprintf('%s-%s-%s, t0:%0.3f s, i:%d-%d, f:%0.0f-%0.0f', ...
                  date_str_this_trial,letter_str_this_trial,syl_name,t0_voc,i_start,i_end, ...
                  f_lo,f_hi);
                
clr_anno=[0 0 0];
clr_mike=[1 0 0 ; ...
          0 0.7 0 ; ...
          0 0 1 ; ...
          0 0.8 0.8 ];
figure_objective_map(x_grid,y_grid,[], ...
                     @jet, ...
                     [], ...
                     title_str, ...
                     'RMSE (mV)', ...
                     clr_mike, ...
                     clr_anno, ...
                     r_est,[], ...
                     R,r_head,r_tail);
n_mice=size(r_head,2);
for i=1:n_mice
  r_mouse_ellipse=polygon_from_ellipse(r_center(:,i),a_vec(:,i),b(i));
  line(100*r_mouse_ellipse(1,:),100*r_mouse_ellipse(2,:),'color','k');
end
% line(100*r_est_in_body(1),100*r_est_in_body(2),0, ...
%      'marker','o','linestyle','none','color','w', ...
%      'markersize',6);
