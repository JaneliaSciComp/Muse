% identifying info for the segment, snippet
date_str='06132012';
letter_str='D';
base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_sys_test';
data_analysis_dir_name='Data_analysis10';
fs=450450;  % Hz, happen to know this a priori
fs_video=29;  % Hz
clr_mike=[1 0 0 ; ...
          0 0.7 0 ; ...
          0 0 1 ; ...
          0 0.8 0.8 ];
associated_video_frame_method = 'close'; %options are begin or close
example_segment_id_string='Voc130';
%i_segment=51; 
i_snippet=16;  % index of the example snippet
%are_positions_on_disk_in_old_style_coords=true;  % uses Josh's coord convention from the pre-Motr days

%session_base_name=sprintf('Test_%s_1',letter_str);
%fs_video=29;  % Hz, ditto

% figure out the name of the stupid ax output file
exp_dir_name=fullfile(base_dir_name, ...
                      sprintf('sys_test_%s',date_str));
demuxed_data_dir_name = fullfile(exp_dir_name,'demux');
ax_output_parent_dir_name=fullfile(demuxed_data_dir_name,'no_merge_only_har');
ax_output_dir_name_pattern = sprintf('*_%s_*',letter_str);
ax_output_parent_dir_listing_struct = dir(fullfile(ax_output_parent_dir_name,ax_output_dir_name_pattern));
ax_output_dir_name=ax_output_parent_dir_listing_struct.name;
ax_output_mat_file_name = sprintf('Test_%s_1_voc_list_no_merge_har.mat',letter_str);
ax_output_mat_file_name_abs=fullfile(ax_output_parent_dir_name,ax_output_dir_name,ax_output_mat_file_name);

% read the video frame pulse data
session_base_name=sprintf('Test_%s_1',letter_str);
yn_load_time_stamps = 'y';
video_pulse_start_ts = ...
  fn_video_pulse_start_ts(demuxed_data_dir_name, ...
                          exp_dir_name, ...
                          session_base_name, ...
                          session_base_name, ...
                          yn_load_time_stamps, ...
                          fs, ...
                          fs_video);

% load the raw ax output, and the figure out the frame index that goes with 
% each snippet
mouse_from_ax=load_ax_segments_and_append_frame_number(ax_output_mat_file_name_abs,video_pulse_start_ts,associated_video_frame_method);

% extract the segment we want
%example_segment_name='Voc130';
is_example_segment= strcmp(example_segment_id_string,{mouse_from_ax.syl_name});
mouse_example_segment=mouse_from_ax(is_example_segment);

% extract the corresponding video frame
i_frame_example=mouse_example_segment.frame_number;
video_file_name=fullfile(exp_dir_name, ...
                         sprintf('Test_%s_1.seq',letter_str));
video_info=fnReadSeqInfo_jpn(video_file_name);
example_frame=fnReadFrameFromSeq(video_info,i_frame_example);
[frame_height_in_pels,frame_width_in_pels]=size(example_frame);
%r_frame_upper_left_pel_center_pels=[1/2 size(example_frame,1)-1/2]';  
  % (non-image-style) cartesian coords of the center of the upper-left pel
%r_frame_lower_right_pel_center_pels=[size(example_frame,2)-1/2 1/2]';  
  % (non-image-style) cartesian coords of the center of the lower-right pel
% this puts one corner of the image at the origin

% get the syl_names, etc for all segments
[snippet_id_string_all,i_start_all,i_end_all,f_lo_all,f_hi_all, ...
 r_head_all,r_tail_all,R,Temp, ...
 dx,x_grid,y_grid,in_cage,r_corners]= ...
  ssl_trial_overhead(base_dir_name, ...
                     data_analysis_dir_name, ...
                     date_str, ...
                     letter_str);
n_mics=size(R,2);

%xl_example_frame_pels=[1 size(example_frame,2)];
%yl_example_frame_pels=[1 size(example_frame,1)];

% load the meters/pixel scaling factor
meters_2_pixels_file_name= ...
  fullfile(exp_dir_name, ...
           'meters_2_pixels.mat');
s=load(meters_2_pixels_file_name);
meters_per_pixel=s.meters_2_pixels;

% convert the image limits to meters
%xl_example_frame_m=meters_per_pixel*xl_example_frame_pels;  % meters
%yl_example_frame_m=meters_per_pixel*yl_example_frame_pels;  % meters
%r_frame_upper_left_pel_center=meters_per_pixel*r_frame_upper_left_pel_center_pels;  % meters
%r_frame_lower_right_pel_center=meters_per_pixel*r_frame_lower_right_pel_center_pels;

% get the r_head and r_tail for that segment
is_example_snippet= strncmp('Voc000051',snippet_id_string_all,9);
r_head_example_segment=r_head_all(:,is_example_snippet);
r_tail_example_segment=r_tail_all(:,is_example_snippet);

% % x and y are swapped for head and tail, need to sort out
% % qucik fix:
% r_head_example_segment=flipud(r_head_example_segment)
% r_tail_example_segment=flipud(r_tail_example_segment)

% make the figure
w_fig=4.5; % in
h_fig=3; % in
w_axes=1;  % in
h_axes=1;  % in
%w_colorbar=0.1;  % in
%w_colorbar_spacer=0.05;  % in

fig_h=figure('color','w');
set_figure_size_explicit(fig_h,[w_fig h_fig]);
%axes_h=axes('parent',fig_h);  %#ok

%fig_h=figure('color','w');
colormap(fig_h,gray(256));
axes_h=axes('parent',fig_h, ...
            'dataaspectratio',[1 1 1], ...
            'clim',[0 160], ...
            'ydir','reverse', ...
            'box','on', ...
            'layer','top');
%            'clim',[0 max(max(example_frame))], ...
set(axes_h,'fontsize',7);
set_axes_size_fixed_center_explicit(axes_h,[w_axes h_axes])
image('parent', axes_h, ...
      'xdata',100*meters_per_pixel*[1 frame_width_in_pels ], ...
      'ydata',100*meters_per_pixel*[1 frame_height_in_pels], ...
      'cdata',example_frame, ...
      'cdatamapping','scaled');
for i_mic=1:n_mics
  line('parent',axes_h, ...
       'xdata',100*R(1,i_mic), ...
       'ydata',100*R(2,i_mic), ...
       'zdata',0, ...
       'marker','.', ...
       'linestyle','none', ...
       'color',clr_mike(i_mic,:), ...
       'markersize',18);       
end
% line('parent',axes_h, ...
%      'xdata',100*[r_corners(1,:) r_corners(1,1)], ...
%      'ydata',100*[r_corners(2,:) r_corners(2,1)], ...
%      'zdata',[0 0 0 0 0], ...
%      'marker','none', ...
%      'color','k');       
%mouse_line_handles= ...
%  draw_mice_given_head_and_tail(axes_h,r_head_example_segment,r_tail_example_segment,1);
% set(mouse_line_handles,'color','w');   

%xlim(axes_h,100*[0 r_frame_lower_right_pel_center(1)+0.5*meters_per_pixel]);
%ylim(axes_h,100*[0 r_frame_upper_left_pel_center(2)+0.5*meters_per_pixel]);
xlim(axes_h,100*[min(R(1,:)) max(R(1,:))]);
ylim(axes_h,100*[min(R(2,:)) max(R(2,:))]);

set(axes_h,'ytick',[]);
set(axes_h,'xtick',[40 60]);
   
   
   
% % extract the example snippet
% mouse_example_snippet=mouse_example_segment(i_snippet);



% % plot the example frame, and everything else
% clr_anno=[0 0 0];
% 
% w_fig=4.5; % in
% h_fig=3; % in
% w_axes=1.8;  % in
% h_axes=1.8;  % in
% w_colorbar=0.1;  % in
% w_colorbar_spacer=0.05;  % in
% 
% fig_h=figure('color','w');
% set_figure_size_explicit(fig_h,[w_fig h_fig]);
% axes_h=axes('parent',fig_h);  %#ok
% set(axes_h,'fontsize',7);
% set_axes_size_fixed_center_explicit(axes_h,[w_axes h_axes])
% 
% 
% 
% 
% 
% place_objective_map_in_axes(axes_h, ...
%                             x_grid,y_grid,10^6*rsrp_grid/N, ...
%                             @bipolar_red_white_blue, ...
%                             10^6*rsrp_abs_max/N*[-1 +1], ...
%                             title_str, ...
%                             clr_mike, ...
%                             clr_anno, ...
%                             r_est,[], ...
%                             R,r_head,r_tail);
% 
% 
% 
% 
% 
%                    
%                    
%                    % directories where to find stuff
%   %base_dir_name='~/egnor_stuff/ssl_vocal_structure_bizarro';
%   % base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_vocal_structure';
% 
%   % call r_est_from_voc_indicators_and_ancillary() for each voc, collect
%   % all the results in r_est_blob_per_voc_per_trial
%   %[r_est_blob_per_voc_per_trial,per_trial_ancillary]= ... 
%   %   scatter_single_segment(base_dir_name, ...
%   %                          data_analysis_dir_name, ...
%   %                          date_str, ...
%   %                          letter_str, ...
%   %                          i_segment_want, ...
%   %                          n_vocs_per_trial_max, ...
%   %                          n_vocs_per_job_max, ...
%   %                          use_cluster, ...
%   %                          @r_est_from_voc_indicators_and_ancillary_for_a_few_vocs, ...
%   %                          args, ...
%   %                          verbosity);
%   [r_est_blobs,trial_overhead]= ...
%     r_est_from_segment_indicators(base_dir_name,...
%                                   data_analysis_dir_name, ...
%                                   date_str, ...
%                                   letter_str, ...
%                                   i_segment, ...
%                                   args, ...
%                                   verbosity);
% 
%   % % save everything           
%   % load('analysis_for_methods_figure_1.mat');
% 
%   % unpack the return blob
%   i_snippet=16;
%   field_names=fieldnames(r_est_blobs);
%   for i=1:length(field_names)
%     if isequal(field_names{i},'rsrp_per_pair_grid')
%       rsrp_per_pair_grid_pretty_snippet=r_est_blobs(i_snippet).rsrp_per_pair_grid;
%     else
%       eval(sprintf('%s_all_snippets={r_est_blobs(:).%s}'';',field_names{i},field_names{i}));
%     end
%   end
%   clear r_est_blobs;
% 
%   % transform things from cell arrays what can
%   n_snippets=length(syl_name_all_snippets);  
%   rsrp_grid_all_snippets=cell2mat(reshape(rsrp_grid_all_snippets,[1 1 n_snippets]));
%   %rsrp_per_pair_grid_all_snippets=cell2mat(reshape(rsrp_per_pair_grid_all_snippets,[1 1 1 n_snippets]));
%   r_est_all_snippets=cell2mat(reshape(r_est_all_snippets,[1 n_snippets]));
%   r_head_all_snippets=cell2mat(reshape(r_head_all_snippets,[1 n_snippets]));
%   r_tail_all_snippets=cell2mat(reshape(r_tail_all_snippets,[1 n_snippets]));
% 
%   % unpack the trial overhead
%   field_name=fieldnames(trial_overhead);
%   for i=1:length(field_name)
%     eval(sprintf('%s=trial_overhead.%s;',field_name{i},field_name{i}));
%   end
%   clear trial_overhead;
%   save(cache_file_name);
% else
%   load(cache_file_name);
% end
% 
% % get dims out
% N=size(v_all_snippets{1},1);  % all same length
% K=size(R,2);
% dt=1/fs;  % s
% n_snippets=length(syl_name_all_snippets)
% n_pairs=nchoosek(K,2);
% 
% % need to do outlier filtering on r_est here
% [is_outlier,~,r_est_trans,Covariance_matrix] = kur_rce(r_est_all_snippets',1);
% is_outlier=logical(is_outlier);
% indices_of_outliers=find(is_outlier)
% r_est=r_est_trans';
% n_outliers=sum(is_outlier)
% 
% % filter out the outliers
% is_keeper=~is_outlier;
% n_keepers=sum(is_keeper);
% rsrp_grid_all_keepers=rsrp_grid_all_snippets(:,:,is_keeper);
% r_est_all_keepers=r_est_all_snippets(:,is_keeper);
% r_est_all_outliers=r_est_all_snippets(:,is_outlier);
% 
% % take the mean of the maps for all the non-outliers
% rsrp_grid=mean(rsrp_grid_all_keepers,3);
% r_head=mean(r_head_all_snippets,2);
% r_tail=mean(r_tail_all_snippets,2);
% 
% % make up some mouse locations
% n_fake_mice=3;
% [r_head_fake,r_tail_fake]= ...
%   random_mouse_locations(R,r_head,r_tail,n_fake_mice);
% 
% % caclulate the density at the real+fake mice, and the posterior
% % probability
% r_head_real_and_fake=[r_head r_head_fake];
% p_head_real_and_fake=mvnpdf(r_head_real_and_fake',r_est',Covariance_matrix)  % density
% P_posterior_head_real_and_fake= ...
%   p_head_real_and_fake/sum(p_head_real_and_fake)  % posterior probability
% 
% % figure setup
% % set the default to be that figures print the same size as on-screen
% set(0,'DefaultFigurePaperPositionMode','auto');
% 
% % make it so that it doesn't change the figure or axes background of
% % printed figures
% set(0,'DefaultFigureInvertHardCopy','off');
% 
% % set up so the default is not to change the figure axis limits and ticks
% % when printing
% newPrintTemplate=printtemplate;
% newPrintTemplate.AxesFreezeTicks=1;
% newPrintTemplate.AxesFreezeLimits=1;
% newPrintTemplate.DriverColor=1;
% set(0,'DefaultFigurePrintTemplate',newPrintTemplate); 
% clear newPrintTemplate
% 
% % set mic colors
% clr_mike=[1 0 0 ; ...
%           0 0.7 0 ; ...
%           0 0 1 ; ...
%           0 0.8 0.8 ];
% 
% %
% % make a figure of the mean objective funtion, and the position estimates
% %
% rsrp_abs_max=max(max(abs(rsrp_grid)));
% title_str=sprintf('%s %s seg %d', ...
%                   date_str,letter_str,i_segment);
% colorbar_label_str='RSRP/sample (mV^2)';         
% clr_anno=[0 0 0];
% 
% w_fig=4.5; % in
% h_fig=3; % in
% w_axes=1.8;  % in
% h_axes=1.8;  % in
% w_colorbar=0.1;  % in
% w_colorbar_spacer=0.05;  % in
% 
% fig_h=figure('color','w');
% set_figure_size_explicit(fig_h,[w_fig h_fig]);
% axes_h=axes('parent',fig_h);
% set(axes_h,'fontsize',7);
% set_axes_size_fixed_center_explicit(axes_h,[w_axes h_axes])
% 
% place_objective_map_in_axes(axes_h, ...
%                             x_grid,y_grid,10^6*rsrp_grid/N, ...
%                             @bipolar_red_white_blue, ...
%                             10^6*rsrp_abs_max/N*[-1 +1], ...
%                             title_str, ...
%                             clr_mike, ...
%                             clr_anno, ...
%                             r_est,[], ...
%                             R,r_head,r_tail);
% fake_mice_handles=draw_mice_given_head_and_tail(axes_h,r_head_fake,r_tail_fake,1);
% line('parent',axes_h, ...
%      'xdata',100*r_est_all_keepers(1,:), ...
%      'ydata',100*r_est_all_keepers(2,:), ...
%      'zdata',ones(1,n_keepers), ...
%      'marker','+', ...
%      'linestyle','none', ...
%      'color','k');
% line('parent',axes_h, ...
%      'xdata',100*r_est_all_outliers(1,:), ...
%      'ydata',100*r_est_all_outliers(2,:), ...
%      'zdata',ones(1,n_outliers), ...
%      'marker','o', ...
%      'linestyle','none', ...
%      'color','k');
% hold on;
% h_error_ellipse=error_ellipse('mu',100*r_est,'C',100^2*Covariance_matrix,'conf',0.95);
% hold off;
% set(h_error_ellipse,'color','k');
% title(axes_h,title_str,'interpreter','none','fontsize',7);
% h=xlabel(axes_h,'');  delete(h);
% h=ylabel(axes_h,'');  delete(h);
% set(axes_h,'ytick',[]);
% set(axes_h,'xtick',[40 60]);
% 
% axes_cb_h=add_colorbar(axes_h,w_colorbar,w_colorbar_spacer);
% set(axes_cb_h,'fontsize',7);
% ylabel(axes_cb_h,colorbar_label_str);
% 
% for i_snippet=i_snippet
%   % unpack stuff for this snippet
%   syl_name=syl_name_all_snippets{i_snippet};
%   v=v_all_snippets{i_snippet};  %#ok
%   V_filt=V_filt_all_snippets{i_snippet};  %#ok
%   a=a_all_snippets{i_snippet};  %#ok
%   rsrp_grid=rsrp_grid_all_snippets(:,:,i_snippet);
%   %rsrp_per_pair_grid=rsrp_per_pair_grid_all_snippets(:,:,:,i_snippet);
%   r_est=r_est_all_snippets(:,i_snippet);
%   r_head=r_head_all_snippets(:,i_snippet);
%   r_tail=r_tail_all_snippets(:,i_snippet);
%   
%   % useful stuff        
%   N=size(v,1);
%   t=dt*(0:(N-1))';  % s
% 
% %   % get the filtered signals back from the f domain
% %   v_filt=real(ifft(V_filt));
% % 
% %   % set up the figure and place all the axes                                   
% %   w_fig=2.5; % in
% %   h_fig=3; % in
% %   n_row=K;
% %   n_col=1;
% %   w_axes=1.5;  % in
% %   h_axes=0.5;  % in
% %   w_space=1;  % in (not used)
% %   h_space=0;  % in                              
% %   [figure_handle,subplot_handles]= ...
% %     layout_axes_grid(w_fig,h_fig,...
% %                      n_row,n_col,...
% %                      w_axes,h_axes,...
% %                      w_space,h_space);
% %   set(figure_handle,'color','w');                               
% % 
% %   % plot the filtered clips with raw in background
% %   white_fraction=0.75;
% %   %figure_handle=figure('color','w');
% %   %set_figure_size_explicit(figure_handle,[3 6]);
% %   for k=1:K
% %     subplot_handle=subplot_handles(k);
% %     axes(subplot_handle);  %#ok
% %     plot(1000*t,1000*v(:,k)     ,'color',(1-white_fraction)*clr_mike(k,:)+white_fraction*[1 1 1]);
% %     hold on
% %     plot(1000*t,1000*v_filt(:,k),'color',clr_mike(k,:));
% %     hold off
% %     set(subplot_handle,'fontsize',7);
% %     ylim(ylim_tight(1000*v(:,k)));
% %     %ylabel(sprintf('Mic %d (mV)',k),'fontsize',7);
% %     if k~=K ,
% %       set(subplot_handle,'xticklabel',{});
% %       set(subplot_handle,'yticklabel',{});
% %     else
% %       set(subplot_handle,'yAxisLocation','right');
% %     end
% %   end
% %   %xlabel('Time (ms)','fontsize',7);
% %   ylim_all_same();
% %   tl(0,1000*t(end));
% % 
% %   % print the per-mic RMS amplitudes
% %   r_in_mV=1000*a  %#ok
% 
%   % % RMSEs are easier to interpret
%   % rmse_grid=sqrt(mse_grid);
%   % rmse_min=sqrt(mse_min);
%   % rmse_crit=sqrt(mse_crit);
%   % rmse_body=sqrt(mse_body);
%   % rms_total=sqrt(ms_total);
% 
%   % % display some of those
%   % rmse_min_disp=1e3*rmse_min
%   % rms_total_disp=1e3*rms_total
% 
%   % % extrapolate a mouse ellipse
%   % r_center=(r_head+r_tail)/2;
%   % a_vec=r_head-r_center;  % vector
%   % b=norm(a_vec)/3;  % scalar, and a guess at the half-width of the mouse
% 
%   % Compute the CR borders from the objective function and critical value
%   %in_cr_grid=(mse_grid<=mse_crit);
%   %r_cr_bounds_ij=bwboundaries(in_cr_grid);
%   %r_cr_bounds=path_xy_from_path_ij(r_cr_bounds_ij,x_grid,y_grid);
% 
%   % make a figure of the objective funtion, estimate
%   rsrp_abs_max=max(max(abs(rsrp_grid)));
%   title_str=sprintf('%s %s %s', ...
%                     date_str,letter_str,syl_name);
%   colorbar_label_str='RSRP/sample (mV^2)';         
%   clr_anno=[0 0 0];
% 
%   w_fig=4.5; % in
%   h_fig=3; % in
%   w_axes=1.8;  % in
%   h_axes=1.8;  % in
%   w_colorbar=0.1;  % in
%   w_colorbar_spacer=0.05;  % in
% 
%   fig_h=figure('color','w');
%   set_figure_size_explicit(fig_h,[w_fig h_fig]);
%   axes_h=axes('parent',fig_h);  %#ok
%   set(axes_h,'fontsize',7);
%   set_axes_size_fixed_center_explicit(axes_h,[w_axes h_axes])
% 
%   place_objective_map_in_axes(axes_h, ...
%                               x_grid,y_grid,10^6*rsrp_grid/N, ...
%                               @bipolar_red_white_blue, ...
%                               10^6*rsrp_abs_max/N*[-1 +1], ...
%                               title_str, ...
%                               clr_mike, ...
%                               clr_anno, ...
%                               r_est,[], ...
%                               R,r_head,r_tail);
%   title(axes_h,title_str,'interpreter','none','fontsize',7);
%   h=xlabel(axes_h,'');  delete(h);
%   h=ylabel(axes_h,'');  delete(h);
%   set(axes_h,'ytick',[]);
%   set(axes_h,'xtick',[40 60]);
% 
%   axes_cb_h=add_colorbar(axes_h,w_colorbar,w_colorbar_spacer);
%   set(axes_cb_h,'fontsize',7);
%   ylabel(axes_cb_h,colorbar_label_str);
% 
% 
% 
%   % normalize the per-pair rsrp images by the RMS amplitude of each signal
%   [M,iMicFirst,iMicSecond]=mixing_matrix_from_n_mics(K);
%   rsrp_per_pair_grid_normed=zeros(size(rsrp_per_pair_grid_pretty_snippet));
%   for i_pair=1:n_pairs
%     keep=logical(abs(M(i_pair,:)));
%     norming_factor=sqrt(prod(a(keep).^2))*N;
%     rsrp_per_pair_grid_normed(:,:,i_pair)=rsrp_per_pair_grid_pretty_snippet(:,:,i_pair)/norming_factor;
%   end
% 
%   % plot the rsrp image for each pair in a single figure
%   max_abs=max(max(max(abs(rsrp_per_pair_grid_normed))));
%   [fig_h,axes_hs]= ...
%     figure_objective_map_per_pair_grid(x_grid,y_grid,rsrp_per_pair_grid_normed, ...
%                                        @bipolar_red_white_blue, ...
%                                        max_abs*[-1 +1], ...
%                                        clr_mike, ...
%                                        clr_anno, ...
%                                        r_est,[], ...
%                                        R,r_head,r_tail);
% 
%   %close all;
% end  % for i_snippets=...
% 
% 
