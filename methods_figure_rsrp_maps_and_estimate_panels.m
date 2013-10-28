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

cache_file_name='methods_figure_1_rsrp_maps_and_estimate_panels_cache.mat';
if ~exist(cache_file_name,'file') ,
%if true ,
  % identifying info for the segment
  date_str='06132012';
  letter_str='D';
  i_segment=51;  % this was voc84 in the old-style

  % directories where to find stuff
  %base_dir_name='~/egnor_stuff/ssl_vocal_structure_bizarro';
  % base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_vocal_structure';
  base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_sys_test';
  data_analysis_dir_name='Data_analysis10';

  % call r_est_from_voc_indicators_and_ancillary() for each voc, collect
  % all the results in r_est_blob_per_voc_per_trial
  %[r_est_blob_per_voc_per_trial,per_trial_ancillary]= ... 
  %   scatter_single_segment(base_dir_name, ...
  %                          data_analysis_dir_name, ...
  %                          date_str, ...
  %                          letter_str, ...
  %                          i_segment_want, ...
  %                          n_vocs_per_trial_max, ...
  %                          n_vocs_per_job_max, ...
  %                          use_cluster, ...
  %                          @r_est_from_voc_indicators_and_ancillary_for_a_few_vocs, ...
  %                          args, ...
  %                          verbosity);
  trial_overhead=ssl_trial_overhead_packaged(base_dir_name,...
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

  % % save everything           
  % load('analysis_for_methods_figure_1.mat');

  % unpack the return blob
  i_snippet_pretty=16;
  field_names=fieldnames(r_est_blobs);
  for i=1:length(field_names)
    if isequal(field_names{i},'rsrp_per_pair_grid')
      rsrp_per_pair_grid_pretty_snippet=r_est_blobs(i_snippet_pretty).rsrp_per_pair_grid;
    else
      eval(sprintf('%s_all_snippets={r_est_blobs(:).%s}'';',field_names{i},field_names{i}));
    end
  end
  clear r_est_blobs;

  % transform things from cell arrays what can
  n_snippets=length(tf_rect_name_all_snippets);  
  rsrp_grid_all_snippets=cell2mat(reshape(rsrp_grid_all_snippets,[1 1 n_snippets]));
  %rsrp_per_pair_grid_all_snippets=cell2mat(reshape(rsrp_per_pair_grid_all_snippets,[1 1 1 n_snippets]));
  r_est_all_snippets=cell2mat(reshape(r_est_all_snippets,[1 n_snippets]));
  r_head_from_video_all_snippets=cell2mat(reshape(r_head_from_video_all_snippets,[1 n_snippets]));
  r_tail_from_video_all_snippets=cell2mat(reshape(r_tail_from_video_all_snippets,[1 n_snippets]));

  % unpack the trial overhead
  R=trial_overhead.R;
  Temp=trial_overhead.Temp;
  dx=trial_overhead.dx;
  x_grid=trial_overhead.x_grid;
  y_grid=trial_overhead.y_grid;
  in_cage=trial_overhead.in_cage;
  fs=trial_overhead.fs;
  clear trial_overhead;
  save(cache_file_name);
else
  load(cache_file_name);
end

% get dims out
N=size(v_all_snippets{1},1);  % all same length
K=size(R,2);
dt=1/fs;  % s
n_snippets=length(tf_rect_name_all_snippets)
n_pairs=nchoosek(K,2);

% figure setup
% set the default to be that figures print the same size as on-screen
set(0,'DefaultFigurePaperPositionMode','auto');

% make it so that it doesn't change the figure or axes background of
% printed figures
set(0,'DefaultFigureInvertHardCopy','off');

% set up so the default is not to change the figure axis limits and ticks
% when printing
newPrintTemplate=printtemplate;
newPrintTemplate.AxesFreezeTicks=1;
newPrintTemplate.AxesFreezeLimits=1;
newPrintTemplate.DriverColor=1;
set(0,'DefaultFigurePrintTemplate',newPrintTemplate); 
clear newPrintTemplate

% set mic colors
% clr_mike=[1 0 0 ; ...
%           0 0.7 0 ; ...
%           0 0 1 ; ...
%           0 0.8 0.8 ];
clr_mike=zeros(4,3);

do_draw_mic_labels=true;                 
for i_snippet=i_snippet_pretty
  % unpack stuff for this snippet
  tf_rect_name=tf_rect_name_all_snippets{i_snippet};
  v=v_all_snippets{i_snippet};  %#ok
  V_filt=V_filt_all_snippets{i_snippet};  %#ok
  a=a_all_snippets{i_snippet};  %#ok
  rsrp_grid=rsrp_grid_all_snippets(:,:,i_snippet);
  %rsrp_per_pair_grid=rsrp_per_pair_grid_all_snippets(:,:,:,i_snippet);
  r_est=r_est_all_snippets(:,i_snippet);
  r_head_from_video=r_head_from_video_all_snippets(:,i_snippet);
  r_tail_from_video=r_tail_from_video_all_snippets(:,i_snippet);
  
  % useful stuff        
  N=size(v,1);
  t=dt*(0:(N-1))';  % s

%   % get the filtered signals back from the f domain
%   v_filt=real(ifft(V_filt));
% 
%   % set up the figure and place all the axes                                   
%   w_fig=2.5; % in
%   h_fig=3; % in
%   n_row=K;
%   n_col=1;
%   w_axes=1.5;  % in
%   h_axes=0.5;  % in
%   w_space=1;  % in (not used)
%   h_space=0;  % in                              
%   [figure_handle,subplot_handles]= ...
%     layout_axes_grid(w_fig,h_fig,...
%                      n_row,n_col,...
%                      w_axes,h_axes,...
%                      w_space,h_space);
%   set(figure_handle,'color','w');                               
% 
%   % plot the filtered clips with raw in background
%   white_fraction=0.75;
%   %figure_handle=figure('color','w');
%   %set_figure_size_explicit(figure_handle,[3 6]);
%   for k=1:K
%     subplot_handle=subplot_handles(k);
%     axes(subplot_handle);  %#ok
%     plot(1000*t,1000*v(:,k)     ,'color',(1-white_fraction)*clr_mike(k,:)+white_fraction*[1 1 1]);
%     hold on
%     plot(1000*t,1000*v_filt(:,k),'color',clr_mike(k,:));
%     hold off
%     set(subplot_handle,'fontsize',7);
%     ylim(ylim_tight(1000*v(:,k)));
%     %ylabel(sprintf('Mic %d (mV)',k),'fontsize',7);
%     if k~=K ,
%       set(subplot_handle,'xticklabel',{});
%       set(subplot_handle,'yticklabel',{});
%     else
%       set(subplot_handle,'yAxisLocation','right');
%     end
%   end
%   %xlabel('Time (ms)','fontsize',7);
%   ylim_all_same();
%   tl(0,1000*t(end));
% 
%   % print the per-mic RMS amplitudes
%   r_in_mV=1000*a  %#ok

  % % RMSEs are easier to interpret
  % rmse_grid=sqrt(mse_grid);
  % rmse_min=sqrt(mse_min);
  % rmse_crit=sqrt(mse_crit);
  % rmse_body=sqrt(mse_body);
  % rms_total=sqrt(ms_total);

  % % display some of those
  % rmse_min_disp=1e3*rmse_min
  % rms_total_disp=1e3*rms_total

  % % extrapolate a mouse ellipse
  % r_center=(r_head_from_video+r_tail)/2;
  % a_vec=r_head_from_video-r_center;  % vector
  % b=norm(a_vec)/3;  % scalar, and a guess at the half-width of the mouse

  % Compute the CR borders from the objective function and critical value
  %in_cr_grid=(mse_grid<=mse_crit);
  %r_cr_bounds_ij=bwboundaries(in_cr_grid);
  %r_cr_bounds=path_xy_from_path_ij(r_cr_bounds_ij,x_grid,y_grid);

  % make a figure of the objective funtion, estimate
  rsrp_abs_max=max(max(abs(rsrp_grid)));
  title_str=sprintf('%s %s %s', ...
                    date_str,letter_str,tf_rect_name);
  colorbar_label_str='RSRP/sample (mV^2)';         
  clr_anno=[0 0 0];

  w_fig=4.5; % in
  h_fig=3; % in
  w_axes=1.8;  % in
  h_axes=1.8;  % in
  w_colorbar=0.1;  % in
  w_colorbar_spacer=0.05;  % in

  fig_h=figure('color','w');
  set_figure_size_explicit(fig_h,[w_fig h_fig]);
  axes_h=axes('parent',fig_h);  %#ok
  set(axes_h,'fontsize',7);
  set_axes_size_fixed_center_explicit(axes_h,[w_axes h_axes])

  place_objective_map_in_axes(axes_h, ...
                              x_grid,y_grid,10^6*rsrp_grid/N, ...
                              @bipolar_red_white_blue, ...
                              10^6*rsrp_abs_max/N*[-1 +1], ...
                              title_str, ...
                              clr_mike, ...
                              clr_anno, ...
                              r_est,[], ...
                              R,r_head_from_video,r_tail_from_video, ...
                              do_draw_mic_labels);
  title(axes_h,title_str,'interpreter','none','fontsize',7);
  h=xlabel(axes_h,'');  delete(h);
  h=ylabel(axes_h,'');  delete(h);
  set(axes_h,'ytick',[]);
  set(axes_h,'xtick',[40 60]);

  axes_cb_h=add_colorbar(axes_h,w_colorbar,w_colorbar_spacer);
  set(axes_cb_h,'fontsize',7);
  ylabel(axes_cb_h,colorbar_label_str);



  % normalize the per-pair rsrp images by the RMS amplitude of each signal
  [M,iMicFirst,iMicSecond]=mixing_matrix_from_n_mics(K);
  rsrp_per_pair_grid_normed=zeros(size(rsrp_per_pair_grid_pretty_snippet));
  for i_pair=1:n_pairs
    keep=logical(abs(M(i_pair,:)));
    norming_factor=sqrt(prod(a(keep).^2))*N;
    rsrp_per_pair_grid_normed(:,:,i_pair)=rsrp_per_pair_grid_pretty_snippet(:,:,i_pair)/norming_factor;
  end

  % plot the rsrp image for each pair in a single figure
  max_abs=max(max(max(abs(rsrp_per_pair_grid_normed))));
  [fig_h,axes_hs]= ...
    figure_objective_map_per_pair_grid(x_grid,y_grid,rsrp_per_pair_grid_normed, ...
                                       @bipolar_red_white_blue, ...
                                       max_abs*[-1 +1], ...
                                       clr_mike, ...
                                       clr_anno, ...
                                       r_est,[], ...
                                       R,r_head_from_video,r_tail_from_video);

  %close all;
end  % for i_snippets=...
















% need to do outlier filtering on r_est here
[is_outlier,~,r_est_trans,Covariance_matrix] = kur_rce(r_est_all_snippets',1);
is_outlier=logical(is_outlier);
indices_of_outliers=find(is_outlier)
r_est=r_est_trans';
n_outliers=sum(is_outlier)

% filter out the outliers
is_keeper=~is_outlier;
n_keepers=sum(is_keeper);
rsrp_grid_all_keepers=rsrp_grid_all_snippets(:,:,is_keeper);
r_est_all_keepers=r_est_all_snippets(:,is_keeper);
r_est_all_outliers=r_est_all_snippets(:,is_outlier);

% take the mean of the maps for all the non-outliers
rsrp_grid=mean(rsrp_grid_all_keepers,3);
rsrp_abs_max=max(max(abs(rsrp_grid)));
r_head_from_video=mean(r_head_from_video_all_snippets,2);
r_tail_from_video=mean(r_tail_from_video_all_snippets,2);

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

% caclulate the density at the real+fake mice, and the posterior
% probability
r_head_from_video_real_and_fake=[r_head_from_video r_head_from_video_fake];
r_tail_from_video_real_and_fake=[r_tail_from_video r_tail_from_video_fake];
r_chest_from_video_real_and_fake= ...
  (3/4)*r_head_from_video_real_and_fake + ...
  (1/4)*r_tail_from_video_real_and_fake ;
p_chest_real_and_fake=mvnpdf(r_chest_from_video_real_and_fake',r_est',Covariance_matrix)  % density
P_posterior_chest_from_video_real_and_fake= ...
  p_chest_real_and_fake/sum(p_chest_real_and_fake)  % posterior probability

% calculate the pdf at all grid points, so we can show it on the next
% figure
x_grid_serial=x_grid(:);
y_grid_serial=y_grid(:);
r_grid_serial=[x_grid_serial';y_grid_serial'];
n_grid_pts=size(r_grid_serial,2);
pdf_grid_serial=mvnpdf(r_grid_serial',r_est',Covariance_matrix);  % density, units: 1/(m^2)
pdf_grid=reshape(pdf_grid_serial,size(x_grid));
pdf_grid_max=max(max(pdf_grid));

%
% make a figure of the pdf funtion, and the position estimates
%
% title_str=sprintf('%s %s seg %d', ...
%                   date_str,letter_str,i_segment);
title_str='';                
colorbar_label_str='Probability density (1/m^2)';         
clr_anno=[0 0 0];

clr_mouse=[0 0 145]/255;  % m1 color: dark blue
clr_mice_fake=0.8*[   0  255  255 ; ...
                   255   69    0 ; ...
                   255  215    0 ]/255;
do_draw_mic_labels=true;                 

w_fig=4.5; % in
h_fig=3; % in
w_axes=1.8;  % in
h_axes=1.8;  % in
w_colorbar=0.1;  % in
w_colorbar_spacer=0.05;  % in

fig_h=figure('color','w');
set_figure_size_explicit(fig_h,[w_fig h_fig]);
axes_h=axes('parent',fig_h,'box','on','ydir','reverse');
set(axes_h,'fontsize',7);
set_axes_size_fixed_center_explicit(axes_h,[w_axes h_axes])

place_objective_map_in_axes(axes_h, ...
                            x_grid,y_grid,pdf_grid, ...
                            @(n)(flipud(0.7+0.3*gray(n))), ...
                            [0 pdf_grid_max], ...
                            title_str, ...
                            clr_mike, ...
                            clr_anno, ...
                            r_est,[], ...
                            R,[],[], ...
                            do_draw_mic_labels);
z_mice=1;                          
draw_mice_given_head_and_tail(axes_h,100*r_head_from_video,100*r_tail_from_video,z_mice,clr_mouse);
draw_mice_given_head_and_tail(axes_h,100*r_head_from_video_fake,100*r_tail_from_video_fake,z_mice,clr_mice_fake);
line('parent',axes_h, ...
     'xdata',100*r_est_all_keepers(1,:), ...
     'ydata',100*r_est_all_keepers(2,:), ...
     'zdata',ones(1,n_keepers), ...
     'marker','+', ...
     'linestyle','none', ...
     'color','k');
line('parent',axes_h, ...
     'xdata',100*r_est(1,:), ...
     'ydata',100*r_est(2,:), ...
     'zdata',1, ...
     'marker','.', ...
     'markersize',6*3, ...
     'linestyle','none', ...
     'color','k');
line('parent',axes_h, ...
     'xdata',100*r_est_all_outliers(1,:), ...
     'ydata',100*r_est_all_outliers(2,:), ...
     'zdata',ones(1,n_outliers), ...
     'marker','o', ...
     'markersize',3, ...
     'linestyle','none', ...
     'color','k');
%hold on;
%h_error_ellipse=error_ellipse('mu',100*r_est,'C',100^2*Covariance_matrix,'conf',0.95);
%hold off;
%set(h_error_ellipse,'color','k');
%title(axes_h,title_str,'interpreter','none','fontsize',7);
h=xlabel(axes_h,'');  delete(h);
h=ylabel(axes_h,'');  delete(h);
set(axes_h,'ytick',[]);
set(axes_h,'xtick',[40 60]);

axes_cb_h=add_colorbar(axes_h,w_colorbar,w_colorbar_spacer);
set(axes_cb_h,'fontsize',7);
ylabel(axes_cb_h,colorbar_label_str);


