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

% identifying info for the voc

% good one
date_str='06132012';
letter_str='D';
syl_name='Voc84';

% % seemingly hard one, but get good estimate
% date_str='06132012';
% letter_str='D';
% syl_name='Voc110';

% % bad one -- seems like mikes 3 & 4 are just giving no information
% date_str='06052012';
% letter_str='D';
% syl_name='Voc572';

% % good one
% date_str='06062012';
% letter_str='E';
% syl_name='Voc266';

% % test voc -- get good localization on three pairs, presumably b/c
% % all signals are similar in amplitude
% date_str='06132012';
% letter_str='D';
% syl_name='Voc133';

% % problematic---mouse head outside microphone rectangle
% date_str='06052012';
% letter_str='D';
% syl_name='Voc685;

% % problematic, but well-localized
% date_str='06052012';
% letter_str='D';
% syl_name='Voc686';

% % problematic one, turns out f_lo>f_hi
% date_str='06052012';
% letter_str='D';
% syl_name='Voc256';

% % problematic one, turns out f_lo>f_hi
% date_str='06062012';
% letter_str='E';
% syl_name='Voc343';

% % problematic one, turns out f_lo>f_hi
% date_str='06062012';
% letter_str='E';
% syl_name='Voc913';

% % problematic one, turns out f_lo>f_hi
% date_str='06062012';
% letter_str='E';
% syl_name='Voc914';

% where the files are stored, etc.
base_dir_name='/groups/egnor/egnorlab/Neunuebel/ssl_sys_test';
%base_dir_name='~/egnor/ssl/ssl_sys_test';
data_analysis_dir_name='Data_analysis';
args.conf_level=0.68;
args.read_from_map_cache=true;
args.write_to_map_cache=false;
args.quantify_confidence=true;
args.return_big_things=true;
verbosity=0;

tic
r_est_blob = ...
  r_est_from_voc_indicators(base_dir_name,...
                            data_analysis_dir_name, ...
                            date_str, ...
                            letter_str, ...
                            syl_name, ...
                            args, ...
                            verbosity);
toc

% unpack the return blob
field_name=fieldnames(r_est_blob);
for i=1:length(field_name)
  eval(sprintf('%s=r_est_blob.%s;',field_name{i},field_name{i}));
end

% get dims out
n_pairs=size(rsrp_per_pair_grid,3);  %#ok
N=size(v,1);  %#ok
K=size(R,2);

% % plot the raw clips
clr_mike=[1 0 0 ; ...
          0 0.7 0 ; ...
          0 0 1 ; ...
          0 0.8 0.8 ];
dt=1/fs;  % s        
t=dt*(0:(N-1))';  % s
% figure('color','w');
% for k=1:K
%   subplot_handle=subplot(K,1,k);
%   plot(1000*t,1000*v(:,k),'color',clr_mike(k,:));
%   ylim(ylim_tight(1000*v(:,k)));
%   ylabel(sprintf('Mic %d (mV)',k));
%   if k~=K ,
%     set(subplot_handle,'xticklabel',{});
%   end
% end
% xlabel('Time (ms)');
% ylim_all_same();
% tl(0,1000*t(end));

% get the filtered signals back from the f domain
v_filt=real(ifft(V_filt));

% set up the figure and place all the axes                                   
w_fig=2.5; % in
h_fig=3; % in
n_row=K;
n_col=1;
w_axes=1.5;  % in
h_axes=0.5;  % in
w_space=1;  % in (not used)
h_space=0;  % in                              
[figure_handle,subplot_handles]= ...
  layout_axes_grid(w_fig,h_fig,...
                   n_row,n_col,...
                   w_axes,h_axes,...
                   w_space,h_space);
set(figure_handle,'color','w');                               

% plot the filtered clips with raw in background
white_fraction=0.75;
%figure_handle=figure('color','w');
%set_figure_size_explicit(figure_handle,[3 6]);
for k=1:K
  subplot_handle=subplot_handles(k);
  axes(subplot_handle);  %#ok
  plot(1000*t,1000*v(:,k)     ,'color',(1-white_fraction)*clr_mike(k,:)+white_fraction*[1 1 1]);
  hold on
  plot(1000*t,1000*v_filt(:,k),'color',clr_mike(k,:));
  hold off
  set(subplot_handle,'fontsize',7);
  ylim(ylim_tight(1000*v(:,k)));
  %ylabel(sprintf('Mic %d (mV)',k),'fontsize',7);
  if k~=K ,
    set(subplot_handle,'xticklabel',{});
    set(subplot_handle,'yticklabel',{});
  else
    set(subplot_handle,'yAxisLocation','right');
  end
end
%xlabel('Time (ms)','fontsize',7);
ylim_all_same();
tl(0,1000*t(end));

% print the per-mic RMS amplitudes
r_in_mV=1000*a  %#ok

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
% r_center=(r_head+r_tail)/2;
% a_vec=r_head-r_center;  % vector
% b=norm(a_vec)/3;  % scalar, and a guess at the half-width of the mouse

% Compute the CR borders from the objective function and critical value
%in_cr_grid=(mse_grid<=mse_crit);
%r_cr_bounds_ij=bwboundaries(in_cr_grid);
%r_cr_bounds=path_xy_from_path_ij(r_cr_bounds_ij,x_grid,y_grid);

% make a figure of the objective funtion, estimate, and CR
rsrp_abs_max=max(max(abs(rsrp_grid)));
title_str=sprintf('%s %s %s', ...
                  date_str,letter_str,syl_name);
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
axes_h=axes('parent',fig_h);
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
                            R,r_head,r_tail);
title(axes_h,title_str,'interpreter','none','fontsize',7);
h=xlabel(axes_h,'');  delete(h);
h=ylabel(axes_h,'');  delete(h);
set(axes_h,'ytick',[]);
set(axes_h,'xtick',[40 60]);
                          
axes_cb_h=add_colorbar(axes_h,w_colorbar,w_colorbar_spacer);
set(axes_cb_h,'fontsize',7);
ylabel(axes_cb_h,colorbar_label_str);



% figure_objective_map(x_grid,y_grid,10^6*rsrp_grid/N, ...
%                      @bipolar_red_white_blue, ...
%                      10^6*rsrp_abs_max/N*[-1 +1], ...
%                      title_str, ...
%                      'RSRP/sample (mV^2)', ...
%                      clr_mike, ...
%                      clr_anno, ...
%                      r_est,[], ...
%                      R,r_head,r_tail);
%r_mouse_ellipse=polygon_from_ellipse(r_center,a_vec,b);
%line(100*r_mouse_ellipse(1,:),100*r_mouse_ellipse(2,:),'color','k');
% line(100*r_est_in_body(1),100*r_est_in_body(2),0, ...
%      'marker','o','linestyle','none','color','w', ...
%      'markersize',6);

% normalize the per-pair rsrp images by the RMS amplitude of each signal
[M,iMicFirst,iMicSecond]=mixing_matrix_from_n_mics(K);
rsrp_per_pair_grid_normed=zeros(size(rsrp_per_pair_grid));
for i_pair=1:n_pairs
  keep=logical(abs(M(i_pair,:)));
  norming_factor=sqrt(prod(a(keep).^2))*N;
  rsrp_per_pair_grid_normed(:,:,i_pair)=rsrp_per_pair_grid(:,:,i_pair)/norming_factor;
end

% % plot the rsrp image for each pair
% for i_pair=1:n_pairs
%   title_str=sprintf('%s %s %s, Mic pair %d,%d', ...
%                     date_str,letter_str,syl_name,iMicFirst(i_pair),iMicSecond(i_pair));
%   rsrp_abs_max_this=max(max(abs(rsrp_per_pair_grid_normed(:,:,i_pair))));
%   figure_objective_map(x_grid,y_grid,rsrp_per_pair_grid_normed(:,:,i_pair), ...
%                        @bipolar_red_white_blue, ...
%                        rsrp_abs_max_this*[-1 +1], ...
%                        title_str, ...
%                        'Normed RSRP', ...
%                        clr_mike, ...
%                        clr_anno, ...
%                        r_est,[], ...
%                        R,r_head,r_tail);
%   r_mouse_ellipse=polygon_from_ellipse(r_center,a_vec,b);
%   line(100*r_mouse_ellipse(1,:),100*r_mouse_ellipse(2,:),'color','k');
% end

% plot the rsrp image for each pair in a single figure
max_abs=max(max(max(abs(rsrp_per_pair_grid_normed))));
% title_str_prefix=sprintf('%s %s %s, ', ...
%                          date_str,letter_str,syl_name);
[fig_h,axes_hs]= ...
  figure_objective_map_per_pair_grid(x_grid,y_grid,rsrp_per_pair_grid_normed, ...
                                     @bipolar_red_white_blue, ...
                                     max_abs*[-1 +1], ...
                                     clr_mike, ...
                                     clr_anno, ...
                                     r_est,[], ...
                                     R,r_head,r_tail);

% % make a histogram of the rsrp
% rsrp_serial=rsrp_grid(:);
% rsrp_edges=1e-6*(floor(1e6*min(rsrp_serial)):0.1:ceil(1e6*max(rsrp_serial)))';
% rsrp_centers=(rsrp_edges(1:end-1)+rsrp_edges(2:end))/2;
% n_pels_per_rsrp_bin=histc(rsrp_serial,rsrp_edges);
% n_pels_per_rsrp_bin(end)=[];
% frac_pels_per_rsrp_bin=n_pels_per_rsrp_bin/length(rsrp_serial);
% 
% figure('color','w');
% h=bar(1e6*rsrp_centers,100*frac_pels_per_rsrp_bin);
% set(h,'edgecolor','none');
% xlabel('RSRP (mV^2)');
% ylabel('Pels per bin (%)');


% % make a separate figure for the CR
% clr_anno=[1 1 1];
% figure_objective_map(x_grid,y_grid,in_cr_grid, ...
%                      @jet, ...
%                      [0 1], ...
%                      title_str, ...
%                      'In CR?', ...
%                      clr_mike, ...
%                      clr_anno, ...
%                      r_est,[], ...
%                      R,r_head,r_tail);
% r_mouse_ellipse=polygon_from_ellipse(r_center,a_vec,b);
% line(100*r_mouse_ellipse(1,:),100*r_mouse_ellipse(2,:),'color','w');
% 
% % calculate the P value as a function of position
% s=load('cdf_dJ_emp_unique.mat');
% cdf_dJ_emp_unique=s.cdf_dJ_emp_unique;
% dJ_line_unique=s.dJ_line_unique;
% dJ_grid=mse_grid/mse_min-1;
% P_grid=1-interp1_fast(dJ_line_unique,cdf_dJ_emp_unique,dJ_grid);
% 
% % plot that
% figure_objective_map(x_grid,y_grid,P_grid, ...
%                      @jet, ...
%                      [0 1], ...
%                      title_str, ...
%                      'P-value', ...
%                      clr_mike, ...
%                      clr_anno, ...
%                      r_est,[], ...
%                      R,r_head,r_tail);
% r_mouse_ellipse=polygon_from_ellipse(r_center,a_vec,b);
% line(100*r_mouse_ellipse(1,:),100*r_mouse_ellipse(2,:),'color','w');
% 
% % calculate log-odds of P value
% loP_grid=log2(P_grid./(1-P_grid));
% 
% % plot that
% figure_objective_map(x_grid,y_grid,loP_grid, ...
%                      @jet, ...
%                      [], ...
%                      title_str, ...
%                      'log-odds of P-value', ...
%                      clr_mike, ...
%                      clr_anno, ...
%                      r_est,[], ...
%                      R,r_head,r_tail);
% r_mouse_ellipse=polygon_from_ellipse(r_center,a_vec,b);
% line(100*r_mouse_ellipse(1,:),100*r_mouse_ellipse(2,:),'color','w');



