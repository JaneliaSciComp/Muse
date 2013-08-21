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
verbosity=1;

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

% show the per-mic RMS amplitudes
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

% extrapolate a mouse ellipse
r_center=(r_head+r_tail)/2;
a_vec=r_head-r_center;  % vector
b=norm(a_vec)/3;  % scalar, and a guess at the half-width of the mouse

% Compute the CR borders from the objective function and critical value
%in_cr_grid=(mse_grid<=mse_crit);
%r_cr_bounds_ij=bwboundaries(in_cr_grid);
%r_cr_bounds=path_xy_from_path_ij(r_cr_bounds_ij,x_grid,y_grid);

% make a figure of the objective funtion, estimate, and CR
rsrp_abs_max=max(max(abs(rsrp_grid)));  %#ok
title_str=sprintf('%s %s %s', ...
                  date_str,letter_str,syl_name);
clr_anno=[0 0 0];
clr_mike=[1 0 0 ; ...
          0 0.7 0 ; ...
          0 0 1 ; ...
          0 0.8 0.8 ];
figure_objective_map(x_grid,y_grid,10^6*rsrp_grid, ...
                     @bipolar_red_white_blue, ...
                     10^6*rsrp_abs_max*[-1 +1], ...
                     title_str, ...
                     'RSRP (mV^2)', ...
                     clr_mike, ...
                     clr_anno, ...
                     r_est,[], ...
                     R,r_head,r_tail);
r_mouse_ellipse=polygon_from_ellipse(r_center,a_vec,b);
line(100*r_mouse_ellipse(1,:),100*r_mouse_ellipse(2,:),'color','k');
% line(100*r_est_in_body(1),100*r_est_in_body(2),0, ...
%      'marker','o','linestyle','none','color','w', ...
%      'markersize',6);

% plot the rsrp image for each pair
n_pairs=size(rsrp_per_pair_grid,3);  %#ok
K=size(R,2);
[M,iMicFirst,iMicSecond]=mixing_matrix_from_n_mics(K);
for i_pair=1:n_pairs
  title_str=sprintf('%s %s %s, Mic pair %d,%d', ...
                    date_str,letter_str,syl_name,iMicFirst(i_pair),iMicSecond(i_pair));
  rsrp_abs_max_this=max(max(abs(rsrp_per_pair_grid(:,:,i_pair))));
  figure_objective_map(x_grid,y_grid,10^6*rsrp_per_pair_grid(:,:,i_pair), ...
                       @bipolar_red_white_blue, ...
                       10^6*rsrp_abs_max_this*[-1 +1], ...
                       title_str, ...
                       'RSRP (mV^2)', ...
                       clr_mike, ...
                       clr_anno, ...
                       r_est,[], ...
                       R,r_head,r_tail);
  r_mouse_ellipse=polygon_from_ellipse(r_center,a_vec,b);
  line(100*r_mouse_ellipse(1,:),100*r_mouse_ellipse(2,:),'color','k');
end  

% make a histogram of the rsrp
rsrp_serial=rsrp_grid(:);
rsrp_edges=1e-6*(floor(1e6*min(rsrp_serial)):0.1:ceil(1e6*max(rsrp_serial)))';
rsrp_centers=(rsrp_edges(1:end-1)+rsrp_edges(2:end))/2;
n_pels_per_rsrp_bin=histc(rsrp_serial,rsrp_edges);
n_pels_per_rsrp_bin(end)=[];
frac_pels_per_rsrp_bin=n_pels_per_rsrp_bin/length(rsrp_serial);

figure('color','w');
h=bar(1e6*rsrp_centers,100*frac_pels_per_rsrp_bin);
set(h,'edgecolor','none');
xlabel('RSRP (mV^2)');
ylabel('Pels per bin (%)');


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



