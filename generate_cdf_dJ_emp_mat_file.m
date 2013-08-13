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

% marshall the values across trials into double arrays (not cell arrays)
[date_str_flat, ...
 letter_str_flat, ...
 i_syl_flat, ...
 r_est_flat, ...
 a_flat, ...
 mse_min_flat, ...
 mse_body_flat, ...
 ms_total_flat, ...
 r_head_flat, ...
 r_tail_flat, ...
 N_flat, ...
 N_filt_flat] = ...
  flatten_trials(r_est_blob_per_voc_per_trial, ...
                 per_trial_ancillary);

% calculate dJ_body_flat
n_voc_total=size(r_est_flat,2);
dJ_body_flat=mse_body_flat./mse_min_flat-1;

% sort the dJ_flat values
[dJ_body_flat,js]=sort(dJ_body_flat);
date_str_flat=date_str_flat(js);
letter_str_flat=letter_str_flat(js);
i_syl_flat=i_syl_flat(js);

% make a histogram of dJ_body values, to see what the real distribution is
dJ_min=0;
dJ_max=8;
ddJ_want=0.02;
[p_dJ_body,dJ_line_coarse,ddJ]= ...
  hist_from_data(dJ_body_flat,dJ_min,dJ_max,ddJ_want,'edges');

% % plot the histogram
% figure('color','w');
% bar(dJ_line,p_dJ_body);
% xlabel('\DeltaSSEN of body (pure)');
% ylabel('Probability');

% attempt to fit a scaled chi-squared distro to the points
x_bar=mean(dJ_body_flat)
s2=var(dJ_body_flat)
a=x_bar
dof=2*x_bar^2/s2

% calc the pdf, cdf of the fit
dJ_line=(dJ_min:ddJ_want/50:dJ_max)';
pdf_dJ=chi2scaledpdf(dJ_line,dof,a);
cdf_dJ_emp=empirical_cdf(dJ_body_flat,dJ_line);
cdf_dJ_fit=chi2scaledcdf(dJ_line,dof,a);

% plot the histogram with the fit
figure('color','w');
bar(dJ_line_coarse,p_dJ_body);
xlabel('dJ of body (pure)');
ylabel('Probability');
line(dJ_line,pdf_dJ*ddJ,'color','r');
xlim([0 2]);
ylim([0 0.2]);
title(sprintf('\\DeltaJ histogram for %d vocalizations',n_voc_total));

% compare the empirical CDF to the fit CDF
figure('color','w');
plot(dJ_line,cdf_dJ_emp,'b',...
     dJ_line,cdf_dJ_fit,'r')
legend('empirical','scaled chi^2 fit','location','southeast');
xlabel('dJ of body (pure)');
ylabel('Probability');
xlim([0 2]);
ylim([0 1]);
title(sprintf('\\DeltaJ CDF for %d vocalizations',n_voc_total));

% make the cdf points unique, so that its inverse function is interpolable
[cdf_dJ_emp_unique,js]=unique(cdf_dJ_fit);  % NB: Using fit, not empirical CDF!
dJ_line_unique=dJ_line(js);

% save the file
save('cdf_dJ_emp_unique.mat','cdf_dJ_emp_unique','dJ_line_unique');
