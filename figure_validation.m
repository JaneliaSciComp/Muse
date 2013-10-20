% load pre-calculated stuff
load('r_est_raw_for_single_mouse_data_snippetized.mat');

% Unpack the blob, and flatten across trials
[date_str, ...
 letter_str, ...
 i_segment_within_trial, ...
 localized, ...
 r_est, ...
 r_head_from_video, ...
 r_tail_from_video] = ...
  flatten_trials_given_r_est_blob_per_segment_per_trial(r_est_blob_per_segment_per_trial);

% Calculate fraction localized
frac_localized=mean(double(localized))  %#ok

% filter out unlocalized segments
date_str=date_str(localized);
letter_str=letter_str(localized);
i_segment_within_trial=i_segment_within_trial(localized);
r_est=r_est(:,localized);
r_head_from_video=r_head_from_video(:,localized);
r_tail_from_video=r_tail_from_video(:,localized);

% calculate errors
n_voc_total=size(r_est,2);
e=r_est-r_head_from_video;
v_mouse=(r_tail_from_video-r_head_from_video);  % head-to-tail vector

% plot errors
figure('color','w');
plot(100*e(1,:),100*e(2,:),'b.');
hold on;
plot3(100*v_mouse(1,:),100*v_mouse(2,:),repmat(-1,[1 n_voc_total]), ...
      '.','color',[0 0.7 0]);
hold off;
line(0,0,'marker','+','color','r','linestyle','none');
xlim([-50 +50]);
ylim([-50 +50]);
axis square;
xlabel('x (cm)');
ylabel('y (cm)');
set(gca,'xtick',-50:10:+50);
set(gca,'ytick',-50:10:+50);
title(sprintf('Vectorial error for %d vocalizations',n_voc_total));
legend({'Estimate','Mouse tail','Mouse head'},'location','southeast');
%legend({'Estimate','Mouse head'},'location','southeast');

% make a histogram of the error magnitudes
e2=sum(e.^2,1);  % m
e_mag=sqrt(e2);  % m
dist_edges=(0:0.005:1)';
dist_bin_counts=histc(e_mag,dist_edges);
dist_bin_counts=dist_bin_counts(1:end-1);
dist_centers=(dist_edges(1:end-1)+dist_edges(2:end))/2;
figure('color','w');
h=bar(100*dist_centers,dist_bin_counts);
set(h,'edgecolor','none');
set(h,'facecolor','k');
xlabel('Error (cm)');
ylabel('Frequency (counts)');
xlim([0 30]);
ylim([0 400]);
title(sprintf('Error histogram for %d vocalizations',n_voc_total));

% calc RMS error, other measures
e_mag_mean=mean(e_mag)  %#ok
e_mag_rms=rms(e_mag)  %#ok
e_mag_median=median(e_mag)  %#ok
text(20,300,sprintf('Median error: %0.1f cm',100*e_mag_median));
text(20,285,sprintf('RMS error: %0.1f cm',100*e_mag_rms));

% make a histogram of the squared errors
e2_bin_width=0.0001;
e2_edges=(0:e2_bin_width:1)';
e2_bin_counts=histc(e2,e2_edges);
e2_bin_counts=e2_bin_counts(1:end-1);
e2_bin_centers=(e2_edges(1:end-1)+e2_edges(2:end))/2;

% fit in a sloppy way to a chi-squared distro
%x_bar=mean(e2)
%s2=var(e2)
%a=x_bar  % scale factor
%dof=2*x_bar^2/s2
a=2*(0.02)^2;  % m^2
dof=2;

% calc the pdf, cdf of the fit
de2=0.0001/50;
e2_line=(0:de2:0.7)';
pdf_e2_fit=chi2scaledpdf(e2_line,dof,a);
cdf_e2_emp=empirical_cdf(e2,e2_line);
cdf_e2_fit=chi2scaledcdf(e2_line,dof,a);

% plot the histogram, with the fit pdf
figure('color','w');
h=bar(1e4*e2_bin_centers,e2_bin_counts);
set(h,'edgecolor','none');
set(h,'facecolor','b');
xlabel('Error^2 (cm^2)');
ylabel('Frequency (counts)');
xlim([0 200]);
%ylim([0 400]);
title(sprintf('Squared error histogram for %d vocalizations',n_voc_total));
line(1e4*e2_line,pdf_e2_fit*e2_bin_width*n_voc_total,'color','r');

% compare the empirical CDF to the fit CDF
figure('color','w');
plot(1e4*e2_line,cdf_e2_emp,'b',...
     1e4*e2_line,cdf_e2_fit,'r')
legend('empirical','scaled \chi^2 fit','location','southeast');
xlabel('Squared error (cm^2)');
ylabel('Probability');
xlim([0 200]);
ylim([0 1]);

% rotate around so mouse butt is to the left
theta=atan2(v_mouse(2,:),v_mouse(1,:));
e_rot=zeros(size(e));
v_mouse_rot=zeros(size(v_mouse));
for k=1:n_voc_total 
  A=[cos(-theta(k)+pi) -sin(-theta(k)+pi) ; ...
     sin(-theta(k)+pi)  cos(-theta(k)+pi)];  % rotatation matrix for -theta rotation, then pi
  e_rot(:,k)=A*e(:,k);
  v_mouse_rot(:,k)=A*v_mouse(:,k);
end

% plot those
figure('color','w');
plot(100*e_rot(1,:),100*e_rot(2,:),'b.');
hold on;
plot(100*v_mouse_rot(1,:),100*v_mouse_rot(2,:),'.','color',[0 0.7 0]);
hold off;
line(0,0,'marker','+','color','r','linestyle','none');
xlim([-50 +50]);
ylim([-50 +50]);
axis square;
xlabel('Distance in front of mouse (cm)');
ylabel('Distance to left of mouse (cm)');
set(gca,'xtick',-50:10:+50);
set(gca,'ytick',-50:10:+50);
title(sprintf('Vectorial error for %d vocalizations, aligned to mouse', ...
              n_voc_total));
legend({'Estimate','Mouse tail','Mouse head'},'location','southeast');
