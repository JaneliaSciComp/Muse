% load pre-calculated stuff
load('r_est_for_single_mouse_data_snippetized.mat');
% r_est_blob_per_segment_per_trial and overhead_per_trial are now in the
% scope

% Unpack the blob, and flatten across trials
[date_str, ...
 letter_str, ...
 i_trial, ...
 i_segment_within_trial, ...
 localized, ...
 r_est, ...
 Covariance_matrix, ...
 p_head, ...
 r_head, ...
 r_tail] = ...
  flatten_trials_given_r_est_blob_per_segment_per_trial(r_est_blob_per_segment_per_trial);
% r_head and r_tail are from the video

% Check that all the segments are there
n_segments_total=length(localized)  %#ok
n_segments_localized=sum(double(localized))  %#ok

% Calculate fraction localized
frac_localized=mean(double(localized))  %#ok

% filter out unlocalized segments
date_str=date_str(localized);
letter_str=letter_str(localized);
i_trial=i_trial(localized);
i_segment_within_trial=i_segment_within_trial(localized);
r_est=r_est(:,localized);
Covariance_matrix=Covariance_matrix(:,:,localized);
p_head=p_head(localized);
r_head=r_head(:,localized);
r_tail=r_tail(:,localized);
clear localized

% calculate errors
r_chest=3/4*r_head+1/4*r_tail;
  % r_chest is what Josh uses in the multi-mouse code
dr_est=r_est-r_chest;

% make a histogram of the error magnitudes
e2=sum(dr_est.^2,1);  % m
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
ylim([0 300]);
title(sprintf('Error histogram for %d localized segments',n_segments_localized));
drawnow;

% calc RMS error, other measures
e_mag_mean=mean(e_mag)  %#ok
e_mag_rms=rms(e_mag)  %#ok
e_mag_median=median(e_mag)  %#ok
text(20,275   ,sprintf('Median error: %0.1f cm',100*e_mag_median));
text(20,275-15,sprintf('RMS error: %0.1f cm',100*e_mag_rms));

% translate so that r_chest is in the center
dr_head=r_head-r_chest;
dr_tail=r_tail-r_chest;

% rotate around so mouse butt is to the bottom
theta=atan2(dr_tail(2,:),dr_tail(1,:));
dr_est_rot=zeros(size(dr_est));
dr_tail_rot=zeros(size(dr_tail));
dr_head_rot=zeros(size(dr_head));
for k=1:n_segments_localized 
  A=[cos(-theta(k)-pi/2) -sin(-theta(k)-pi/2) ; ...
     sin(-theta(k)-pi/2)  cos(-theta(k)-pi/2)];  % rotation matrix for -theta rotation, then pi
  dr_est_rot(:,k)=A*dr_est(:,k);
  dr_tail_rot(:,k)=A*dr_tail(:,k);
  dr_head_rot(:,k)=A*dr_head(:,k);
end

% plot those
figure('color','w');
plot(100*dr_est_rot(1,:),100*dr_est_rot(2,:),'b.');
hold on;
plot(100*dr_tail_rot(1,:),100*dr_tail_rot(2,:),'.','color',[0 0.7 0]);
plot(100*dr_head_rot(1,:),100*dr_head_rot(2,:),'+','color','r');
hold off;
%line(0,0,'marker','+','color','r','linestyle','none');
xlim([-50 +50]);
ylim([-50 +50]);
axis square;
ylabel('Distance in front of mouse (cm)');
xlabel('Distance to right of mouse (cm)');
set(gca,'xtick',-50:10:+50);
set(gca,'ytick',-50:10:+50);
title(sprintf('Vectorial error for %d localized segments, aligned to mouse', ...
              n_segments_localized));
legend({'Estimate','Mouse tail','Mouse head'},'location','southeast');
drawnow;



%
% Generate three random mouse locations for each localized segment, figure 
% out what fraction are "assignable", and of those, what fraction are
% correct
%

% make up some fake mouse locations
n_trials=max(i_trial);
n_fake_mice=3;
r_head_fake=zeros(2,n_fake_mice,0);
r_tail_fake=zeros(2,n_fake_mice,0);
for i_trial_this=1:n_trials
  is_segment_in_this_trial=(i_trial==i_trial_this);
  n_segments_this=sum(is_segment_in_this_trial);
  r_ht_this=r_head(:,is_segment_in_this_trial)-r_tail(:,is_segment_in_this_trial);
  mouse_length_this=hypot(r_ht_this(1,:,:),r_ht_this(2,:,:));
  mean_mouse_length_this=mean(mean(mouse_length_this));  
  r_corners_this=overhead_per_trial(i_trial).r_corners;
  [r_head_fake_this,r_tail_fake_this]= ...
    random_mouse_locations(r_corners_this,mean_mouse_length_this,n_fake_mice,n_segments_this);
  r_head_fake=cat(3,r_head_fake,r_head_fake_this);
  r_tail_fake=cat(3,r_tail_fake,r_tail_fake_this);  
end

% combine real and fake mouse locations
r_head_with_fake=zeros(2,1+n_fake_mice,n_segments_localized);
r_tail_with_fake=zeros(2,1+n_fake_mice,n_segments_localized);
r_head_with_fake(:,1,:)=r_head;
r_head_with_fake(:,2:end,:)=r_head_fake;
r_tail_with_fake(:,1,:)=r_tail;
r_tail_with_fake(:,2:end,:)=r_tail_fake;

% calculate chest location for real and fake
r_chest_with_fake=3/4*r_head_with_fake+1/4*r_tail_with_fake;

% calculate the density at the real+fake mice, and the posterior
% probability
r_est_big=repmat(reshape(r_est,[2 1 n_segments_localized]),[1 1+n_fake_mice 1]);
Covariance_matrix_big= ...
  repmat(reshape(Covariance_matrix,[2 2 1 n_segments_localized]),[1 1 1+n_fake_mice 1]);

r_chest_with_fake_serial= ...
  reshape(r_chest_with_fake,[2 (1+n_fake_mice)*n_segments_localized]);
r_est_big_serial=reshape(r_est_big,[2 (1+n_fake_mice)*n_segments_localized]);
Covariance_matrix_big_serial=reshape(Covariance_matrix_big,[2 2 (1+n_fake_mice)*n_segments_localized]);

p_chest_serial=mvnpdf(r_chest_with_fake_serial',r_est_big_serial',Covariance_matrix_big_serial);  % density, units: 1/(m^2)
p_chest=reshape(p_chest_serial,[(1+n_fake_mice) n_segments_localized]);
P_posterior_chest= ...
  bsxfun(@rdivide,p_chest,sum(p_chest,1));  % posterior probability

[P_posterior_assigned_maybe,i_mouse_assigned_maybe]=max(P_posterior_chest,[],1);
p_chest_assigned_maybe=max(p_chest,[],1);
p_chest_assigned_thresh=0.7;  % 1/(m^2)
is_assigned= (P_posterior_assigned_maybe>0.95) & (p_chest_assigned_maybe>p_chest_assigned_thresh);
n_assigned=sum(is_assigned)
frac_assigned=n_assigned/n_segments_localized
is_assigned_correctly=is_assigned & (i_mouse_assigned_maybe==1);
n_assigned_correctly=sum(is_assigned_correctly)
frac_assigned_correctly=n_assigned_correctly/n_assigned

