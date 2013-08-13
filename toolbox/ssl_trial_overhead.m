function [syl_name,i_start,i_end,f_lo,f_hi,r_head,r_tail,R,Temp, ...
          dx,x_grid,y_grid,in_cage]= ...
  ssl_trial_overhead(base_dir_name,data_analysis_dir_name,date_str,letter_str)

% construct the experiment dir name
% try this variant first
exp_dir_name=fullfile(base_dir_name,...
                      sprintf('sys_test_%s',date_str));
if ~exist(exp_dir_name,'dir')
  % if that didn't work, try this
  exp_dir_name=fullfile(base_dir_name,date_str);
end
                      
% read the vocalization index file                      
voc_index_file_name= ...
  fullfile(exp_dir_name, ...
           sprintf('%s/Test_%s_1_Mouse.mat', ...
                   data_analysis_dir_name, ...
                   letter_str));
[syl_name,i_start,i_end,f_lo,f_hi,r_head_pels,r_tail_pels]= ...
  read_voc_index(voc_index_file_name);
%n_voc=length(i_syl);

% get the meters per pel
meters_per_pel= ...
  load_anonymous(sprintf('%s/meters_2_pixels.mat',exp_dir_name));  % m/pel

% convert the head and tail positions
r_head=meters_per_pel*r_head_pels;  % m
r_tail=meters_per_pel*r_tail_pels;  % m

% read the cage bounds
corner_file_name= ...
  fullfile(exp_dir_name, ...
           sprintf('Test_%s_1_mark_corners.mat',letter_str));
r_corners=load_corner_file(corner_file_name);  % m

% load the microphone positions
positions_out=load_anonymous(sprintf('%s/positions_out.mat',exp_dir_name));
n_mike=4;
R=zeros(3,n_mike);  % Mike positions in cols (in meters)
R(1,:)=[positions_out.x_m];  % m
R(2,:)=[positions_out.y_m];  % m
R(3,:)=[positions_out.z_m];  % m
clear positions_out;

% velocity of sound in air
Temp = load_anonymous(sprintf('%s/temps.mat',exp_dir_name));  % C
Temp=mean(Temp);  % should fix at some point

% set the grid resolution
%dx=0.001*1;  % m
dx=0.001*0.25;  % m, the resolution we really want

% figure grid bounds
x_min=dx*floor(min(R(1,:))/dx);
x_max=dx*ceil(max(R(1,:))/dx);
y_min=dx*floor(min(R(2,:))/dx);
y_max=dx*ceil(max(R(2,:))/dx);

% make some grids and stuff
xl=[x_min x_max];  % m
yl=[y_min y_max];  % m
x_line=(xl(1):dx:xl(2))';
y_line=(yl(1):dx:yl(2))';
n_x=length(x_line);
n_y=length(y_line);
x_grid=repmat(x_line ,[1 n_y]);
y_grid=repmat(y_line',[n_x 1]);

% move the corners out by a certain amount, since the mice are sometimes
% a little bit out of the quadrilateral defined by r_corners
dx_corner=0.00;  % m
r_center=mean(r_corners,2);  % m, 2x1
v_out=bsxfun(@minus,r_corners,r_center);  % m, 2x4
v_out_sign=sign(v_out);
r_corners_nudged=bsxfun(@plus,r_corners,dx_corner*v_out_sign);

% make a mask that indicates when a grid point is within the cage
in_cage=inside_convex_poly(x_grid,y_grid,r_corners_nudged);

end
