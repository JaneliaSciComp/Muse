function trial_overhead= ...
  ssl_trial_overhead_packaged(base_dir_name,data_analysis_dir_name,date_str,letter_str, ...
                              are_positions_in_old_style_coords, ...
                              frame_height_in_pels)

% Call the function that sorts all this stuff out                            
[tf_rect_name,i_start,i_end,f_lo,f_hi,r_head_from_video,r_tail_from_video,R,Temp, ...
 dx,x_grid,y_grid,in_cage,r_corners,fs]= ...
  ssl_trial_overhead(base_dir_name,data_analysis_dir_name,date_str,letter_str, ...
                     are_positions_in_old_style_coords, ...
                     frame_height_in_pels);

% Package it all up in a scalar struct
trial_overhead.tf_rect_name=tf_rect_name;
trial_overhead.i_start=i_start;
trial_overhead.i_end=i_end;
trial_overhead.f_lo=f_lo;
trial_overhead.f_hi=f_hi;
trial_overhead.r_head_from_video=r_head_from_video;
trial_overhead.r_tail_from_video=r_tail_from_video;
trial_overhead.R=R;
trial_overhead.Temp=Temp;
trial_overhead.dx=dx;
trial_overhead.x_grid=x_grid;
trial_overhead.y_grid=y_grid;
trial_overhead.in_cage=in_cage;
trial_overhead.r_corners=r_corners;
trial_overhead.fs=fs;

end
