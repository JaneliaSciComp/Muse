function [syl_name,i_start,i_end,f_lo,f_hi,r_head_from_video,r_tail_from_video,R,Temp, ...
          dx,x_grid,y_grid,in_cage,r_corners]= ...
  ssl_trial_overhead_packaged(base_dir_name,data_analysis_dir_name,date_str,letter_str, ...
                              are_positions_in_old_style_coords, ...
                              frame_height_in_pels)

[syl_name,i_start,i_end,f_lo,f_hi,r_head_from_video,r_tail_from_video,R,Temp, ...
 dx,x_grid,y_grid,in_cage,r_corners]= ...
  ssl_trial_overhead(base_dir_name,data_analysis_dir_name,date_str,letter_str, ...
                     are_positions_in_old_style_coords, ...
                     frame_height_in_pels);

trial_overhead.syl_name=syl_name;
trial_overhead.i_start=i_start;
trial_overhead.i_end=i_end;
trial_overhead.f_lo=f_lo;
trial_overhead.f_hi=f_hi;
trial_overhead.r_head_from_video=r_head_from_video;
trial_overhead.syl_name=syl_name;
trial_overhead.syl_name=syl_name;
trial_overhead.syl_name=syl_name;
trial_overhead.syl_name=syl_name;
trial_overhead.syl_name=syl_name;
trial_overhead.syl_name=syl_name;
trial_overhead.syl_name=syl_name;



end
