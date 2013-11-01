function axes_cb_h= ...
  place_gridded_function_in_axes(axes_h, ...
                                 x_grid,y_grid,f_grid, ...
                                 z, ...
                                 color_function, ...
                                 color_limits, ...
                                 do_add_colorbar, ...
                                 w_colorbar, ...
                                 w_colorbar_spacer, ...
                                 color_axis_label, ...
                                 r_corners)

if ~exist('do_add_colorbar','var') || isempty(do_add_colorbar) ,
  do_add_colorbar=false;
end
if ~exist('w_colorbar','var') 
  w_colorbar=[];
end
if ~exist('w_colorbar_spacer','var')
  w_colorbar_spacer=[];
end
if ~exist('color_axis_label','var') ,
  color_axis_label='';
end
if ~exist('r_corners','var') ,
  r_corners=[];
end

fig_h=get(axes_h,'parent');
colormap(fig_h,feval(color_function,256));

if ~isempty(color_limits)
  set(axes_h,'clim',color_limits);
end

xd=[x_grid(1,1) x_grid(end,1)];
yd=[y_grid(1,1) y_grid(1,end)];
image('parent',axes_h, ...
      'cdata',f_grid', ...
      'xdata',100*xd, ...
      'ydata',100*yd, ...
      'zdata',100*z, ...
      'cdatamapping','scaled');

% draw the colorbar
if do_add_colorbar ,
  axes_cb_h=add_colorbar(axes_h,w_colorbar,w_colorbar_spacer);
  set(axes_cb_h,'fontsize',7);
  ylabel(axes_cb_h,color_axis_label);
  set(axes_cb_h,'ytick',color_limits);
  if ~isempty(r_corners) ,
    scale_colorbar_to_corners(axes_cb_h, ...
                              axes_h, ...
                              100*r_corners);
  end
else
  axes_cb_h=[];
end
    
end
