function [h_line,h_circle,h_text]=mic_symbol(axes_h,r,v,radius,str,varargin)

% Draw a microphone symbol.
%
% r the position where the circle and line touch, 2x1
% v the direction it's pointing, 2x1
% radius the radius of the circle

v=v/norm(v);  % normalize
r_center=r-radius*v;
v_perp=[v(2);-v(1)];

r_line=bsxfun(@plus,r,radius*[-v_perp v_perp]);
h_line=line('parent',axes_h, ...
            'xdata',r_line(1,:), ...
            'ydata',r_line(2,:), ...
            'color','k', ...
            'linewidth',1, ...
            varargin{:});

n_pts=361;
theta=linspace(0,2*pi,n_pts);
r_circle=bsxfun(@plus,r_center,radius*[cos(theta);sin(theta)]);
h_circle=line('parent',axes_h, ...
              'xdata',r_circle(1,:), ...
              'ydata',r_circle(2,:), ...
              'color','k', ...
              'linewidth',0.5, ...
              varargin{:});

h_text=text('parent',axes_h, ...
            'position',r_center', ...
            'horizontalalignment','center', ...
            'verticalalignment','middle', ...
            'fontsize',6, ...
            'string',str);
           
end
