function [nose_x nose_y tail_x tail_y xout yout aout bout thetaout] = find_nose_tail_jpn3(m1,indx,mouse_num)
% find_nose: returns x,y position of nose for each frame in indx
%
% form: m1=find_nose(m1,indx)
%
% m1 is a struct with fields x,y,a,b, and theta
% x,y are the x/y position of the center of the mouse ellipse in meters
% a,b are the major and minor axes of the ellipse (not radii)
% theta is the clockwise angle from the positive x-axis

x = m1(1,mouse_num).m_afX;
y = m1(1,mouse_num).m_afY; %y position
theta = m1(1,mouse_num).m_afTheta; %theta;
a = m1(1,mouse_num).m_afA;
b = m1(1,mouse_num).m_afB;

h=a(indx(:))/2;
xout = x(indx(:));
yout = y(indx(:));
aout = a(indx(:));
bout = b(indx(:));
thetaout = theta(indx(:));
%find nose
xn=(h.*cos(theta(indx(:))));
yn=-(h.*sin(theta(indx(:))));
nose_x=x(indx(:))+xn;
nose_y=y(indx(:))+yn;

%find tail
xn=-(h.*cos(theta(indx(:))));
yn=(h.*sin(theta(indx(:))));
tail_x=x(indx(:))+xn;
tail_y=y(indx(:))+yn;


