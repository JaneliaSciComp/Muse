function clr=wycmgrb_of_x(x)

%  x is a col vector, on [0,1]

% calculates a color based on x, which goes:
% white
% yellow
% cyan
% magenta
% green
% red
% blue

% the colors
C=[1  1  1  ;
   1  1  0  ;...
   0  1  0  ;...
   0  1  1  ;...
   1  0  0  ;...
   1  0  1  ; ...
   0  0  1  ];
n_clr=size(C,1);
n_span=n_clr-1;
 
% the four cases to be handled
x_frac=mod(n_span*x,1);
x_floor=floor(n_span*x)+1;
special=(x==1);
x_frac(special)=1;
x_floor(special)=n_span;

% make weights
n_x=length(x);
w=zeros(n_x,n_clr);
for i=1:n_x
  w(i,x_floor(i)  )=1-x_frac(i);
  w(i,x_floor(i)+1)=  x_frac(i);
end

% do the matrix mult
clr=w*C;
