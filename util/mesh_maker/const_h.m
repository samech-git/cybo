%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = const_h(x,y,x1,x2,y1,y2,h0)

% User defined size function specifying a constant size of h0 within the
% rectangle bounded by [x1,y1] & [x2,y2].

h = inf*ones(size(x,1),1);

in = (x>=x1)&(x<=x2)&(y>=y1)&(y<=y2);
h(in) = h0;