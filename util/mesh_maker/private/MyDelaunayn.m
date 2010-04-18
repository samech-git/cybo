function t = MyDelaunayn(p)

% My version of the MATLAB delaunayn function that attempts to deal with
% some of the compatibility and robustness problems.
%
% Darren Engwirda - 2007

% Translate to the origin and scale the min xy range onto [-1,1]
% This is absolutely critical to avoid precision issues for large problems!
maxxy = max(p);
minxy = min(p);
p(:,1) = p(:,1)-0.5*(minxy(1)+maxxy(1));
p(:,2) = p(:,2)-0.5*(minxy(2)+maxxy(2));
p = p/(0.5*min(maxxy-minxy));

try
   % Use the default settings. This will be 'Joggled Input' prior to
   % MATLAB 7.0 or 'Triangulated Output' in newer versions
   t = delaunayn(p);
catch
   if ~(str2double(version('-release'))<=13)
      % 'Triangulated Output' is generally more accurate, but less robust.
      % If Qhull crashes with 'Triangulated Output', redo with 'Joggled
      % Input'
      t = delaunayn(p,{'QJ','Pp'});
   end
end
