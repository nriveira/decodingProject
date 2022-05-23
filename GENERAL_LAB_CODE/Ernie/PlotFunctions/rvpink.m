function p = rvpink(m)
%PINK   Pastel shades of pink color map
%   PINK(M) returns an M-by-3 matrix containing a "pink" colormap.
%   PINK, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(pink)
%
%   See also HSV, GRAY, HOT, COOL, BONE, COPPER, FLAG, 
%   COLORMAP, RGBPLOT.

%   C. Moler, 5-11-91, 8-19-92.
%   Copyright 1984-2004 The MathWorks, Inc.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
p = sqrt((2*gray(m) + hot(m))/3);
p = flip(p,1);
