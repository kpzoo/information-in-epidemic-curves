% Converts data to piecewise-constant form for plotting
function [x, y, yL, yU] = getStairs(x, y, yL, yU)

% Original x vector
xorig = x;
% Transformed y anc x vectors
[x, y] = stairs(xorig, y);
[~, yL] = stairs(xorig, yL);
[~, yU] = stairs(xorig, yU);