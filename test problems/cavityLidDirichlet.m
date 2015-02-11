function [uv] = cavityLidDirichlet(x, y)
% lid velocity of one
uv = [[0, ones(1, length(x) - 2), 0]; zeros(1, length(x))];
end