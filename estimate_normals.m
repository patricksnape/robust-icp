function [dx, dy, dz] = estimate_normals(data)

n = lsqnormest(data', 4);

dx = n(1, :);
dy = n(2, :);
dz = n(3, :);