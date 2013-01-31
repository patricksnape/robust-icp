function dt = pointclouddt(cloud, N)

% N is the size of the volume NxNxN
% NB the x and y coordinate are swapped so that when the derivatives
%    are taken they correctly represent the order of (x, y, z). Matlab's
%    gradient function takes the first output, Fx, as the second dimension,
%    so we swap the x and y to remain consistent

B = round(cloud);
M = zeros(N,N,N);
M(sub2ind(size(M), B(:, 2), B(:, 1), B(:, 3))) = 1;
dt = bwdistX(M);

end