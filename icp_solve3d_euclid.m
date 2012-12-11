function [R, t] = icp_solve3d_euclid(M, P)

% ICP_SOLVE3D_EUCLID A function
%               ...

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 10 Apr 01

% find model centroid and deviations from centroid
mmark = M - repmat(mean(M), size(M,1), 1);;

% find data centroid and deviations from centroid
dmark = P - repmat(mean(P), size(P,1), 1);

N = dmark*transpose(mmark);

[U,D,V] = svd(N); % singular value decomposition

R = V*diag([1 1 det(U*V')])*transpose(U);
c = tr(D * diag([1 1 det(U*V')])) / var(P);
t = mean(M) - c*R*mean(P);