function [x, R, t] = icp_3dlm(Model, Data, initial_p)

% RUN_ICP3D     A function
%               ...

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 31 Aug 01

VOLUME_SIZE = 100;

[Model, Data] = scale_and_centre(Model, Data, VOLUME_SIZE);

% Make model DT
D = pointclouddt(Model, VOLUME_SIZE);

% Compute derivatives of D
[dDd_x, dDd_y, dDd_z] = gradient(D);
[dDd_xx, dDd_xy, dDd_xz] = gradient(dDd_x);
[dDd_yx, dDd_yy, dDd_yz] = gradient(dDd_y);
[dDd_zx, dDd_zy, dDd_zz] = gradient(dDd_z);
dDd_yx = dDd_xy;
dDd_yz = dDd_zy;
dDd_xz = dDd_zx;

clf
BOX = [
  0 0 0;
  0 1 0;
  1 1 0;
  1 0 0;
  0 0 0;
  0 0 1;
  0 1 1; 
  0 1 0; 
  0 1 1;
  1 1 1; 
  1 1 0; 
  1 1 1;
  1 0 1; 
  1 0 0; 
  1 0 1;
  0 0 1;
  ];

scatter(BOX * (VOLUME_SIZE+1), 'b-');
hold on
set(scatter(Model, 'r.'), 'markersize', 0.5)
camlight
hold on
h = scatter(Data, '.');
view(3)
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Run minimization
%%
global run_icp3d_iter
run_icp3d_iter = 0;

icp.Data = Data;
icp.ModelDistanceTransform = D;
icp.ModelDx = dDd_x;
icp.ModelDy = dDd_y;
icp.ModelDz = dDd_z;
icp.ModelDxx = dDd_xx;
icp.ModelDxy = dDd_xy;
icp.ModelDxz = dDd_xz;
icp.ModelDyx = dDd_yx;
icp.ModelDyy = dDd_yy;
icp.ModelDyz = dDd_yz;
icp.ModelDzx = dDd_zx;
icp.ModelDzy = dDd_zy;
icp.ModelDzz = dDd_zz;
icp.handle = h;
icp.volume_size = VOLUME_SIZE;
icp.paramscale = [1 1 1 1 VOLUME_SIZE VOLUME_SIZE VOLUME_SIZE];

options = optimset('lsqnonlin');
options.TypicalX = [1 1 1 1 1 1 1];
options.TolFun = 0.0001;
options.TolX = 0.00001;
options.DiffMinChange = .001;
options.LargeScale = 'on';
options.maxFunEvals = 1000;
options.Jacobian = 'on';
options.DerivativeCheck = 'off';
params = [0 0 0 1 0 0 0]; % quat, tx, ty, tz

if nargin > 3
    params = initial_p;
end

x = lsqnonlin(@(X) icp_error_with_derivs(X, icp), params, [], [], options);

[R, t] = icp_deparam(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dists, J] = icp_error_with_derivs(params, icp)
% In:
% icp.Data                    Nx7 data points
% icp.ModelDistanceTransform  DT of model points
%
% Out:
%   dists = Nx1
%       J = Nx7
%% Derivations:
%
%          TX := T(params, icp.Data(i))
% dists(i) = DT(TX);
%  J(i, j) = d[ DT(TX) ]/d[params(j)]
%          = ddX[ DT ](TX) * d[TX]/d[params(j)]
%
% [TX, Jx, Jy, Jz] = icp_3d_err_transformed(params, icp.Data);
%      Jx = N x 7
%  T(i,j) = Jx .* Dx(TX) + Jy .* Dy(TX) +  Jz .* Dz(TX);

%% Compte Jacobian and transform data
%
PARAMSCALE = icp.paramscale;
% SCALE up
params = params .* PARAMSCALE;

% Compute transformed data points and Jacobians
[Tdata, Jx, Jy, Jz] = icp_3d_err_transformed(params, icp.Data);

% Index each row of Tdata into the DT
i = (Tdata(:, 1));
j = (Tdata(:, 2));
k = (Tdata(:, 3));

%% Interpolate distances and calculate m-estimate

dists = interp3(icp.ModelDistanceTransform, i, j, k) / icp.volume_size;
ddists_dx = interp3(icp.ModelDx, i, j, k) / icp.volume_size;
ddists_dy = interp3(icp.ModelDy, i, j, k) / icp.volume_size;
ddists_dz = interp3(icp.ModelDz, i, j, k) / icp.volume_size;

% Calculate df(g1,x[0]) / dg1,x[0]
ddists_dx2 = ddists_dx .^ 2;
ddists_dy2 = ddists_dy .^ 2;
ddists_dz2 = ddists_dz .^ 2;

ddists_denom = sqrt(ddists_dx2 + ddists_dy2 + ddists_dz2);
ddists_denom = median_adjusted(ddists_denom);
df_g1_denom = sqrt(ddists_dx2 + ddists_dy2 + ddists_dz2) .^ 3;

% Prevent division by zero by adding the median
% Take the median ignoring the dead voxels
df_g1_denom = median_adjusted(df_g1_denom);

dF_g1x = (ddists_dy2 + ddists_dz2) ./ df_g1_denom;
dF_g1y = (ddists_dx2 + ddists_dz2) ./ df_g1_denom;
dF_g1z = (ddists_dx2 + ddists_dy2) ./ df_g1_denom;

ddists_dxx = interp3(icp.ModelDxx, i, j, k);
ddists_dxy = interp3(icp.ModelDxy, i, j, k);
ddists_dxz = interp3(icp.ModelDxz, i, j, k);
ddists_dyx = interp3(icp.ModelDyx, i, j, k);
ddists_dyy = interp3(icp.ModelDyy, i, j, k);
ddists_dyz = interp3(icp.ModelDyz, i, j, k);
ddists_dzx = interp3(icp.ModelDzx, i, j, k);
ddists_dzy = interp3(icp.ModelDzy, i, j, k);
ddists_dzz = interp3(icp.ModelDzz, i, j, k);

% [g1,xx, g1,xy, g1,xz] * dW_dp
N_p = length(params);
L = ones(1, N_p);
dg1x_dp = ddists_dxx(:, L) .* Jx + ddists_dxy(:, L) .* Jy + ddists_dxz(:, L) .* Jz;
dg1y_dp = ddists_dyx(:, L) .* Jx + ddists_dyy(:, L) .* Jy + ddists_dyz(:, L) .* Jz;
dg1z_dp = ddists_dzx(:, L) .* Jx + ddists_dzy(:, L) .* Jy + ddists_dzz(:, L) .* Jz;

ddists_dx = ddists_dx ./ ddists_denom;
ddists_dy = ddists_dy ./ ddists_denom;
ddists_dz = ddists_dz ./ ddists_denom;
dF_g1x = repmat(dF_g1x, 1, N_p);
dF_g1y = repmat(dF_g1y, 1, N_p);
dF_g1z = repmat(dF_g1z, 1, N_p);
Jx = dF_g1x .* dg1x_dp;
Jy = dF_g1y .* dg1y_dp;
Jz = dF_g1z .* dg1z_dp;

%% Scale Jacobian by distance transform
J = Jx + Jy + Jz;

% Scale down J
J = J .* PARAMSCALE(ones(size(J, 1), 1),:);

J = double(J);
dists = double(dists);

%% Print iteration information and draw scatter

global run_icp3d_iter
fprintf('Iter %3d ', run_icp3d_iter);
run_icp3d_iter = run_icp3d_iter + 1;

fprintf('%5.2f ', params);
fprintf('err %g\n', norm(dists));

set(icp.handle, ...
  'xdata', Tdata(:, 1), ...
  'ydata', Tdata(:, 2), ...
  'zdata', Tdata(:, 3));
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R,t] = icp_deparam(p)

p1 = p(1);
p2 = p(2);
p3 = p(3);
p4 = p(4);
p5 = p(5);
p6 = p(6);
p7 = p(7);

R = quat2rot([p1 p2 p3 p4]) / sum([p1 p2 p3 p4].^2);
t = [p5 p6 p7]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = median_adjusted(b)
bc    = b(:);
bc(bc == 0) = NaN;
m_b   = nanmedian(bc);

res = b + m_b;