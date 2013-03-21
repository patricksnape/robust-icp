function [x, R, t] = icp_3dlm(Model, Data, varargin)

% RUN_ICP3D     A function
%               ...

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 31 Aug 01
%% Parse input
inp = inputParser;

inp.addRequired('Model', @(x)isstruct(x) && isfield(x, 'vertices') && isreal(x.vertices) && size(x.vertices, 2) == 3);
inp.addRequired('Data', @(x)isstruct(x) && isfield(x, 'vertices') && isreal(x.vertices) && size(x.vertices, 2) == 3);

inp.addOptional('EstimateNormals', false, @(x)islogical(x));
inp.addOptional('initial_p', [0 0 0 1 0 0 0], @(x)isreal(x) && length(x) == 7);

inp.parse(Model, Data, varargin{:});
estimate_normals = inp.Results.EstimateNormals;
initial_p = inp.Results.initial_p;
clear('inp');

if (estimate_normals)
    Model.normals = lsqnormest(Model.vertices, 4);
else
    if ~isfield(Model, 'triangles') || ~isfield(Data, 'triangles')
        error('Must pass triangle list')
    end
    
    Model.normals = compute_normal(Model.vertices, Model.triangles);
end

%% Setup plot
set(scatter(Model.vertices, 'r.'), 'markersize', 0.5)
camlight
hold on
h = scatter(Data.vertices, '.');
view(3)
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Run minimization
%%
global run_icp3d_iter
run_icp3d_iter = 0;

m_dx = Model.normals(1, :);
m_dy = Model.normals(2, :);
m_dz = Model.normals(3, :);

m_x2 = m_dx .^ 2;
m_y2 = m_dy .^ 2;
m_z2 = m_dz .^ 2;

denom = sqrt(m_x2 + m_y2 + m_z2);
denom = median_adjusted(denom);

icp.M_dx = m_dx ./ denom;
icp.M_dy = m_dy ./ denom;
icp.M_dz = m_dz ./ denom;

icp.Data = Data;
icp.handle = h;
icp.kdObj = KDTreeSearcher(Model.vertices);
icp.EstimateNormals = estimate_normals;

options = optimset('lsqnonlin');
options.TypicalX = [1 1 1 1 1 1 1];
options.TolFun = 0.0001;
options.TolX = 0.00001;
options.DiffMinChange = .001;
options.Algorithm = 'levenberg-marquardt';
options.maxFunEvals = 1000;
options.Jacobian = 'on';
options.DerivativeCheck = 'off';

x = lsqnonlin(@(X) icp_error_with_derivs(X, icp), initial_p, [], [], options);

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

data = icp.Data.vertices;
m_dx = icp.M_dx;
m_dy = icp.M_dy;
m_dz = icp.M_dz;
N_p = length(params);

% Compute transformed data points and Jacobians
[Tdata, Jx, Jy, Jz] = icp_3d_err_transformed(params, data);

%% Interpolate distances and calculate m-estimate
idx = knnsearch(icp.kdObj, Tdata);
closest_Tdata = Tdata(idx, :);

if (icp.EstimateNormals)
    [d_dx, d_dy, d_dz] = estimate_normals(closest_Tdata);
else
    normals = compute_normal(closest_Tdata, icp.Data.triangles);
    d_dx = normals(1, :);
    d_dy = normals(2, :);
    d_dz = normals(3, :);
end

% Calculate df(g1,x[0]) / dg1,x[0]
dx2 = d_dx .^ 2;
dy2 = d_dy .^ 2;
dz2 = d_dz .^ 2;

denom = sqrt(dx2 + dy2 + dz2);
denom = median_adjusted(denom);
% df_g1_denom = sqrt(ddists_dx2 + ddists_dy2 + ddists_dz2) .^ 3;

% Prevent division by zero by adding the median
% Take the median ignoring the dead voxels
% df_g1_denom = median_adjusted(df_g1_denom);

% dF_g1x = (ddists_dy2 + ddists_dz2) ./ df_g1_denom;
% dF_g1y = (ddists_dx2 + ddists_dz2) ./ df_g1_denom;
% dF_g1z = (ddists_dx2 + ddists_dy2) ./ df_g1_denom;
% 
d_dx = d_dx ./ denom;
d_dy = d_dy ./ denom;
d_dz = d_dz ./ denom;

dists = (m_dx - d_dx) + (m_dy - d_dy) + (m_dz - d_dz);
% dF_g1x = repmat(dF_g1x, 1, N_p);
% dF_g1y = repmat(dF_g1y, 1, N_p);
% dF_g1z = repmat(dF_g1z, 1, N_p);
% Jx = dF_g1x .* dg1x_dp;
% Jy = dF_g1y .* dg1y_dp;
% Jz = dF_g1z .* dg1z_dp;
d_dx = repmat(d_dx', 1, N_p);
d_dy = repmat(d_dy', 1, N_p);
d_dz = repmat(d_dz', 1, N_p);
Jx = Jx .* d_dx;
Jy = Jy .* d_dy;
Jz = Jz .* d_dz;

%% Scale Jacobian by distance transform
J = Jx + Jy + Jz;

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