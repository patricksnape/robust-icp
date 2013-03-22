function [x, R, t, error] = icp_3dlm_normals(Model, Data, varargin)

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
    Model.normals = lsqnormest(Model.vertices', 4);
else
    if ~isfield(Model, 'triangles') || ~isfield(Data, 'triangles')
        error('Must pass triangle list')
    end
    
    Model.normals = compute_normal(Model.vertices, Model.triangles);
end

%% Setup plot
clf
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

m_dx = Model.normals(1, :)';
m_dy = Model.normals(2, :)';
m_dz = Model.normals(3, :)';

m_x2 = m_dx .^ 2;
m_y2 = m_dy .^ 2;
m_z2 = m_dz .^ 2;

denom = sqrt(m_x2 + m_y2 + m_z2);

m_dx = m_dx ./ denom;
m_dy = m_dy ./ denom;
m_dz = m_dz ./ denom;
icp.Mnormals = [m_dx, m_dy, m_dz]';

icp.Data = Data;
icp.handle = h;
icp.kdObj = KDTreeSearcher(Model.vertices);
icp.EstimateNormals = estimate_normals;

options = optimset('lsqnonlin');
options.TypicalX = [1 1 1 1 1 1 1];
options.TolFun = 0.0001;
options.TolX = 0.00001;
options.DiffMinChange = .0001;
options.Algorithm = 'levenberg-marquardt';
options.maxFunEvals = 1000;
options.Jacobian = 'on';
options.DerivativeCheck = 'off';

x = lsqnonlin(@(X) icp_error_with_derivs(X, icp), initial_p, [], [], options);

[R, t] = icp_deparam(x);

[~, error] = knnsearch(Model.vertices, (R * Data.vertices' + repmat(t, [1 size(Data.vertices, 1)]))');
error = norm(error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dists, J] = icp_error_with_derivs(params, icp)
%% Compte Jacobian and transform data
%

data = icp.Data.vertices;
N_p = length(params);

% Compute transformed data points and Jacobians
[Tdata, Jx, Jy, Jz] = icp_3d_err_transformed(params, data);

%% Interpolate distances and calculate m-estimate
if (icp.EstimateNormals)
    normals = lsqnormest(Tdata', 4);
else
    normals = compute_normal(Tdata, icp.Data.triangles);
end

[idx, pdists] = knnsearch(icp.kdObj, Tdata);
closest_norms = icp.Mnormals(:, idx);
m_dx = closest_norms(1, :)';
m_dy = closest_norms(2, :)';
m_dz = closest_norms(3, :)';

d_dx = normals(1, :)';
d_dy = normals(2, :)';
d_dz = normals(3, :)';

% Calculate df(g1,x[0]) / dg1,x[0]
d_x2 = d_dx .^ 2;
d_y2 = d_dy .^ 2;
d_z2 = d_dz .^ 2;

denom = sqrt(d_x2 + d_y2 + d_z2);
df_g1_denom = denom .^ 3;

dF_g1x = (d_y2 + d_z2) ./ df_g1_denom;
dF_g1y = (d_x2 + d_z2) ./ df_g1_denom;
dF_g1z = (d_x2 + d_y2) ./ df_g1_denom;

d_dx = d_dx ./ denom;
d_dy = d_dy ./ denom;
d_dz = d_dz ./ denom;

% Calculate residuals
dists = (m_dx - d_dx) + (m_dy - d_dy) + (m_dz - d_dz);

dx = m_dx - d_dx;
dy = m_dy - d_dy;
dz = m_dz - d_dz;

dF_g1x = repmat(dx, 1, N_p);
dF_g1y = repmat(dy, 1, N_p);
dF_g1z = repmat(dz, 1, N_p);
Jx = dF_g1x .* Jx;
Jy = dF_g1y .* Jy;
Jz = dF_g1z .* Jz;

%% Scale Jacobian by distance transform
J = Jx + Jy + Jz;

J = double(J);
dists = double(dists);

%% Print iteration information and draw scatter

global run_icp3d_iter
fprintf('Iter %3d ', run_icp3d_iter);
run_icp3d_iter = run_icp3d_iter + 1;

fprintf('%5.2f ', params);
fprintf('err %g\n', norm(pdists));

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