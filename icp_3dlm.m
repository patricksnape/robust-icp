function [x, R, t] = icp_3dlm(Model, Data, ROBUST_THRESH, initial_p)

% RUN_ICP3D     A function
%               ...

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 31 Aug 01

if nargin < 3
    ROBUST_THRESH = 1000; % i.e. non-robust
end

icp.ROBUST_THRESH = ROBUST_THRESH;
N = 100;

[Model, Data] = scale_and_centre(Model, Data, N);

% Make model DT
D = pointclouddt(Model, N);

% Compute derivatives of D
dDdX = convn(D, [1,0,-1]/2, 'same');
dDdY = convn(D, [1,0,-1]'/2, 'same');
dDdZ = convn(D, cat(3, 1, 0, -1)/2, 'same');

m_estimator = 'huber';
D = sqrt(awf_m_estimator(m_estimator, D, N));
Efactor = 1 ./ (2 * D);
dDdX = Efactor .* awf_m_estimator(['D' m_estimator], dDdX, N);
dDdY = Efactor .* awf_m_estimator(['D' m_estimator], dDdY, N);
dDdZ = Efactor .* awf_m_estimator(['D' m_estimator], dDdZ, N);

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

scatter(BOX * (N+1), 'b-')
hold on
set(scatter(Model(:,[2 1 3]), 'r.'), 'markersize', 0.5)
camlight
hold on
h = scatter(Data(:,[2 1 3]), '.');
if size(Data,1) > 1000
  set(h, 'markersize', 0.1);
end
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
icp.ModelDx = dDdX;
icp.ModelDy = dDdY;
icp.ModelDz = dDdZ;
icp.handle = h;

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

%% Code:
%
PARAMSCALE = [1 1 1 1 100 100 100];

% SCALE up
params = params .* PARAMSCALE;

% Compute transformed data points and Jacobians
[Tdata, Jx, Jy, Jz] = icp_3d_err_transformed(params, icp.Data);

set(icp.handle, ...
  'xdata', Tdata(:,2), ...
  'ydata', Tdata(:,1), ...
  'zdata', Tdata(:,3));
drawnow

% Index each row of Tdata into the DT.  It's oddly ordered...
i = (Tdata(:,2));
j = (Tdata(:,1));
k = (Tdata(:,3));

dists = interp3(icp.ModelDistanceTransform, i, j, k, 'linear');
ddists_dx = interp3(icp.ModelDy, i, j, k, 'linear'); % x,y swapped
ddists_dy = interp3(icp.ModelDx, i, j, k, 'linear');
ddists_dz = interp3(icp.ModelDz, i, j, k, 'linear');

ddists_dx(~isfinite(ddists_dx)) = 0;
ddists_dy(~isfinite(ddists_dy)) = 0;
ddists_dz(~isfinite(ddists_dz)) = 0;
[a, b, c] = size(icp.ModelDistanceTransform);
dists(~isfinite(dists)) = a+b+c;

Grad_pdists = [ddists_dx ddists_dy ddists_dz];

dists = dists / 100;
Grad_pdists = Grad_pdists / 100;

% d/dp f(T(p)) = [d/dx f(x)](T(p)) d/dp T(p)
%    1 x 7             1 x 3       3 x 7

L = ones(1,7);
% Replicate each column 7 times
J = Grad_pdists(:, L*1) .* Jx + ...
    Grad_pdists(:, L*2) .* Jy + ...
    Grad_pdists(:, L*3) .* Jz;

% Scale down J
J = J .* PARAMSCALE(ones(size(J,1), 1),:);

J = double(J);
dists = double(dists);

% Print iteration information and draw scatter

global run_icp3d_iter
fprintf('Iter %3d ', run_icp3d_iter);
run_icp3d_iter = run_icp3d_iter + 1;

fprintf('%5.2f ', params);
fprintf('err %g\n', norm(dists));

set(icp.handle, ...
  'xdata', Tdata(:,2), ...
  'ydata', Tdata(:,1), ...
  'zdata', Tdata(:,3));
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

R = coolquat2mat([p1 p2 p3 p4]) / sum([p1 p2 p3 p4].^2);
t = [p5 p6 p7]';
