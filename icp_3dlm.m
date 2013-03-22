function [x, R, t, error] = icp_3dlm(model, data, initial_p)

% RUN_ICP3D     A function
%               ...

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 31 Aug 01

VOLUME_SIZE = 100;

[Model, Data] = scale_and_centre(model, data, VOLUME_SIZE);

% Make model DT
D = pointclouddt(Model, VOLUME_SIZE);

% Compute derivatives of D
[dDdX, dDdY, dDdZ] = gradient(D);

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
set(scatter(Model, 'r.'), 'markersize', 0.5);
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
icp.ModelDx = dDdX;
icp.ModelDy = dDdY;
icp.ModelDz = dDdZ;
icp.handle = h;
icp.volume_size = VOLUME_SIZE;
icp.paramscale = [1 1 1 1 VOLUME_SIZE VOLUME_SIZE VOLUME_SIZE];

options = optimset('lsqnonlin');
options.TypicalX = [1 1 1 1 1 1 1];
options.TolFun = 0.0001;
options.TolX = 0.00001;
options.DiffMinChange = .001;
options.Algorithm = 'levenberg-marquardt';
options.Jacobian = 'on';
options.DerivativeCheck = 'off';
params = [0 0 0 1 0 0 0]; % quat, tx, ty, tz

if nargin > 3
    params = initial_p;
end

x = lsqnonlin(@(X) icp_error_with_derivs(X, icp), params, [], [], options);

[R, t] = icp_deparam(x);

fit = (R * data' + repmat(t, [1 size(data, 1)]))';
[~, error] = knnsearch(model, fit);
error = norm(error);

clf;
scatter(model, 'r.');
hold on;
scatter(fit, 'b.');

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

%% Compute distance penalties near boundaries
distpenalty = zeros(size(Tdata, 1), 1);
Grad_distpenalty = zeros(size(Tdata, 1), 3);

LO = 3;
HI = icp.volume_size - 3;

% i < 1 .. clip and add pythagorean penalty
I = find(i < LO);
if ~isempty(I)
  distpenalty(I) = distpenalty(I) + (LO - i(I)).^2;
  Grad_distpenalty(I, 1) = Grad_distpenalty(I, 1) + 2 * (LO - i(I));
  i(I) = LO;
end

% j < LO
I = find(j < LO);
if ~isempty(I)
  distpenalty(I) = distpenalty(I) + (LO - j(I)).^2;
  Grad_distpenalty(I, 2) = Grad_distpenalty(I, 2) + 2 * (LO - j(I));
  j(I) = LO;
end

% k < LO
I = find(k < LO);
if ~isempty(I)
  distpenalty(I) = distpenalty(I) + (LO - k(I)).^2;
  Grad_distpenalty(I, 3) = Grad_distpenalty(I, 3) + 2 * (LO - k(I));
  k(I) = LO;
end

% i >= HI .. clip and add pythagorean penalty
I = find(i >= HI);
if ~isempty(I)
  distpenalty(I) = distpenalty(I) + (i(I) - HI).^2;
  Grad_distpenalty(I, 1) = Grad_distpenalty(I, 1) + 2 * (i(I) - HI);
  i(I) = HI;
end

% j >= HI
I = find(j >= HI);
if ~isempty(I)
  distpenalty(I) = distpenalty(I) + (j(I) - HI).^2;
  Grad_distpenalty(I, 2) = Grad_distpenalty(I, 2) + 2 * (j(I) - HI);
  j(I) = HI;
end

% k >= HI
I = find(k >= HI);
if ~isempty(I)
  distpenalty(I) = distpenalty(I) + (k(I) - HI).^2;
  Grad_distpenalty(I, 3) = Grad_distpenalty(I, 3) + 2 * (k(I) - HI);
  k(I) = HI;
end

%% Interpolate distances and calculate m-estimate

dists = interp3(icp.ModelDistanceTransform, i, j, k, 'linear');
ddists_dx = interp3(icp.ModelDx, i, j, k, 'linear');
ddists_dy = interp3(icp.ModelDy, i, j, k, 'linear');
ddists_dz = interp3(icp.ModelDz, i, j, k, 'linear');

% Add penalties to distances, so points outside the DT get
% an upper-bound distance to the model.
dists = dists + 1e-4;
pdists = sqrt(awf_m_estimator('ls', dists, icp.volume_size) + distpenalty);
mprime = awf_m_estimator('Dls', dists, icp.volume_size);
c = mprime ./ (2 * pdists);

dpdists_dx = c .* (ddists_dx + Grad_distpenalty(:, 1));
dpdists_dy = c .* (ddists_dy + Grad_distpenalty(:, 2));
dpdists_dz = c .* (ddists_dz + Grad_distpenalty(:, 3));
Grad_pdists = [dpdists_dx dpdists_dy dpdists_dz];

% Scale back distances
dists = pdists / icp.volume_size;
Grad_pdists = Grad_pdists / icp.volume_size;

%% Scale Jacobian by distance transform

% d/dp f(T(p)) = [d/dx f(x)](T(p)) d/dp T(p)
%    1 x 7             1 x 3       3 x 7

L = ones(1,7);
% Replicate each column 7 times
J = Grad_pdists(:, L*1) .* Jx + ...
    Grad_pdists(:, L*2) .* Jy + ...
    Grad_pdists(:, L*3) .* Jz;

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
