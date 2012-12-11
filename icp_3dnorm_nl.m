function [phi, theta] = icp_3dnorm_nl(model, data)

% Center model
model_centre = mean(model);
model = awf_translate_pts(model, -model_centre);

% Center data
data_centre = mean(data);
data = awf_translate_pts(data, -data_centre);

model_normals = lsqnormest(model', 4)';
data_normals = lsqnormest(data', 4)';

options = optimset('lsqnonlin');
options.TypicalX = [1 1];
options.TolFun = 0.0001;
options.TolX = 0.00001;
options.DiffMinChange = .001;
options.LargeScale = 'on';
options.maxFunEvals = 1000;
% options.Jacobian = 'on';
% options.DerivativeCheck = 'off';

[data_phi, data_theta, ~] = cart2sph(data_normals(:, 1), data_normals(:, 2), data_normals(:, 3));
[model_phi, model_theta, ~] = cart2sph(model_normals(:, 1), model_normals(:, 2), model_normals(:, 3));

params.model_normals = model_normals;
params.model_phi = model_phi;
params.model_theta = model_theta;
params.data_phi = data_phi;
params.data_theta = data_theta;

x = lsqnonlin(@(X) icp_norm(X, params), [0 0], [], [], options);

phi = x(1);
theta = x(2);

function [dists] = icp_norm(estimate, params)

model_normals = params.model_normals;
model_phi = params.model_phi;
model_theta = params.model_theta;
data_phi = params.data_phi + estimate(1);
data_theta = params.data_theta + estimate(2);

[x,y,z] = sph2cart(data_phi, data_theta, ones(length(data_phi), 1));

[matches] = knnsearch(model_normals, [x,y,z]);

mp = model_phi(matches, :);
dp = data_phi;

mt = model_theta(matches, :);
dt = data_theta;

dists = abs(cos(dp) - cos(mp)) + abs(sin(dp) - sin(mp)) + ...
        abs(cos(dt) - cos(mt)) + abs(sin(dt) - sin(mt));