function [phi, theta] = icp_3dnorm(model, data)

% phi = azimuth
% theta = elevation

% Center model
model_centre = mean(model);
model = awf_translate_pts(model, -model_centre);

% Center data
data_centre = mean(data);
data = awf_translate_pts(data, -data_centre);

model_normals = lsqnormest(model', 4)';
data_normals = lsqnormest(data', 4)';

[data_phi, data_theta, ~] = cart2sph(data_normals(:, 1), data_normals(:, 2), data_normals(:, 3));
[model_phi, model_theta, ~] = cart2sph(model_normals(:, 1), model_normals(:, 2), model_normals(:, 3));

%% Build the kd-tree
m = [cos(model_phi),     sin(model_phi),  ... 
     cos(model_theta), sin(model_theta)];
kdOBJ = KDTreeSearcher(m);

%% Setup variables

eps = 0.00001;

phi = 0;
theta = 0;

delta_phi = inf;
delta_theta = inf;

iterations = 0;

%% Perform search

while abs(delta_phi) > eps && abs(delta_theta) > eps && iterations < 30
    iterations = iterations + 1;
    
    new_phi = data_phi + phi;
    new_theta = data_theta + theta;
    
    d = [cos(new_phi),   sin(new_phi),  ... 
         cos(new_theta), sin(new_theta)];
    [matches] = knnsearch(kdOBJ, d);
    
    model_phi_matches = model_phi(matches);
    model_theta_matches = model_theta(matches);
    
    N = size(model_phi_matches, 1);
    delta_phi = sum(sin(data_phi + repmat(phi, [N 1])) - model_phi_matches) / N;
    delta_theta = sum(sin(data_theta + repmat(theta, [N 1])) - model_theta_matches) / N;
    
    phi = phi + delta_phi;
    theta = theta + delta_theta;
    
    err = norm(cos(data_phi + phi) - cos(model_phi))^2 + norm(sin(data_phi + phi) - sin(model_phi))^2 ...
          + norm(cos(data_theta + theta) - cos(model_theta))^2 + norm(sin(data_theta + theta) - sin(model_theta))^2;
      
    fprintf('Iteration %d - Phi: %d Theta: %d err: %d \n', iterations, phi, theta, err); 
end