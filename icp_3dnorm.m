function [phi, theta] = icp_3dnorm(model, data)

% Center model
model_centre = mean(model);
model = awf_translate_pts(model, -model_centre);

% Center data
data_centre = mean(data);
data = awf_translate_pts(data, -data_centre);

model_normals = lsqnormest(model', 4)';
data_normals = lsqnormest(data', 4)';

[data_phi, data_theta, r] = cart2sph(data_normals(:, 1), data_normals(:, 2), data_normals(:, 3));

eps = 0.00001;

phi = 0;
theta = 0;

delta_phi = inf;
delta_theta = inf;

iterations = 0;

while abs(delta_phi) > eps && abs(delta_theta) > eps
    iterations = iterations + 1;
    
    [x,y,z] = sph2cart(data_phi + phi, data_theta + theta, r);
    [matches] = knnsearch(model_normals, [x,y,z]);
    
    model_matches = model_normals(matches, :);
    
    [model_phi, model_theta, ~] = cart2sph(model_matches(:, 1), model_matches(:, 2), model_matches(:, 3));
    
    [N, ~] = size(model_matches);
    delta_phi = sum(sin(data_phi + repmat(phi, [size(model, 1) 1])) - model_phi) / N;
    delta_theta = sum(sin(data_theta + repmat(theta, [size(model, 1) 1])) - model_theta) / N;
    
    phi = phi + delta_phi;
    theta = theta + delta_theta;
    
    err = norm(cos(data_phi + phi) - cos(model_phi))^2 + norm(sin(data_phi + phi) - sin(model_phi))^2 ...
          + norm(cos(data_theta + theta) - cos(model_theta))^2 + norm(sin(data_theta + theta) - sin(model_theta))^2;
      
    fprintf('Iteration %d - Phi: %d Theta: %d err: %d \n', iterations, phi, theta, err); 
end

[x,y,z] = sph2cart(data_phi + phi, data_theta + theta, r);
[matches] = knnsearch(model_normals, [x,y,z]);
model_matches = model_normals(matches, :);

plot_matches([x,y,z], model_matches);