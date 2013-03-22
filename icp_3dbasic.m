function [params, R,t] = icp_3dbasic(model, data)

% ICP_3DLM   A function
%               ...

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 10 Apr 01

model_centre = mean(model);
model = awf_translate_pts(model, -model_centre);

data_centre = mean(data);
data = awf_translate_pts(data, -data_centre);

lmdata.kdObj = KDTreeSearcher(model);
lmdata.data = data;

figure(1)
hold off
set(scatter(model, 'b.'), 'markersize', .001);
set(gcf, 'renderer', 'opengl')
hold on
axis off
lmdata.h = scatter(data, 'r+');
set(lmdata.h, 'markersize', 2);
axis equal
axis vis3d

global run_icp3d_iter
run_icp3d_iter = 0;

if nargout < 2
  disp('No return values, returning....');
  return
end

% Set up levmarq and tallyho
options = optimset('lsqnonlin');
options.TypicalX = [1 1 1 1 1 1 1];
options.TolFun = 0.0001;
options.TolX = 0.00001;
options.DiffMinChange = .001;
options.LargeScale = 'on';
params = [0 0 0 1 0 0 0]; % quat, tx, ty, tz

params = lsqnonlin(@(X) icp_3derror(X, lmdata), params, [], [], options);
[R,t] = icp_deparam(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dists = icp_3derror(params, lm)
% 1. Extract R, t
[R,t] = icp_deparam(params);

t = t(:)'; % colvec

% 2. Evaluate

D = lm.data;
D = D * R' + t(ones(1, size(D, 1)), :);

[~, dists] = knnsearch(lm.kdObj, D);

stdDists = std(dists);
dists = awf_m_estimator('ls', dists, stdDists);

global run_icp3d_iter
fprintf('Iter %3d ', run_icp3d_iter);
run_icp3d_iter = run_icp3d_iter + 1;

fprintf('%5.2f ', params);
fprintf('err %.2f\n', norm(dists));

set(lm.h, ...
  'xdata', D(:, 1), ...
  'ydata', D(:, 2), ...
  'zdata', D(:, 3));
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

