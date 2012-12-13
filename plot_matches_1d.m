function plot_matches_1d(A, B, varargin)

%% Parse Input
inp = inputParser;

inp.addRequired('A', @(x)isreal(x) && size(x, 2) == 1);
inp.addRequired('B', @(x)isreal(x) && size(x, 2) == 1);
inp.addOptional('tol', 0.2, @(x)x > 0 && x < 2 * pi);

inp.parse(A,B,varargin{:});
arg = inp.Results;
clear('inp');

%% Plot figure

hndl = figure; % creates a plotting window and stores the handle in hndl
set(hndl,'Renderer','OpenGL');

idx = find(abs(A - B) > arg.tol);
N = size(idx, 1);

plot(idx, A(idx), 'r.');
hold on;
plot(idx, B(idx), 'g.');

for i = 1:N
    line([idx(i), idx(i)], ...
         [A(idx(i)), B(idx(i))]);
end
hold off;