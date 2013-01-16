function h_out = scatter(cloud, style)
% SCATTER Plot a point cloud
%
% Input: 
%       cloud   (required) 3xN or 4xN matrix representing a cloud of points
%       style   (optional) style character. Defaults to '+'
% Output:
%       h_out   handle to a figure containing a plot3
%
% Author: Patrick Snape
% Date: 16 Jan 13

%% Parse input
inp = inputParser;

inp.addRequired('cloud', @(x)isreal(x) &&  ...
                         (size(x, 1) == 3 || size(x, 1) == 4 || ... 
                          size(x, 2) == 3 || size(x, 2) == 4));
inp.addOptional('style', '+', @(x)ischar(x));

inp.parse(cloud, style);
arg = inp.Results;
clear('inp');

%% Plot
[r, ~] = size(cloud);

if r == 3 || r == 4
    h_out = plot3(cloud(1, :), cloud(2, :), cloud(3, :), arg.style);
else
    h_out = plot3(cloud(:, 1), cloud(:, 2), cloud(:, 3), arg.style);
end

end