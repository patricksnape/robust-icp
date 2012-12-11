function [model, data] = scale_and_centre(m, d, n)

% Center model
model_centre = mean(m);
m = awf_translate_pts(m, -model_centre);

% Center data
data_centre = mean(d);
d = awf_translate_pts(d, -data_centre);

bbox_min = min(min(m), min(d));
bbox_max = max(max(m), max(d));

maxdiff = max(bbox_max - bbox_min);

% Map to integer grid, but keep EDGE_DIST pixels off the edge
% Map min -> 3
% Map max -> N - 3
EDGE_DIST = n / 10;

scale = (n - (2 * EDGE_DIST)) / maxdiff;
t = EDGE_DIST - bbox_min * scale;

data = awf_translate_pts(d * scale, t);
model = awf_translate_pts(m * scale, t);

end