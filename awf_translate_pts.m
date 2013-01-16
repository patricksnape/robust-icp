function xpt = awf_translate_pts(x, t)
% AWF_TRANSLATE_PTS translate matrix of vectors by t
%
% Author: Patrick Snape
% Date: 16 Jan 2013

[r, c] = size(t);

if r == 3 && c == 1
    xpt = x + repmat(t, 1, size(x, 2));
elseif c == 3 && r == 1
    xpt = x + repmat(t, size(x, 1), 1);
else
    fprintf('Failed to translate points, t must be either 1x3 or 3x1');
end