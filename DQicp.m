function [T, E] = DQicp(M, D, varargin)
% DQICP Closed-form dual quaternion ICP
%
% Input: 
%       M   (required) 3xK matrix representing a cloud of points
%       D   (required) 3xK matrix representing a cloud of points
%       Mv  (optional) 3xL matrix representing feature vectors for the
%           point cloud M. Defaults to [0; 0; 0; 1]xK
%       Dv  (optional) 3xL matrix representing feature vectors for the
%           point cloud D. Defaults to [0; 0; 0; 1]xK
% Output:
%       T   returns the homogenous transformation matrix T
%           that minimizes the distances from (T * D) to M.
%       E   the total error ||T*D - M||
%
% Walker, M. W., Shao, L. & Volz, R. A. (1991)
% Estimating 3-D location parameters using dual number quaternions.
% CVGIP: Image Understanding. 54 (3), 358-367.
%
% Author: Patrick Snape
% Date: 16 Jan 2013

%% Parse input
inp = inputParser;


inp.addRequired('M', @(x)isreal(x) && size(x,1) == 3);
inp.addRequired('D', @(x)isreal(x) && size(x,1) == 3);

% Default feature vectors
default_vectors = repmat([0; 0; 0], 1, size(M, 2));

inp.addOptional('Mv', default_vectors, @(x)isreal(x) && size(x,1) == 3);
inp.addOptional('Dv', default_vectors, @(x)isreal(x) && size(x,1) == 3);

inp.parse(M, D, varargin{:});
arg = inp.Results;
clear('inp');

%% Convert to quaternions

% Number of points - assume size(D) == size(M)
k = size(arg.M, 2);

M = 0.5 * [arg.M; zeros(1, k)];
D = 0.5 * [arg.D; zeros(1, k)];

% Number of feature vectors - assume size(Dv) == size(Mv)
l = size(arg.Mv, 2);

Mv = [arg.Mv; zeros(1, l)];
Dv = [arg.Dv; zeros(1, l)];

%% Closed-form solution
%   NOTE: We have swapped the order of M and D from the paper as the paper
%   calculates D = R*M and we want M = R*D to be consistent with other ICP
%   methods

C1 = zeros(4,4);
% Points
for i=1:k
    C1 = C1 + Q(M(:, i))' * W(D(:, i));
end
% Feature Vectors
for i=1:l
    C1 = C1 + Q(Mv(:, i))' * W(Dv(:, i));
end
C1 = C1 * -2;

% Assumes all weight = 1 therefore sum(weights) = number of points
C2 = k * eye(4);

C3 = zeros(4,4);
for i=1:k
    C3 = C3 + (W(D(:, i)) - Q(M(:, i)));
end
C3 = C3 * 2;

A = 0.5 * (C3' * inv(C2 + C2') * C3 - C1 - C1');

%% Results and build homogenous transformation matrix, T

[r, ~] = eigs(A, 1);

s = -inv(C2 + C2') * C3 * r;
R = quat2rot(r);

t = W(r)' * s;
t = 0.5 * t(1:3);

T = [
        R     t;
        0 0 0 1
    ];

E = norm(T * D - M);

end

%% Builds a skew matrix from quaternion q = [ vx vy vz s ]'
function skew = K(q)
vx = q(1);
vy = q(2);
vz = q(3);

skew = [ 
            0    -vz   vy;
            vz    0   -vx;
           -vy    vx   0
       ];
end

%% Builds the Q matrix as defined in (15) from real part r = [ vx vy vz s ]'
function q = Q(r)
s = r(4);

q = zeros(4, 4);
q(1:3, 1:3) = s * eye(3,3) + K(r);
q(:, 4) = r;
q(4, :) = (-r)';
end

%% Builds the W matrix as defined in (16) from real part r = [ vx vy vz s ]'
function q = W(r)
s = r(4);

q = zeros(4, 4);
q(1:3, 1:3) = s * eye(3,3) - K(r);
q(:, 4) = r;
q(4, :) = (-r)';
end