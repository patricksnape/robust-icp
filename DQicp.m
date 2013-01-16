function [T] = DQicp(M, D, varargin)
% DQICP Closed-form dual quaternion ICP
%
% Walker, M. W., Shao, L. & Volz, R. A. (1991) 
% Estimating 3-D location parameters using dual number quaternions. 
% CVGIP: Image Understanding. 54 (3), 358-367.
%
% [T] = DQicp(M,D)   returns the homogenous transformation matrix T
% that minimizes the distances from (T * p + TT) to q.
%
% [T] = DQicp(M, D, Mv, Dv)   uses the two sets of feature vectors Mv and 
% Dv in conjunction with the point data
%
% Author: Patrick Snape
% Date: 16 Jan 2013

%% Parse input
inp = inputParser;

% Model point cloud - 3xK
inp.addRequired('M', @(x)isreal(x) && size(x,1) == 3);
% Data point cloud - 3xK
inp.addRequired('D', @(x)isreal(x) && size(x,1) == 3);
% Model feature vectors - 3xL
inp.addOptional('Mv', 10, @(x)isreal(x) && size(x,1) == 3);
% Data feature vectors - 3xL
inp.addOptional('Dv', 10, @(x)isreal(x) && size(x,1) == 3);

inp.parse(M, D, varargin{:});
arg = inp.Results;
clear('inp');

%% Convert to position quaternion

% Number of points - assume size(D) == size(M)
k = size(M, 2);

M = 0.5 * [arg.M; zeros(1, k)];
D = 0.5 * [arg.D; zeros(1, k)];

%% Closed-form solution

C1 = zeros(4,4);
for i=1:k
    C1 = C1 + Q(D(:, i))' * W(M(:, i));
end
C1 = C1 * -2;

% Assumes all weight = 1 therefore sum(weights) = number of points
C2 = k * eye(4);

C3 = zeros(4,4);
for i=1:k
    C3 = C3 + (W(M(:, i)) - Q(D(:, i)));
end
C3 = C3 * 2;

A = 0.5 * (C3' * inv(C2 + C2') * C3 - C1 - C1');

%% Results and build homogenous transformation matrix, T

[r, ~] = eigs(A, 1);

s = -inv(C2 + C2') * C3 * r;

R = coolquat2mat(r);

t = W(r)' * s;
t = 0.5 * t(1:3);

T = [ 
      R     t;
      0 0 0 1
    ];

end

% Builds a skew matrix from quaternion q = [ vx vy vz s ]'
function skew = K(q)
    vx = q(1);
    vy = q(2);
    vz = q(3);
    
    skew = [ 0    -vz   vy;
             vz    0    vx;
            -vy    vx   0 
           ];
end

% Builds the Q matrix as defined in (15) from real part r = [ vx vy vz s ]'
function q = Q(r)
    s = r(4);
    
    q = zeros(4, 4);
    q(1:3, 1:3) = s * eye(3,3) + K(r);
    q(:, 4) = r;
    q(4, :) = (-r)';
end

% Builds the W matrix as defined in (16) from real part r = [ vx vy vz s ]'
function q = W(r)
    s = r(4);
    
    q = zeros(4, 4);
    q(1:3, 1:3) = s * eye(3,3) - K(r);
    q(:, 4) = r;
    q(4, :) = (-r)';
end