function T = DQicp(M, D)
%UNTITLED4 Dual quaternion ICP
%   Detailed explanation goes here

%% Convert to position quaternion

M = 0.5 * [M; zeros(1, size(M, 2))];
D = 0.5 * [D; zeros(1, size(D, 2))];

% Number of points
n = size(M, 2);

%% Closed-form solution

C1 = zeros(4,4);
for i=1:n
    C1 = C1 + Q(D(:, i))' * W(M(:, i));
end
C1 = C1 * -2;

% Assumes all weight = 1 therefore sum(weights) = number of points
C2 = n * eye(4);

C3 = zeros(4,4);
for i=1:n
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