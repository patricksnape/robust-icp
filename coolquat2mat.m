function R = coolquat2mat(q)

% COOLQUAT2MAT  Convert quaternion to 3x3 rotation matrix
% Input: q = [ vx vy vz s ]
%
% Output: R =
%     [ s2 + x2 - y2 - z2    2 * (xy - sz)        2 * (zx + sy)     ]
%     [ 2 * (xy + sz)        s2 - x2 + y2 - z2    2 * (yz - sx)     ]
%     [ 2 * (zx - sy)        2 * (yz + sx)        s2 - x2 - y2 + z2 ]

% Author: Andrew Fitzgibbon <awf@robots.ox.ac.uk>
% Date: 23 Mar 99
%
% Edited: Patrick Snape
% Date: 16 Jan 13

% Cache squares
x2 = q(1) * q(1);
y2 = q(2) * q(2);
z2 = q(3) * q(3);
s2 = q(4) * q(4);

% Fill diagonal terms
R(1, 1) = s2 + x2 - y2 - z2;		
R(2, 2) = s2 - x2 + y2 - z2;
R(3, 3) = s2 - x2 - y2 + z2;

% Cache multiplications
xy = q(1) * q(2);
yz = q(2) * q(3);
zx = q(3) * q(1);
sx = q(4) * q(1);
sy = q(4) * q(2);
sz = q(4) * q(3);

% Fill off diagonal terms
% x terms
R(1, 2) = 2 * (xy - sz);			
R(1, 3) = 2 * (zx + sy);

% y terms
R(2, 1) = 2 * (xy + sz);
R(2, 3) = 2 * (yz - sx);

% z terms
R(3, 1) = 2 * (zx - sy);
R(3, 2) = 2 * (yz + sx);
