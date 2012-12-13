function [R] = spherical_rotation(az, el)

[x, y, z] = sph2cart(0, 0, 1);
X = [x, y, z];
[x, y, z] = sph2cart(az, el, 1);
Y = [x, y, z];

C = cross(X, Y);
CX = cross(C, X);

R = [X/norm(X); C/norm(C); CX/norm(CX)];

end