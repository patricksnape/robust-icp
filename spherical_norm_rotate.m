function [ rotated ] = spherical_norm_rotate(phi, theta, data)
%spherical_norm_rotate Given an Azimuth and Elevation rotate a point cloud
%   Rotate a point cloud given a pair of spherical coordinates representing
%   a rotation amount

    data_centre = mean(data);
    rotated = awf_translate_pts(data, -data_centre);

    [a, e, r] = cart2sph(rotated(:, 1), rotated(:, 2), rotated(:, 3));

    a = a + phi;
    e = e + theta;

    [x, y, z] = sph2cart(a, e, r);
    rotated = [x, y, z];

    rotated = awf_translate_pts(rotated, data_centre);

end

