function [R] = spherical_rotation(az, el)

R = [sin(el)*cos(az) -sin(az) cos(el)*cos(az);
     sin(el)*sin(az) cos(az)  cos(el)*sin(az);
     -cos(el)        0        sin(el)];

end