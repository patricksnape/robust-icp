function [a v] = angleaxis(x)
% angleaxis quaternion direction cosine matrix angle axis
%*******************************************************************
%
% angleaxis calculates the rotation angle and rotation axis of the
% input quaternion or direction cosine matrix.
%
% Input: x = quaternion, x(1) = scalar, x(2:4) = vector
% Rotation sense = Successive rotations are right multiplies.
% Assumes x is normalized.
%
% or
%
% x = direction cosine matrix.
% Assumes x is orthonormalized.
%
% Output: a = rotation angle (radians)
% v = rotation axis (1x3 unit vector)
%
% Programmer: James Tursa
%
%*******************************************************************

if( numel(x) == 9 )
    q = dc2quat(x);
else
    q = x;
end
a = 2 * acos(q(1));
if( nargout == 2 )
    if( a == 0 )
        v = [1 0 0];
    else
        v = q(2:4)/sin(a/2);
    end
end

return
end