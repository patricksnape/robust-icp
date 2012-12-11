function dc = quat2dc(q)
% quat2dc quaternion direction cosine matrix angle axis
%*******************************************************************
%
% quat2dc calculates the dirction cosine matrix corresponding to a
% quaternion. Assumes input quaternion is normalized.
%
% Input: q = quaternion, q(1) = scalar, q(2:4) = vector
% Rotation sense = Successive rotations are right multiplies.
%
% Output: dc = 3x3 direction cosine matrix
%
% Programmer: James Tursa
%
%*******************************************************************

q11 = q(1)^2;
q12 = q(1)*q(2);
q13 = q(1)*q(3);
q14 = q(1)*q(4);
q22 = q(2)^2;
q23 = q(2)*q(3);
q24 = q(2)*q(4);
q33 = q(3)^2;
q34 = q(3)*q(4);
q44 = q(4)^2;

dc = zeros(3);
dc(1,1) = q11 + q22 - q33 - q44;
dc(2,1) = 2 * (q23 - q14);
dc(3,1) = 2 * (q24 + q13);
dc(1,2) = 2 * (q23 + q14);
dc(2,2) = q11 - q22 + q33 - q44;
dc(3,2) = 2 * (q34 - q12);
dc(1,3) = 2 * (q24 - q13);
dc(2,3) = 2 * (q34 + q12);
dc(3,3) = q11 - q22 - q33 + q44;

return
end