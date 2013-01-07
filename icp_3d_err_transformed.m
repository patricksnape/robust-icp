function [T, Jx, Jy, Jz] = icp_3d_err_transformed(p, D)

% ICP_3D_ERR    Error function for icp3d
% In:
%               p 7x1 params
%               D Nx3 data points
% Out:
%               T Nx3 transformed data points R(q(p1..4)) * (D + p(5:7))
%               Jx Nx7 jacobian dT(1)/dp
%               Jy Nx7 jacobian dT(2)/dp
%               Jz Nx7 jacobian dT(3)/dp

% Author: Patrick Snape <pts08@ic.ac.uk>
% Date: 07 Jan 13

s = p(1);
vx = p(2);
vy = p(3);
vz = p(4);
tx = p(5);
ty = p(6);
tz = p(7);

dx = D(:,1);
dy = D(:,2);
dz = D(:,3);

vz2 = vz.*vz;
s2 = s.*s;
vx2 = vx.*vx;
vy2 = vy.*vy;

% My failed attempt at a derivation from Castellani and Bartoli

% denom = s2 + vx2 + vy2 + vz2;
% 
% T(:,1) = ((s2 + vx2 - vy2 - vz2) .* dx) / denom + ...
%          (2.0 * (vx.*vy - s.*vz) .* dy) / denom + ...
%          (2.0 * (vx.*vz - s.*vy) .* dz) / denom + ...
%          tx;
%      
% T(:,2) = (2.0 * (vy.*vx + s.*vz) .* dx) / denom + ...
%          ((s2 - vx2 + vy2 - vz2) .* dy) / denom + ...
%          (2.0 * (vy.*vz - s.*vx) .* dz) / denom + ...
%          ty;
%      
% T(:,3) = (2.0 * (vz.*vx + s.*vy) .* dx) / denom + ...
%          (2.0 * (vz.*vy + s.*vx) .* dy) / denom + ...
%          ((s2 - vx2 - vy2 + vz2) .* dz) / denom + ...
%          tz;
% 
% % dTx/ds
% Jx(:, 1) = ((((2.0 .* s .* dx) - (2.0 .* vz .* dy) - (2.0 .* vy .* dz)) .* denom) - ...
%            ((2.0 .* s) .* ((s2 + vx2 - vy2 - vz2) .* dx) + (2.0 * (vx.*vy - s.*vz) .* dy) + (2.0 * (vx.*vz - s.*vy) .* dz))) / ...
%            denom.^2;
% % dTx/dvx
% Jx(:, 2) = ((((2.0 .* vx .* dx) + (2.0 .* vy .* dy) + (2.0 .* vz .* dz)) .* denom) - ...
%            ((2.0 .* vx) .* ((s2 + vx2 - vy2 - vz2) .* dx) + (2.0 * (vx.*vy - s.*vz) .* dy) + (2.0 * (vx.*vz - s.*vy) .* dz))) / ...
%            denom.^2;
% % dTx/dvy
% Jx(:, 3) = ((((2.0 .* vy .* dx) + (2.0 .* vx .* dy) - (2.0 .* s .* dz)) .* denom) - ...
%            ((2.0 .* vy) .* ((s2 + vx2 - vy2 - vz2) .* dx) + (2.0 * (vx.*vy - s.*vz) .* dy) + (2.0 * (vx.*vz - s.*vy) .* dz))) / ...
%            denom.^2;
% % dTx/dvz
% Jx(:, 4) = ((((2.0 .* vz .* dx) - (2.0 .* s .* dy) + (2.0 .* vx .* dz)) .* denom) - ...
%            ((2.0 .* vz) .* ((s2 + vx2 - vy2 - vz2) .* dx) + (2.0 * (vx.*vy - s.*vz) .* dy) + (2.0 * (vx.*vz - s.*vy) .* dz))) / ...
%            denom.^2;
% % dTx/tx
% Jx(:, 5) = 1;
% % dTx/ty
% Jx(:, 6) = 0;
% % dTx/tz
% Jx(:, 7) = 0;
% 
% % dTy/ds
% Jy(:, 1) = ((((2.0 .* vz .* dx) + (2.0 .* s .* dy) - (2.0 .* vx .* dz)) .* denom) - ...
%            ((2.0 .* s) .* (2.0 * (vy.*vx + s.*vz) .* dx) + ((s2 - vx2 + vy2 - vz2) .* dy) + (2.0 * (vy.*vz - s.*vx) .* dz))) / ...
%            denom.^2;
% % dTy/dvx
% Jy(:, 2) = ((((2.0 .* vy .* dx) - (2.0 .* vx .* dy) - (2.0 .* s .* dz)) .* denom) - ...
%            ((2.0 .* dx) .* (2.0 * (vy.*vx + s.*vz) .* dx) + ((s2 - vx2 + vy2 - vz2) .* dy) + (2.0 * (vy.*vz - s.*vx) .* dz))) / ...
%            denom.^2;
% % dTy/dvy
% Jy(:, 3) = ((((2.0 .* vx .* dx) + (2.0 .* vy .* dy) + (2.0 .* vz .* dz)) .* denom) - ...
%            ((2.0 .* dy) .* (2.0 * (vy.*vx + s.*vz) .* dx) + ((s2 - vx2 + vy2 - vz2) .* dy) + (2.0 * (vy.*vz - s.*vx) .* dz))) / ...
%            denom.^2;
% % dTy/dvz
% Jy(:, 4) = ((((2.0 .* s .* dx) - (2.0 .* vz .* dy) + (2.0 .* vy .* dz)) .* denom) - ...
%            ((2.0 .* dz) .* (2.0 * (vy.*vx + s.*vz) .* dx) + ((s2 - vx2 + vy2 - vz2) .* dy) + (2.0 * (vy.*vz - s.*vx) .* dz))) / ...
%            denom.^2;
% % dTy/tx
% Jy(:, 5) = 0;
% % dTy/ty
% Jy(:, 6) = 1;
% % dTy/tz
% Jy(:, 7) = 0;
% 
% % dTz/ds
% Jz(:, 1) = ((((2.0 .* vy .* dx) + (2.0 .* vx .* dy) + (2.0 .* s .* dz)) .* denom) - ...
%            ((2.0 .* s) .* (2.0 * (vz.*vx + s.*vy) .* dx) + (2.0 * (vz.*vy + s.*vx) .* dy) + ((s2 - vx2 - vy2 + vz2) .* dz))) / ...
%            denom.^2;
% % dTz/dvx
% Jz(:, 2) = ((((2.0 .* vz .* dx) + (2.0 .* s .* dy) - (2.0 .* vx .* dz)) .* denom) - ...
%            ((2.0 .* vx) .* (2.0 * (vz.*vx + s.*vy) .* dx) + (2.0 * (vz.*vy + s.*vx) .* dy) + ((s2 - vx2 - vy2 + vz2) .* dz))) / ...
%            denom.^2;
% % dTz/dvy
% Jz(:, 3) = ((((2.0 .* s .* dx) + (2.0 .* vz .* dy) - (2.0 .* vy .* dz)) .* denom) - ...
%            ((2.0 .* vy) .* (2.0 * (vz.*vx + s.*vy) .* dx) + (2.0 * (vz.*vy + s.*vx) .* dy) + ((s2 - vx2 - vy2 + vz2) .* dz))) / ...
%            denom.^2;
% % dTz/dvz
% Jz(:, 4) = ((((2.0 .* vx .* dx) + (2.0 .* vy .* dy) + (2.0 .* vz .* dz)) .* denom) - ...
%            ((2.0 .* vz) .* (2.0 * (vz.*vx + s.*vy) .* dx) + (2.0 * (vz.*vy + s.*vx) .* dy) + ((s2 - vx2 - vy2 + vz2) .* dz))) / ...
%            denom.^2;
% % dTz/tx
% Jz(:, 5) = 0;
% % dTz/ty
% Jz(:, 6) = 0;
% % dTz/tz
% Jz(:, 7) = 1;

% Andrew Fitzgibbon's derivation

t5 = s2+vz2-vx2-vy2;
t6 = s2+vx2+vy2+vz2;
t7 = 1/t6;
t8 = t5.*t7;
t9 = dx+tx;
t11 = s.*vx;
t12 = vz.*vy;
t13 = t11+t12;
t14 = 2.0.*t13.*t7;
t15 = dy+ty;
t17 = vy.*s;
t18 = vz.*vx;
t19 = t17-t18;
t20 = 2.0.*t19.*t7;
t21 = dz+tz;
t24 = s.*t7;
t25 = t24.*t9;
t26 = t6.*t6;
t27 = 1/t26;
t28 = t5.*t27;
t29 = t9.*s;
t31 = vx.*t7;
t32 = t31.*t15;
t33 = 2.0.*t13.*t27;
t34 = t15.*s;
t36 = vy.*t7;
t37 = t36.*t21;
t38 = 2.0.*t19.*t27;
t39 = t21.*s;
t42 = t31.*t9;
t43 = t9.*vx;
t45 = t24.*t15;
t46 = t15.*vx;
t48 = vz.*t7;
t49 = t48.*t21;
t50 = t21.*vx;
t53 = t36.*t9;
t54 = t9.*vy;
t56 = t48.*t15;
t57 = t15.*vy;
t59 = t24.*t21;
t60 = t21.*vy;
t63 = t48.*t9;
t64 = t9.*vz;
t66 = t36.*t15;
t67 = t15.*vz;
t69 = t31.*t21;
t70 = t21.*vz;
t73 = t11-t12;
t74 = 2.0.*t73.*t7;
t76 = vz2-s2+vx2-vy2;
t77 = t76.*t7;
t79 = vx.*vy;
t80 = vz.*s;
t81 = t79+t80;
t82 = 2.0.*t81.*t7;
t85 = 2.0.*t73.*t27;
t87 = t76.*t27;
t89 = 2.0.*t81.*t27;
t104 = t17+t18;
t105 = 2.0.*t104.*t7;
t107 = t79-t80;
t108 = 2.0.*t107.*t7;
t110 = vz2-s2-vx2+vy2;
t111 = t110.*t7;
t114 = 2.0.*t104.*t27;
t116 = 2.0.*t107.*t27;
t118 = t110.*t27;

T(:,1) = t8.*t9+t14.*t15+t20.*t21;
T(:,2) = t74.*t9+t77.*t15+t82.*t21;
T(:,3) = t105.*t9+t108.*t15+t111.*t21;

Jx(:,1) = 2.0.*t25-2.0.*t28.*t29+2.0.*t32-2.0.*t33.*t34+2.0.*t37-2.0.*t38.*t39;
Jx(:,2) = -2.0.*t42-2.0.*t28.*t43+2.0.*t45-2.0.*t33.*t46-2.0.*t49-2.0.*t38.*t50;
Jx(:,3) = -2.0.*t53-2.0.*t28.*t54+2.0.*t56-2.0.*t33.*t57+2.0.*t59-2.0.*t38.*t60;
Jx(:,4) = 2.0.*t63-2.0.*t28.*t64+2.0.*t66-2.0.*t33.*t67-2.0.*t69-2.0.*t38.*t70;
Jx(:,5) = t8;
Jx(:,6) = t14;
Jx(:,7) = t20;

Jy(:,1) = 2.0.*t42-2.0.*t85.*t29-2.0.*t45-2.0.*t87.*t34+2.0.*t49-2.0.*t89.*t39;
Jy(:,2) = 2.0.*t25-2.0.*t85.*t43+2.0.*t32-2.0.*t87.*t46+2.0.*t37-2.0.*t89.*t50;
Jy(:,3) = -2.0.*t63-2.0.*t85.*t54-2.0.*t66-2.0.*t87.*t57+2.0.*t69-2.0.*t89.*t60;
Jy(:,4) = -2.0.*t53-2.0.*t85.*t64+2.0.*t56-2.0.*t87.*t67+2.0.*t59-2.0.*t89.*t70;
Jy(:,5) = t74;
Jy(:,6) = t77;
Jy(:,7) = t82;

Jz(:,1) = 2.0.*t53-2.0.*t114.*t29-2.0.*t56-2.0.*t116.*t34-2.0.*t59-2.0.*t118.*t39;
Jz(:,2) = 2.0.*t63-2.0.*t114.*t43+2.0.*t66-2.0.*t116.*t46-2.0.*t69-2.0.*t118.*t50;
Jz(:,3) = 2.0.*t25-2.0.*t114.*t54+2.0.*t32-2.0.*t116.*t57+2.0.*t37-2.0.*t118.*t60;
Jz(:,4) = 2.0.*t42-2.0.*t114.*t64-2.0.*t45-2.0.*t116.*t67+2.0.*t49-2.0.*t118.*t70;
Jz(:,5) = t105;
Jz(:,6) = t108;
Jz(:,7) = t111;

end
