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

vx = p(1);
vy = p(2);
vz = p(3);
s = p(4);
tx = p(5);
ty = p(6);
tz = p(7);

dx = D(:,1);
dy = D(:,2);
dz = D(:,3);

s2 = s.*s;
vx2 = vx.*vx;
vy2 = vy.*vy;
vz2 = vz.*vz;

squares = s2 + vx2 + vy2 + vz2;
denom = 1 / squares;
denom2 = 1 / (squares .* squares);

dxtx = (dx + tx);
dyty = (dy + ty);
dztz = (dz + tz);

Tx_num = ((s2 + vx2 - vy2 - vz2) .* dxtx) + ...
         (2.0 .* (vx.*vy - s.*vz) .* dyty) + ...
         (2.0 .* (vx.*vz + s.*vy) .* dztz);
     
Ty_num = (2.0 .* (vy.*vx + s.*vz) .* dxtx) + ...
         ((s2 - vx2 + vy2 - vz2) .* dyty) + ...
         (2.0 .* (vy.*vz - s.*vx) .* dztz);

Tz_num = (2.0 .* (vz.*vx - s.*vy) .* dxtx) + ...
         (2.0 .* (vz.*vy + s.*vx) .* dyty) + ...
         ((s2 - vx2 - vy2 + vz2) .* dztz);

T(:,1) = (Tx_num .* denom);
     
T(:,2) = (Ty_num .* denom);
     
T(:,3) = (Tz_num .* denom);

% dTx/dvx
Jx(:, 1) = (((2.0 .* vx .* dxtx) + (2.0 .* vy .* dyty) + (2.0 .* vz .* dztz)) .* denom) - ...
           ((2.0 .* vx) .* Tx_num) .* denom2;
% dTx/dvy
Jx(:, 2) = ((-(2.0 .* vy .* dxtx) + (2.0 .* vx .* dyty) + (2.0 .* s .* dztz)) .* denom) - ...
           ((2.0 .* vy) .* Tx_num) .* denom2;
% dTx/dvz
Jx(:, 3) = ((-(2.0 .* vz .* dxtx) - (2.0 .* s .* dyty) + (2.0 .* vx .* dztz)) .* denom) - ...
           ((2.0 .* vz) .* Tx_num) .* denom2;   
% dTx/ds
Jx(:, 4) = (((2.0 .* s .* dxtx) - (2.0 .* vz .* dyty) + (2.0 .* vy .* dztz)) .* denom) - ...
           ((2.0 .* s) .* Tx_num) .* denom2;
% dTx/tx
Jx(:, 5) = (s2 + vx2 - vy2 - vz2) .* denom;
% dTx/ty
Jx(:, 6) = (2.0 .* (vx.*vy - s.*vz)) .* denom;
% dTx/tz
Jx(:, 7) = (2.0 .* (vx.*vz + s.*vy)) .* denom;

% dTy/dvx
Jy(:, 1) = (((2.0 .* vy .* dxtx) - (2.0 .* vx .* dyty) - (2.0 .* s .* dztz)) .* denom) - ...
           ((2.0 .* vx) .* Ty_num) .* denom2;
% dTy/dvy
Jy(:, 2) = (((2.0 .* vx .* dxtx) + (2.0 .* vy .* dyty) + (2.0 .* vz .* dztz)) .* denom) - ...
           ((2.0 .* vy) .* Ty_num) .* denom2;
% dTy/dvz
Jy(:, 3) = (((2.0 .* s .* dxtx) - (2.0 .* vz .* dyty) + (2.0 .* vy .* dztz)) .* denom) - ...
           ((2.0 .* vz) .* Ty_num) .* denom2;
% dTy/ds
Jy(:, 4) = (((2.0 .* vz .* dxtx) + (2.0 .* s .* dyty) - (2.0 .* vx .* dztz)) .* denom) - ...
           ((2.0 .* s) .* Ty_num) .* denom2;
% dTy/tx
Jy(:, 5) = (2.0 .* (vy.*vx + s.*vz)) .* denom;
% dTy/ty
Jy(:, 6) = (s2 - vx2 + vy2 - vz2) .* denom;
% dTy/tz
Jy(:, 7) = (2.0 .* (vy.*vz - s.*vx)) .* denom;

% dTz/dvx
Jz(:, 1) = (((2.0 .* vz .* dxtx) + (2.0 .* s .* dyty) - (2.0 .* vx .* dztz)) .* denom) - ...
           ((2.0 .* vx) .* Tz_num) .* denom2;
% dTz/dvy
Jz(:, 2) = ((-(2.0 .* s .* dxtx) + (2.0 .* vz .* dyty) - (2.0 .* vy .* dztz)) .* denom) - ...
           ((2.0 .* vy) .* Tz_num) .* denom2;
% dTz/dvz
Jz(:, 3) = (((2.0 .* vx .* dxtx) + (2.0 .* vy .* dyty) + (2.0 .* vz .* dztz)) .* denom) - ...
           ((2.0 .* vz) .* Tz_num) .* denom2;  
% dTz/ds
Jz(:, 4) = ((-(2.0 .* vy .* dxtx) + (2.0 .* vx .* dyty) + (2.0 .* s .* dztz)) .* denom) - ...
           ((2.0 .* s) .* Tz_num) .* denom2;
% dTz/tx
Jz(:, 5) = (2.0 .* (vz.*vx - s.*vy)) .* denom;
% dTz/ty
Jz(:, 6) = (2.0 .* (vz.*vy + s.*vx)) .* denom;
% dTz/tz
Jz(:, 7) = (s2 - vx2 - vy2 + vz2) .* denom;

end
