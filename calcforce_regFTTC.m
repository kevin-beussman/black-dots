% This program (regularized fourier transform traction force
% reconstruction) was produced at the University of Heidelberg, BIOMS
% group of Ulrich Schwarz. It calculates traction from a gel displacement
% field.
%
% Benedikt Sabass 13-10-2008
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
%
% This file is part of TFM_Package.
% 
% TFM_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TFM_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TFM_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
function  [trac,resid_norm,soln_norm,disp_e] = calcforce_regFTTC(displacement,E,nu,dot_spacing,lambdas,real_points)
% added by Achim:
% Input : grid_mat, u, cluster_size have to be in the same units, namely
%         pixels. If they are given in different units e.g. meters, then a
%         different scaling factor for the elastic energy has to be
%         applied! The force value remains, however, unaffected, see below.
% Output: The output force is actually a surface stress with the same units
%         as the input E! In particular, the unit of the output force is
%         independent of the units of the input grid_mat,u and cluster_size
%         The reason for this is essentially that the elastic stress is
%         only dependent on the non-dimensional strain which is given by
%         spatial derivatives of the displacements, that is du/dx. If u and
%         dx (essentially cluster_size) are in the same units, then the
%         resulting force has the same dimension as the input E.
% updated by Sangyoon Han for usage for L1 regularization
    
npx = 2 + mod(size(displacement(:,:,1),2),2);
npy = 2 + mod(size(displacement(:,:,1),1),2);
zero_pad = zeros([size(displacement(:,:,1)) + [npy,npx], 2]);
u = zero_pad;
u(2:end-npy+1,2:end-npx+1,1) = displacement(:,:,1);
u(2:end-npy+1,2:end-npx+1,2) = displacement(:,:,2);
[i_max,j_max] = size(u(:,:,1));

V = 2*(1+nu)/E;

kx_vec = 2*pi/i_max/dot_spacing.*[0:(i_max/2-1) (-i_max/2:-1)];
ky_vec = 2*pi/j_max/dot_spacing.*[0:(j_max/2-1) (-j_max/2:-1)];
kx = repmat(kx_vec',1,j_max);
ky = repmat(ky_vec,i_max,1);
if nargin<10
    Rx=ones(size(kx));
    Ry=ones(size(ky));
end


kx(1,1) = 1;
ky(1,1) = 1;
X = i_max*dot_spacing/2;
Y = j_max*dot_spacing/2;
g0x = pi.^(-1).*V.*((-1).*Y.*log((-1).*X+sqrt(X.^2+Y.^2))+Y.*log( ...
      X+sqrt(X.^2+Y.^2))+((-1)+nu).*X.*(log((-1).*Y+sqrt(X.^2+Y.^2) ...
      )+(-1).*log(Y+sqrt(X.^2+Y.^2))));
g0y = pi.^(-1).*V.*(((-1)+nu).*Y.*(log((-1).*X+sqrt(X.^2+Y.^2))+( ...
      -1).*log(X+sqrt(X.^2+Y.^2)))+X.*((-1).*log((-1).*Y+sqrt( ...
      X.^2+Y.^2))+log(Y+sqrt(X.^2+Y.^2))));

for i_lambda = 1:length(lambdas)
    L = lambdas(i_lambda);
    
    % inverse solution
    Ginv_xx = (kx.^2+ky.^2).^(-1/2).*V.*(kx.^2.*L.*Rx+ky.^2.*L.*Ry+V.^2).^(-1).*(kx.^2.* ...
              L.*Rx+ky.^2.*L.*Ry+((-1)+nu).^2.*V.^2).^(-1).*(kx.^4.*(L.*Rx+(-1).*L.*nu.*Rx)+ ...
              kx.^2.*((-1).*ky.^2.*L.*Ry.*((-2)+nu)+(-1).*((-1)+nu).*V.^2)+ky.^2.*( ...
              ky.^2.*L.*Ry+((-1)+nu).^2.*V.^2));
    Ginv_yy = (kx.^2+ky.^2).^(-1/2).*V.*(kx.^2.*L+ky.^2.*L.*Ry+V.^2).^(-1).*(kx.^2.* ...
              L.*Rx+ky.^2.*L.*Ry+((-1)+nu).^2.*V.^2).^(-1).*(kx.^4.*L+(-1).*ky.^2.*((-1)+ ...
              nu).*(ky.^2.*L.*Rx+V.^2)+kx.^2.*((-1).*ky.^2.*L.*Ry.*((-2)+nu)+((-1)+nu).^2.* ...
              V.^2));
    Ginv_xy = (-1).*kx.*ky.*(kx.^2+ky.^2).^(-1/2).*nu.*V.*(kx.^2.*L.*Rx+ky.^2.*L.*Ry+ ...
              V.^2).^(-1).*(kx.^2.*L.*Rx+ky.^2.*L.*Ry+((-1)+nu).*V.^2).*(kx.^2.*L.*Rx+ky.^2.* ...
              L.*Ry+((-1)+nu).^2.*V.^2).^(-1);

    Ginv_xx(1,1) = 1/g0x; % effectively, you can just set this equal to 0
    Ginv_yy(1,1) = 1/g0y;
    Ginv_xy(1,1) = 0;

    Ginv_xy(i_max/2+1,:) = 0;
    Ginv_xy(:,j_max/2+1) = 0;

    Ftu(:,:,1) = fft2(u(:,:,1));
    Ftu(:,:,2) = fft2(u(:,:,2));

    Ftf(:,:,1) = Ginv_xx.*Ftu(:,:,1) + Ginv_xy.*Ftu(:,:,2);
    Ftf(:,:,2) = Ginv_xy.*Ftu(:,:,1) + Ginv_yy.*Ftu(:,:,2);

    f(:,:,1) = ifft2(Ftf(:,:,1),'symmetric');
    f(:,:,2) = ifft2(Ftf(:,:,2),'symmetric');

    % remove padding
    trac_x = f(2:end-npy+1,2:end-npx+1,1);
    trac_y = f(2:end-npy+1,2:end-npx+1,2);
    trac(:,:,1,i_lambda) = trac_x;
    trac(:,:,2,i_lambda) = trac_y;
    % UNITS are traction stress, same as E (pN*um^-2)
    
    soln_norm(i_lambda) = sum(sqrt(sum([trac_x(:); trac_y(:)].^2,2)));

    % forward solution
%     kx = repmat(kx_vec',1,j_max);
%     ky = repmat(ky_vec,i_max,1);
    k = sqrt(kx.^2+ky.^2);
    
    G_xx = 2*(1 + nu)./(E*k.^3).*((1 - nu)*k.^2 + nu*ky.^2);
    G_yy = 2*(1 + nu)./(E*k.^3).*((1 - nu)*k.^2 + nu*kx.^2);
    G_xy = 2*(1 + nu)./(E*k.^3).*(-nu*kx.*ky);

    Ftu_e(:,:,1) = G_xx.*Ftf(:,:,1) + G_xy.*Ftf(:,:,2);
    Ftu_e(:,:,2) = G_xy.*Ftf(:,:,1) + G_yy.*Ftf(:,:,2);

    u_e(:,:,1) = ifft2(Ftu_e(:,:,1),'symmetric');
    u_e(:,:,2) = ifft2(Ftu_e(:,:,2),'symmetric');

    % remove padding
    disp_e(:,:,:,i_lambda) = u_e(2:end-npy+1,2:end-npx+1,:);
    % UNITS are displacement (same as input displacement)
    
    resid_x = disp_e(:,:,1,i_lambda) - displacement(:,:,1);
    resid_y = disp_e(:,:,2,i_lambda) - displacement(:,:,2);
    resid = sqrt(resid_x.^2 + resid_y.^2);
    resid(~real_points) = 0;
%     max(resid)
%     min(resid)
    resid_norm(i_lambda) = sum(resid(real_points));
%     if i_lambda == 100
%         resid_norm(i_lambda)
%         1;
%     end
%     subplot(1,2,1); imagesc(displacement(:,:,1)); subplot(1,2,2); imagesc(disp_e(:,:,1,i_lambda));

%     fnorm = (force(:,1).^2 + force(:,2).^2).^0.5;
%     energie = 1/2*sum(sum(u(2:end-1,2:end-1,1).*f(2:end-1,2:end-1,1) + u(2:end-1,2:end-1,2).*f(2:end-1,2:end-1,2)))*(cluster_size)^2*pix_durch_my^3/10^6; 
end