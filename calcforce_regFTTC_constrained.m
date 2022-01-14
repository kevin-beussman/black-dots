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
function  [trac,resid_norm,soln_norm,disp_e] = calcforce_regFTTC_constrained(displacement,E,nu,dot_spacing,lambdas,real_points,cell_points,positions)
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

cp = zero_pad;
cp(2:end-npy+1,2:end-npx+1,1) = cell_points;
cp(2:end-npy+1,2:end-npx+1,2) = cell_points;
cp = logical(cp);

d = zero_pad;
d(2:end-npy+1,2:end-npx+1,1) = positions(:,:,1);
d(2:end-npy+1,2:end-npx+1,2) = positions(:,:,2);
x = d(:,:,1);
y = d(:,:,2);

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
    
    % inverse solution:
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

    Ginv_xx(1,1) = 1/g0x;
    Ginv_yy(1,1) = 1/g0y;
    Ginv_xy(1,1) = 0;

    Ginv_xy(i_max/2+1,:) = 0;
    Ginv_xy(:,j_max/2+1) = 0;
    
    % forward solution:
    k = sqrt(kx.^2+ky.^2);
    
    G_xx = 2*(1 + nu)./(E*k.^3).*((1 - nu)*k.^2 + nu*ky.^2);
    G_yy = 2*(1 + nu)./(E*k.^3).*((1 - nu)*k.^2 + nu*kx.^2);
    G_xy = 2*(1 + nu)./(E*k.^3).*(-nu*kx.*ky);
    
    % constrained FTTC iterative procedure
    n_iter = 0;
    itermax = 1000;
    u_constrained = u;
    while true
        n_iter = n_iter + 1;

        Ftu(:,:,1) = fft2(u_constrained(:,:,1));
        Ftu(:,:,2) = fft2(u_constrained(:,:,2));

        Ftf(:,:,1) = Ginv_xx.*Ftu(:,:,1) + Ginv_xy.*Ftu(:,:,2);
        Ftf(:,:,2) = Ginv_xy.*Ftu(:,:,1) + Ginv_yy.*Ftu(:,:,2);

        f_constrained(:,:,1) = ifft2(Ftf(:,:,1),'symmetric');
        f_constrained(:,:,2) = ifft2(Ftf(:,:,2),'symmetric');
        
        if n_iter == itermax
            break
%         elseif n_iter == (itermax - 1)
%             figure(1)
%             quiver(x,y,u(:,:,1),u(:,:,2),1,'-k','linewidth',2)
%             hold on
%             quiver(x,y,u_constrained(:,:,1),u_constrained(:,:,2),1,'-r')
%             hold off
%             title(num2str(n_iter))
%             drawnow
% %             pause(0.1)
        end
        
        if n_iter == 1
            f_init = f_constrained;
        end
        
        % step 2: set f = 0 outside the cell
        f_constrained(~cp) = 0;
        
        % step 3: forward solution
        Ftf_2(:,:,1) = fft2(f_constrained(:,:,1));
        Ftf_2(:,:,2) = fft2(f_constrained(:,:,2));
        
        Ftu_2(:,:,1) = G_xx.*Ftf_2(:,:,1) + G_xy.*Ftf_2(:,:,2);
        Ftu_2(:,:,2) = G_xy.*Ftf_2(:,:,1) + G_yy.*Ftf_2(:,:,2);
        
        u_constrained(:,:,1) = ifft2(Ftu_2(:,:,1),'symmetric');
        u_constrained(:,:,2) = ifft2(Ftu_2(:,:,2),'symmetric');
        
        % step 4: set u = meas inside the cell
        u_constrained(cp) = u(cp);
    end
    
    f = f_constrained;
    
    figure(1)
    quiver(x,y,u(:,:,1),u(:,:,2),1,'-k','linewidth',2)
    hold on
    quiver(x,y,u_constrained(:,:,1),u_constrained(:,:,2),1,'-r')
    hold off
    axis image
    title(['Displacement ' num2str(n_iter)])
    drawnow
    
    figure(2)
%     quiver(x,y,f_init(:,:,1),f_init(:,:,2),1,'-k','linewidth',2)
    hold on
    quiver(x,y,f_constrained(:,:,1),f_constrained(:,:,2),1,'-r')
    hold off
    axis image
    title(['Force ' num2str(n_iter)])
    drawnow

    % remove padding
    trac_x = f(2:end-npy+1,2:end-npx+1,1);
    trac_y = f(2:end-npy+1,2:end-npx+1,2);
    trac_x(~real_points) = 0;
    trac_y(~real_points) = 0;
    trac(:,:,1,i_lambda) = trac_x;
    trac(:,:,2,i_lambda) = trac_y;
    % UNITS are traction stress, same as E (pN*um^-2)
    
    soln_norm(i_lambda) = sum(sqrt(trac_x(:).^2 + trac_y(:).^2));

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
    
    disp_e(:,:,:,i_lambda) = u_e(2:end-npy+1,2:end-npx+1,:);

% %     % danuser way
% %     r = sqrt(x.^2 + y.^2);
% %     preFactor = (1 + nu)./(pi*E*r.^3);
% %     G_xx = preFactor.*((1 - nu)*r.^2 + nu*x.^2);
% %     G_xy = preFactor.*nu.*x.*y;
% %     G_yy = preFactor.*((1 - nu)*r.^2 + nu*y.^2);
% %     
% %     G_xx(isnan(G_xx)) = 0;
% %     G_xy(isnan(G_xy)) = 0;
% %     G_yy(isnan(G_yy)) = 0;
% %     
% %     FtG_xx = fft2(G_xx);
% %     FtG_xy = fft2(G_xy);
% %     FtG_yy = fft2(G_yy);
% %     
% %     Ftu_e(:,:,1) = FtG_xx.*Ftf(:,:,1) + FtG_xy.*Ftf(:,:,2);
% %     Ftu_e(:,:,2) = FtG_xy.*Ftf(:,:,1) + FtG_yy.*Ftf(:,:,2);
% % 
% %     u_e(:,:,1) = ifft2(Ftu_e(:,:,1),'symmetric');
% %     u_e(:,:,2) = ifft2(Ftu_e(:,:,2),'symmetric');
% %     
% %     % remove padding
% %     u_e_unpad = u_e(2:end-npy+1,2:end-npx+1,:);
% %     
% %     % something better
% %     dx = dot_spacing;
% %     dy = dot_spacing;
% %     
% %     int_x2_over_r3 = 2*dy*log((dy^2+2*dx*(dx+sqrt(dx^2+dy^2)))/(dy^2));    
% %     int_y2_over_r3 = 2*dx*log((dx^2+2*dy*(dy+sqrt(dx^2+dy^2)))/(dx^2));    
% %     int_1_over_r = int_x2_over_r3 + int_y2_over_r3;
% %         
% %     corrTerm_11 = (1 + nu)/(pi*E)*((1 - nu)*int_1_over_r + nu*int_x2_over_r3);
% %     corrTerm_22 = (1 + nu)/(pi*E)*((1 - nu)*int_1_over_r + nu*int_y2_over_r3);
% %     
% %     u_e_unpad(:,:,1) = u_e_unpad(:,:,1) + trac_x*corrTerm_11;
% %     u_e_unpad(:,:,2) = u_e_unpad(:,:,2) + trac_y*corrTerm_22;
% %     
% % %     x0 = positions(:,:,1);
% % %     y0 = positions(:,:,2);
% % %     Nx_F=size(x0,2); % 2^10 is the densest sampling possible.
% % %     Ny_F=size(y0,1);
% % %     xRange=(max(max(x0))-min(min(x0)));
% % %     yRange=(max(max(y0))-min(min(y0)));
% % %     scalingFactor=(xRange*yRange)/(Nx_F*Ny_F);
% %     scalingFactor = dx*dy;
% %     
% %     disp_e(:,:,:,i_lambda) = scalingFactor*u_e_unpad;
    
% % %     figure(1)
% % %     sf = 5;
% % %     quiver(x,y,u(:,:,1)*sf,u(:,:,2)*sf,0,'-k')
% % %     hold on
% % %     quiver(positions(:,:,1),positions(:,:,2),disp_e(:,:,1,i_lambda)*sf,disp_e(:,:,2,i_lambda)*sf,0,'-r','linewidth',2)
% % %     hold off
% % %     title(sprintf('%i',i_lambda));
% % %     drawnow
% % %     pause(0.1)
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