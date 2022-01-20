function [px,py,real_points] = extend_dots(px,py,img,celldata,meta)

[num_Y,num_X] = size(px);

R = [cos(-celldata.rot_angle) -sin(-celldata.rot_angle); sin(-celldata.rot_angle) cos(-celldata.rot_angle)];

temp_p = (R*([px(:), py(:)] - [celldata.N/2,celldata.M/2])')' + [celldata.N/2,celldata.M/2];
px_rot = temp_p(:,1);
py_rot = temp_p(:,2);
px_rot = reshape(px_rot,num_Y,num_X);
py_rot = reshape(py_rot,num_Y,num_X);

% flip things around if you selected the reverse direction
if mean(px_rot(:,1)) > mean(px_rot(:,end))
    px_rot = px_rot(:,end:-1:1);
    py_rot = py_rot(:,end:-1:1);
end
if mean(py_rot(1,:)) > mean(py_rot(end,:))
    px_rot = px_rot(end:-1:1,:);
    py_rot = py_rot(end:-1:1,:);
end

lim = [1,1; 1,celldata.M; celldata.N,1; celldata.N,celldata.M];
lim_rot = (R*(lim - [celldata.N/2,celldata.M/2])')' + [celldata.N/2,celldata.M/2];
lim_rot_x = [min(lim_rot(:,1)), max(lim_rot(:,1))];
lim_rot_y = [min(lim_rot(:,2)), max(lim_rot(:,2))];

xgrid_rot = mean(px_rot,1);
ygrid_rot = mean(py_rot,2)';

xspacing = abs(mean(diff(xgrid_rot)));
xgride_rot = [fliplr(min(xgrid_rot):-xspacing:lim_rot_x(1)), sort(xgrid_rot(2:end-1)), max(xgrid_rot):xspacing:lim_rot_x(2)];
num_Xe = length(xgride_rot);

yspacing = abs(mean(diff(ygrid_rot)));
ygride_rot = [fliplr(min(ygrid_rot):-yspacing:lim_rot_y(1)), sort(ygrid_rot(2:end-1)), max(ygrid_rot):yspacing:lim_rot_y(2)];
num_Ye = length(ygride_rot);

% [pxe_rot,pye_rot] = meshgrid(xgride_rot,ygride_rot);
% 
% [num_Ye,num_Xe] = size(pxe_rot);
% 
% temp_pe = (R2*([pxe_rot(:), pye_rot(:)] - [vM.N/2,vM.M/2])')' + [vM.N/2,vM.M/2];
% pxe = temp_pe(:,1);
% pye = temp_pe(:,2);
% pxe = reshape(pxe,num_Ye,num_Xe);
% pye = reshape(pye,num_Ye,num_Xe);

for iy = 1:num_Y
    jx = 1:num_X;
    pca_coeff = pca([px_rot(iy,jx)', py_rot(iy,jx)']);
    slope = pca_coeff(2,1)/pca_coeff(1,1);
    intercept = -slope*mean(px_rot(iy,jx)) + mean(py_rot(iy,jx));
    spacing = mean(sqrt(sum(diff([px_rot(iy,jx)', py_rot(iy,jx)'],1).^2,2)));
    xlines(iy,:) = [slope, intercept, spacing];
    xspacing = abs(spacing*cos(atan(slope)));
    yspacing = abs(spacing*sin(atan(slope)));
    
    temp_x = px_rot(iy,jx);
    temp_y = py_rot(iy,jx);
    
    [~,sorti] = sort(temp_x);
    
    if iy == 1
        pxe1_rot(iy,:) = [fliplr(min(temp_x):-xspacing:lim_rot_x(1)), temp_x(sorti(2:end-1)), max(temp_x):xspacing:lim_rot_x(2)];
        L1 = [length(fliplr(min(temp_x):-xspacing:lim_rot_x(1)))-1, length(max(temp_x):xspacing:lim_rot_x(2))-1];
    else
        pxe1_rot(iy,:) = [fliplr(min(temp_x):-xspacing:min(temp_x)-xspacing*L1(1)), temp_x(sorti(2:end-1)), max(temp_x):xspacing:max(temp_x)+xspacing*L1(2)];
    end
    if sum(diff(temp_y)) > 0
        pye1_rot(iy,:) = [fliplr(temp_y(sorti(1)):-yspacing:temp_y(sorti(1))-yspacing*L1(1)), temp_y(sorti(2:end-1)), temp_y(sorti(end)):yspacing:temp_y(sorti(end))+yspacing*L1(2)];
    else
        pye1_rot(iy,:) = [fliplr(temp_y(sorti(1)):yspacing:temp_y(sorti(1))+yspacing*L1(1)), temp_y(sorti(2:end-1)), temp_y(sorti(end)):-yspacing:temp_y(sorti(end))-yspacing*L1(2)];
    end
end

check1 = true;
while check1
    if any(pxe1_rot(:,1) < lim_rot_x(1))
        pxe1_rot(:,1) = [];
        pye1_rot(:,1) = [];
        L1(1) = L1(1) - 1;
        xlines(1,:) = [];
        check1 = true;
    else
        check1 = false;
    end
    if any(pxe1_rot(:,end) > lim_rot_x(2))
        pxe1_rot(:,end) = [];
        pye1_rot(:,end) = [];
        L1(2) = L1(2) - 1;
        xlines(end,:) = [];
        check1 = check1 | true;
    else
        check1 = check1 | false;
    end
end

for jx = 1:num_X
    iy = 1:num_Y;
    pca_coeff = pca([px_rot(iy,jx), py_rot(iy,jx)]);
    slope = pca_coeff(2,1)/pca_coeff(1,1);
    intercept = -slope*mean(px_rot(iy,jx)) + mean(py_rot(iy,jx));
    spacing = mean(sqrt(sum(diff([px_rot(iy,jx), py_rot(iy,jx)],1).^2,2)));
    ylines(jx,:) = [slope, intercept, spacing];
    xspacing = abs(spacing*cos(atan(slope)));
    yspacing = abs(spacing*sin(atan(slope)));
    
    temp_x = px_rot(iy,jx);
    temp_y = py_rot(iy,jx);
    
    [~,sorti] = sort(temp_y);
    
    if jx == 1
        pye2_rot(:,jx) = [fliplr(min(temp_y):-yspacing:lim_rot_y(1))'; temp_y(sorti(2:end-1)); (max(temp_y):yspacing:lim_rot_y(2))'];
        L2 = [length(fliplr(min(temp_y):-yspacing:lim_rot_y(1)))-1; length(max(temp_y):yspacing:lim_rot_y(2))-1];
    else
        pye2_rot(:,jx) = [fliplr(min(temp_y):-yspacing:min(temp_y)-yspacing*L2(1))'; temp_y(sorti(2:end-1)); (max(temp_y):yspacing:max(temp_y)+yspacing*L2(2))'];
    end
    if sum(diff(temp_x)) < 0
        pxe2_rot(:,jx) = [fliplr(temp_x(sorti(1)):xspacing:temp_x(sorti(1))+xspacing*L2(1))'; temp_x(sorti(2:end-1)); (temp_x(sorti(end)):-xspacing:temp_x(sorti(end))-xspacing*L2(2))'];
    else
        pxe2_rot(:,jx) = [fliplr(temp_x(sorti(1)):-xspacing:temp_x(sorti(1))-xspacing*L2(1))'; temp_x(sorti(2:end-1)); (temp_x(sorti(end)):xspacing:temp_x(sorti(end))+xspacing*L2(2))'];
    end
end

check1 = true;
while check1
    if any(pye2_rot(1,:) < lim_rot_y(1))
        pxe2_rot(1,:) = [];
        pye2_rot(1,:) = [];
        L2(1) = L2(1) - 1;
        ylines(1,:) = [];
        check1 = true;
    else
        check1 = false;
    end
    if any(pye2_rot(end,:) > lim_rot_y(2))
        pxe2_rot(end,:) = [];
        pye2_rot(end,:) = [];
        L2(2) = L2(2) - 1;
        ylines(end,:) = [];
        check1 = check1 | true;
    else
        check1 = check1 | false;
    end
end

num_Ye = size(pxe2_rot,1);
num_Xe = size(pxe1_rot,2);

pxe_rot = zeros(num_Ye,num_Xe);
pxe_rot(L2(1)+1:L2(1)+size(pxe1_rot,1),:) = pxe1_rot;
pxe_rot(:,L1(1)+1:L1(1)+size(pxe2_rot,2)) = pxe2_rot;
pye_rot = zeros(num_Ye,num_Xe);
pye_rot(L2(1)+1:L2(1)+size(pxe1_rot,1),:) = pye1_rot;
pye_rot(:,L1(1)+1:L1(1)+size(pxe2_rot,2)) = pye2_rot;

for jx = 1:num_Xe
    for iy = 1:num_Ye
        if pxe_rot(iy,jx) == 0
            xpoints = pxe_rot(iy,L1(1)+1:num_Xe-L1(2));
            ypoints = pye_rot(iy,L1(1)+1:num_Xe-L1(2));
            pca_coeff = pca([xpoints(:), ypoints(:)]);
            slope_x = pca_coeff(2,1)/pca_coeff(1,1);
            intercept_x = -slope_x*mean(xpoints) + mean(ypoints);
            
            xpoints = pxe_rot(L2(1)+1:num_Ye-L2(2),jx);
            ypoints = pye_rot(L2(1)+1:num_Ye-L2(2),jx);
            pca_coeff = pca([xpoints(:), ypoints(:)]);
            slope_y = pca_coeff(2,1)/pca_coeff(1,1);
            intercept_y = -slope_y*mean(xpoints) + mean(ypoints);
            
            pxe_rot(iy,jx) = (intercept_x - intercept_y)/...
                (slope_y - slope_x);
            pye_rot(iy,jx) = slope_x*pxe_rot(iy,jx) + intercept_x;
        end
    end
end

R2 = [cos(celldata.rot_angle) -sin(celldata.rot_angle); sin(celldata.rot_angle) cos(celldata.rot_angle)];

temp_p = (R2*([pxe_rot(:), pye_rot(:)] - [celldata.N/2,celldata.M/2])')' + [celldata.N/2,celldata.M/2];
pxe = temp_p(:,1);
pye = temp_p(:,2);
pxe = reshape(pxe,num_Ye,num_Xe);
pye = reshape(pye,num_Ye,num_Xe);

[px,py,real_points] = find_centroids_new(pxe,pye,img,celldata,meta); % about 3x faster than find_centroids
% [px,py,real_points] = find_centroids(pxe,pye,img,celldata,meta);

end