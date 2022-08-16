function [px,py,real_points,img_bw] = find_centroids(px,py,img,celldata,meta)
% px,py are a grid of points where we think a dot should exist
% this code takes those initial guesses and finds the actual centroid of
% each dot
% real_points tells us if those points are actual dots or imaginary ones

dotspacing_px = meta.DotSpacing/meta.Calibration;
% bpass_noise = 1;
dotsize_px = 2*round((meta.DotSize/meta.Calibration + 1)/2) - 1;

% img = mat2gray(imcomplement(bpass_kb2(imcomplement(img),bpass_noise,dotsize_px)));

[num_Y,num_X] = size(px);

% figure
% imagesc(img2)
% hold on
% p2 = plot(px(:),py(:),'.r');
% p3 = plot(px(1,1),py(1,1),'.-k');

real_points = false(num_Y,num_X);

img_bw = false(size(img));

x_ind = [1:floor(num_X/2), num_X:-1:ceil((num_X)/2+0.25)];
y_ind = [1:floor(num_Y/2), num_Y:-1:ceil((num_Y)/2+0.25)];

pxorig = px;
pyorig = py;

% loop over all the points
for j = 1:num_X
    jx = x_ind(j);
    for i = 1:num_Y
        iy = y_ind(i);
        
        check1 = true;
        area_old = 0;
        area_old2 = 0;
        
        % figure out how the previous dot "left of" this one was moved
        disp_y = [];
        if (y_ind(i) ~= 1) && (y_ind(i) ~= num_Y)
            disp_y = [(px(y_ind(i-1),jx) - pxorig(y_ind(i-1),jx)), (py(y_ind(i-1),jx) - pyorig(y_ind(i-1),jx))];
        end
        
        % figure out how the previous dot "above" this one was moved
        disp_x = [];
        if (x_ind(j) ~= 1) && (x_ind(j) ~= num_X)
            disp_x = [(px(iy,x_ind(j-1)) - pxorig(iy,x_ind(j-1))), (py(iy,x_ind(j-1)) - pyorig(iy,x_ind(j-1)))];
        end
        
        if isempty(disp_x) && isempty(disp_y)
            disp = [];
        elseif isempty(disp_x)
            disp = disp_y;
        elseif isempty(disp_y)
            disp = disp_x;
        else
            disp = mean([disp_y; disp_x],1);
        end
        
        % move the current dot according to how the previous left/above ones moved
        if ~isempty(disp)
            px(iy,jx) = px(iy,jx) + disp(1);
            py(iy,jx) = py(iy,jx) + disp(2);
        end
        
        % do the object tracking
        count = 0;
        while check1
            % rectangular bounding area to check
            rect = [px(iy,jx), py(iy,jx), 0, 0] + [-0.5, -0.5, 1, 1]*dotspacing_px;

            % circular bounding area to check
            radius = dotspacing_px;
            circ = [px(iy,jx), py(iy,jx)] + radius*[cos(linspace(0,2*pi,10))', sin(linspace(0,2*pi,10))'];
            
            if ~any(any(circ < 1) | any(circ > [celldata.N, celldata.M])) % check that circle is still in the cell's crop boundary
                mask = poly2mask(circ(:,1),circ(:,2),celldata.M,celldata.N);
                img_crop = imcrop((img  + 1).*mask - 1,rect);
                img_crop(img_crop < 0) = max(img_crop(:));
                [M,N] = size(img_crop);
            
                img_crop_filt = imgaussfilt(img_crop,1);
                img_crop_filt_bw = imfill(imcomplement(imbinarize(img_crop_filt)),'holes');
                D = -bwdist(~img_crop_filt_bw);
                mask = imextendedmin(D,2);
                D2 = imimposemin(D,mask);
                D2(~img_crop_filt_bw) = Inf;
                L = watershed(D2,8);
                L(~img_crop_filt_bw) = 0;
                img_crop_filt_bw_watershed = L > 0;
                
                ind = [round(rect(2)) round(rect(2))+size(img_crop_filt_bw_watershed,1)-1,...
                    round(rect(1)) round(rect(1))+size(img_crop_filt_bw_watershed,2)-1];

                img_bw(ind(1):ind(2),ind(3):ind(4)) = img_crop_filt_bw_watershed;

                temp = regionprops(L,'Centroid','Area','Eccentricity','EquivDiameter');

                L_good = find([temp.EquivDiameter] > 0.25*dotsize_px); % this checks that it is a dot at all
                temp = temp(L_good);
                if length(L_good) > 1
                    dist = sqrt(sum((vertcat(temp.Centroid) - [N/2,M/2]).^2,2));
                    L_good = L_good(dist == min(dist));
                    temp = temp(dist == min(dist));
                end
                if ~isempty(L_good)
                    img_final = bwmorph((L == L_good),'close');
                    temp = regionprops(img_final,'Centroid','Area','Eccentricity','EquivDiameter');
                    
                else
                    check1 = false;
                    real_points(iy,jx) = false;
                    continue
                end
                
                count = count + 1;
                if count > 10 % stuck in a loop, probably not a real dot.
                    check1 = false;
                    real_points(iy,jx) = false;
                    continue
                end
                
                if (temp.Area == area_old) || (temp.Area == area_old2) % sometimes gets in a loop between 2 values
                    check1 = false;
                end
                
                area_old2 = area_old;
                area_old = temp.Area;

                px(iy,jx) = px(iy,jx) - 0.5*dotspacing_px + temp.Centroid(1) - 1;
                py(iy,jx) = py(iy,jx) - 0.5*dotspacing_px + temp.Centroid(2) - 1;

                real_points(iy,jx) = true;
            else
                % dot is too close to image boundary, do not analyze it
                check1 = false;
                real_points(iy,jx) = false;
            end
        end
    end
end

% remove rows if there aren't at least 2*uPoints points
check2 = true;
while check2
    check2 = false;
    Lremove = [];
    for jx = 1:num_X
        if nnz(real_points(:,jx)) < 2
            Lremove = [Lremove; jx];
            check2 = true;
            continue
        end
        
        if isinf(meta.uPoints)
            % do nothing
        elseif nnz(real_points(:,jx)) < 2*meta.uPoints
            Lremove = [Lremove; jx];
            check2 = true;
        end
    end
    px(:,Lremove) = [];
    py(:,Lremove) = [];
    real_points(:,Lremove) = [];
    [num_Y,num_X] = size(px);

    Lremove = [];
    for iy = 1:num_Y
        if nnz(real_points(iy,:)) < 2
            Lremove = [Lremove; iy];
            check2 = true;
            continue
        end
        
        if isinf(meta.uPoints)
            % do nothing
        elseif nnz(real_points(iy,:)) < 2*meta.uPoints
            Lremove = [Lremove; iy];
            check2 = true;
        end
    end
    px(Lremove,:) = [];
    py(Lremove,:) = [];
    real_points(Lremove,:) = [];
    [num_Y,num_X] = size(px);
end