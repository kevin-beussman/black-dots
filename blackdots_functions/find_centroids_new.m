function [px,py,real_points,img_bw] = find_centroids_new(px,py,img,celldata,meta)
% px,py are a grid of points where we think a dot should exist
% this code takes those initial guesses and finds the actual centroid of
% each dot
% real_points tells us if those points are actual dots or imaginary ones

dotspacing_px = meta.DotSpacing/meta.Calibration;
dotsize_px = 2*round((meta.DotSize/meta.Calibration + 1)/2) - 1;

[num_Y,num_X] = size(px);

if ~isfield(celldata,'real_points') || isempty(celldata.real_points)
    real_points = true(num_Y,num_X);
else
    real_points = celldata.real_points;
end

x_ind = [1:floor(num_X/2), num_X:-1:ceil((num_X)/2+0.25)];
y_ind = [1:floor(num_Y/2), num_Y:-1:ceil((num_Y)/2+0.25)];

pxorig = px;
pyorig = py;

img_filt = imgaussfilt(img,1);
img_filt_bw = imfill(imcomplement(imbinarize(img_filt)),'holes');
D = bwdist(~img_filt_bw);
D2 = -bwdist(img_filt_bw);
hills = mat2gray(D+D2);

% loop over all the points
for j = 1:num_X
    jx = x_ind(j);
    for i = 1:num_Y
        iy = y_ind(i);
        if real_points(iy,jx)
            
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
            pxr = round(px(iy,jx));
            pyr = round(py(iy,jx));
            if pxr < 0.5*dotspacing_px || pxr > (celldata.N - 0.5*dotspacing_px) || pyr < 0.5*dotspacing_px || pyr > (celldata.M - 0.5*dotspacing_px)
                real_points(iy,jx) = false;
                continue
            end
            % neighborhood is a 3x3 area of pixels
            neighborhood = hills(pyr-1:pyr+1,pxr-1:pxr+1);
            while neighborhood(2,2) ~= max(neighborhood(:))
                if pxr < 0.5*dotspacing_px || pxr > (celldata.N - 0.5*dotspacing_px) || pyr < 0.5*dotspacing_px || pyr > (celldata.M - 0.5*dotspacing_px)
                    real_points(iy,jx) = false;
                    break
                else
    %                 wx = [1 0 -1;2 0 -2;1 0 -1]; % sobel gradient weights
    %                 wy = [1 2 1;0 0 0;-1 -2 -1];
    %                 Gx = sum(neighborhood(:).*wx(:));
    %                 Gy = sum(neighborhood(:).*wy(:));
    %                 Gmag = sqrt(Gx^2+Gy^2);
    %                 Gdir = -atan2(Gy,-Gx)*180/pi;
    % %                 [Gmag2,Gdir2] = imgradient(neighborhood);
    % %                 if Gdir2(2,2) ~= Gdir
    % %                     1;
    % %                 end
    %                 
    %                 pxr = pxr + cosd(Gdir)/abs(cosd(Gdir))*(abs(cosd(Gdir)) > sqrt(2)/2);
    %                 pyr = pyr + sind(Gdir)/abs(sind(Gdir))*(abs(sind(Gdir)) > sqrt(2)/2);
                    
                    [~,ind_max] = max(neighborhood,[],'all','linear');
                    [y_max,x_max] = ind2sub(size(neighborhood),ind_max);
                    pxr = pxr - 2 + x_max;
                    pyr = pyr - 2 + y_max;
                    neighborhood = hills(pyr-1:pyr+1,pxr-1:pxr+1);
                end
            end
            px(iy,jx) = pxr;
            py(iy,jx) = pyr;
        end
    end
end
% fprintf('Done')

% remove rows/columns if there aren't at least 2*uPoints points
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

% now update the points to the actual centroid (instead of the rounded
% pixel value)
img_bw = bwselect(img_filt_bw,px(real_points),py(real_points));
CC = bwconncomp(img_bw);

% get centroid of each object
CC.Centroids = [];
for i = 1:CC.NumObjects
    [y,x] = ind2sub(size(img_bw),CC.PixelIdxList{i});
    CC.Centroids(i,:) = [mean(x),mean(y)];
end

% find closest centroid to each px,py, set px,py = centroid
[x_ind_real,y_ind_real] = find(real_points);
ind_real = sub2ind(size(real_points),x_ind_real,y_ind_real);
for np = ind_real'
    dist = sqrt(sum((CC.Centroids - [px(np),py(np)]).^2,2));
    [~,idx_min] = min(dist);
    px(np) = CC.Centroids(idx_min,1);
    py(np) = CC.Centroids(idx_min,2);
end

% figure
% imagesc(img_bw)
% hold on
% plot(x_ind_bw,y_ind_bw,'.r')
% hold off

% THIS BELOW IS VERY TIME CONSUMING
% x_ind = [1:floor(num_X/2), num_X:-1:ceil((num_X)/2+0.25)];
% y_ind = [1:floor(num_Y/2), num_Y:-1:ceil((num_Y)/2+0.25)];
% for j = 1:num_X
%     jx = x_ind(j);
%     for i = 1:num_Y
%         iy = y_ind(i);
%         if real_points(iy,jx)
%             
%             [~,idx_bw] = bwselect(img_filt_bw,px(iy,jx),py(iy,jx));
%             [y_bw,x_bw] = ind2sub(size(img),idx_bw);
%     
%             px(iy,jx) = mean(x_bw);
%             py(iy,jx) = mean(y_bw);
%         end
%     end
% end