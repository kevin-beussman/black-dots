function [px,py,real_points,img_bw] = find_centroids(px,py,img,vM)

dotspacing_px = vM.DotSpacing/vM.Calibration;
bpass_noise = 1;
dotsize_px = 2*round((vM.DotSize/vM.Calibration + 1)/2) - 1;

% img2 = imrotate(mat2gray(imcomplement(bpass_kb2(imcomplement(img),bpass_noise,dotsize_px))) + 1,vM.rot_angle*180/pi,'loose') - 1;
img2 = mat2gray(imcomplement(bpass_kb2(imcomplement(img),bpass_noise,dotsize_px)));

[num_Y,num_X] = size(px);

% figure
% imagesc(img2)
% hold on
% p2 = plot(px(:),py(:),'.r');
% p3 = plot(px(1,1),py(1,1),'.-k');

real_points = false(num_Y,num_X);

img_bw = false(size(img));

% x_ind = [ceil((num_X)/2+0.25):num_X, floor(num_X/2):-1:1];
% y_ind = [ceil((num_Y)/2+0.25):num_Y, floor(num_Y/2):-1:1];

x_ind = [1:floor(num_X/2), num_X:-1:ceil((num_X)/2+0.25)];
y_ind = [1:floor(num_Y/2), num_Y:-1:ceil((num_Y)/2+0.25)];

% % c = 0;
% % x_ind = zeros(1,num_X);
% % while c < num_X/2
% %     c = c + 1;
% %     c2 = 2*c - 1;
% %     x_ind(c2) = 1 + (c - 1);
% %     x_ind(c2 + 1) = num_X - (c - 1);
% % end
% % if mod(num_X,2)
% %     x_ind(end) = [];
% % end
% % c = 0;
% % y_ind = zeros(1,num_Y);
% % while c < num_Y/2
% %     c = c + 1;
% %     c2 = 2*c - 1;
% %     y_ind(c2) = 1 + (c - 1);
% %     y_ind(c2 + 1) = num_Y - (c - 1);
% % end
% % if mod(num_Y,2)
% %     y_ind(end) = [];
% % end

pxorig = px;
pyorig = py;

for j = 1:num_X
    jx = x_ind(j);
    for i = 1:num_Y
        iy = y_ind(i);
        
        
        check1 = true;
        area_old = 0;
        area_old2 = 0;
        
%         if (y_ind(i) ~= ceil((num_Y)/2+0.25)) && (y_ind(i) ~= floor(num_Y/2))
%             disp_y = [(px(y_ind(i-1),jx) - pxorig(y_ind(i-1),jx)), (py(y_ind(i-1),jx) - pyorig(y_ind(i-1),jx))];
%         else
%             disp_y = [NaN, NaN];
%         end
%         
%         if (x_ind(j) ~= ceil((num_X)/2+0.25)) && (x_ind(j) ~= floor(num_X/2))
%             disp_x = [(px(iy,x_ind(j-1)) - pxorig(iy,x_ind(j-1))), (py(iy,x_ind(j-1)) - pyorig(iy,x_ind(j-1)))];
%         else
%             disp_x = [NaN, NaN];
%         end

        if (y_ind(i) ~= 1) && (y_ind(i) ~= num_Y)
            disp_y = [(px(y_ind(i-1),jx) - pxorig(y_ind(i-1),jx)), (py(y_ind(i-1),jx) - pyorig(y_ind(i-1),jx))];
        else
            disp_y = [NaN, NaN];
        end
        
        if (x_ind(j) ~= 1) && (x_ind(j) ~= num_X)
            disp_x = [(px(iy,x_ind(j-1)) - pxorig(iy,x_ind(j-1))), (py(iy,x_ind(j-1)) - pyorig(iy,x_ind(j-1)))];
        else
            disp_x = [NaN, NaN];
        end
        
        disp = nanmean([disp_y; disp_x],1);
        if all(~isnan(disp))
            px(iy,jx) = px(iy,jx) + disp(1);
            py(iy,jx) = py(iy,jx) + disp(2);
        end
        
        count = 0;
        while check1
            rect = [px(iy,jx), py(iy,jx), 0, 0] + [-0.5, -0.5, 1, 1]*dotspacing_px;
%             img_crop = imcrop(img2,rect);
            radius = dotspacing_px;
            circ = [px(iy,jx), py(iy,jx)] + radius*[cos(linspace(0,2*pi,10))', sin(linspace(0,2*pi,10))'];
            
            if ~any(any(circ < 1) | any(circ > [vM.N, vM.M])) 
                mask = poly2mask(circ(:,1),circ(:,2),vM.M,vM.N);
                img_crop = imcrop((img2  + 1).*mask - 1,rect);
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
% %                 while ind(1) < 1
% %                     ind(1) = 1;
% %                     ind(2) = ind(2) + 1;
% %                 end
% %                 while ind(3) < 1
% %                     ind(3) = 1;
% %                     ind(4) = ind(4) + 1;
% %                 end
                img_bw(ind(1):ind(2),ind(3):ind(4)) = img_crop_filt_bw_watershed;

%                 temp = regionprops(bwselect(img_crop_filt_bw_watershed,dotspacing_px/2,dotspacing_px/2),'Centroid','Area','Eccentricity','EquivDiameter');
%                 if isempty(temp)
%                     temp = regionprops(img_crop_filt_bw_watershed,'Centroid','Area','Eccentricity','EquivDiameter');
% %                     temp = temp([temp.Area] == max([temp.Area]));
%                     temp = temp([temp.Eccentricity] < 0.6 & [temp.EquivDiameter] > 0.5*dotsize_px);
%                     if length(temp) > 1
%                         dist = sqrt(sum((vertcat(temp.Centroid) - [N/2,M/2]).^2,2));
%                         temp = temp(dist == min(dist));
%                     end
%                 end
                temp = regionprops(L,'Centroid','Area','Eccentricity','EquivDiameter');
%                 temp = temp([temp.Eccentricity] < 0.8 & [temp.EquivDiameter] > 0.75*dotsize_px);
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
                    
%                     if temp.Eccentricity > 0.7
%                         % more processing required to extract dot
%                         temp2 = regionprops(bwulterode(img_final),'Centroid');
%                         if length(temp2) > 1
%                             dist = sqrt(sum((vertcat(temp2.Centroid) - [N/2,M/2]).^2,2));
%                             temp2 = temp2(dist == min(dist));
%                         end
%                         temp.Centroid = temp2.Centroid;
%                     end
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
                
%                 [Gmag,~] = imgradient(mat2gray(img_crop));
%                 meangrad = mean(Gmag(:));
%                 if meangrad > 0.5 % this value is trial-and-error
%                     check1 = false;
%                     real_points(iy,jx) = false;
%                     continue
%                 end
                
%                 if isempty(temp)
% %                     imagesc(img_crop)
% %                     hold on
% %                     plot(temp.Centroid(1),temp.Centroid(2),'.r')
% %                     hold off
%                     check1 = false;
%                     real_points(iy,jx) = false;
%                     continue
%                 end
                
                if (temp.Area == area_old) || (temp.Area == area_old2) % sometimes gets in a loop between 2 values
                    check1 = false;
                end
                
                area_old2 = area_old;
                area_old = temp.Area;

                px(iy,jx) = px(iy,jx) - 0.5*dotspacing_px + temp.Centroid(1) - 1;
                py(iy,jx) = py(iy,jx) - 0.5*dotspacing_px + temp.Centroid(2) - 1;
                
%                 if area_old2 == 0
%                     set(p2,'XData',px(:),'YData',py(:))
%                     set(p3,'XData',[pxold, px(iy,jx)],'YData',[pyold, py(iy,jx)])
%                     pause
%                 end

                real_points(iy,jx) = true;
            else
                % do nothing
                check1 = false;
                real_points(iy,jx) = false;
            end
        end
    end
end

% remove rows if there aren't enough points
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
        
        if isinf(vM.uPoints)
            % do nothing
        elseif nnz(real_points(:,jx)) < 2*vM.uPoints
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
        
        if isinf(vM.uPoints)
            % do nothing
        elseif nnz(real_points(iy,:)) < 2*vM.uPoints
            Lremove = [Lremove; iy];
            check2 = true;
        end
    end
    px(Lremove,:) = [];
    py(Lremove,:) = [];
    real_points(Lremove,:) = [];
    [num_Y,num_X] = size(px);
end