function [px_full, py_full, real_points] = track_dots_across_frames(images,px,py,real_points,celldata,meta)

dotspacing_px = meta.DotSpacing/meta.Calibration;
dotsize_px = 2*round((meta.DotSize/meta.Calibration + 1)/2) - 1;

[num_Y,num_X] = size(px);

px_full = repmat(px,1,1,meta.nFrames);
py_full = repmat(py,1,1,meta.nFrames);

ind_k = [(meta.uFrame-1:-1:1) (meta.uFrame+1:meta.nFrames)];

% figure
% i1 = imagesc(images(:,:,1));
% hold on
% plot(px,py,'.r')
% p1 = plot(px(1,1),py(1,1),'or');
% hold off

pct = 0;
dnum = 0;
nnum = meta.nFrames;
nlen = length(num2str(nnum));
fprintf(['[%-20s] %' num2str(nlen) '.0f/%-' num2str(nlen) '.0f\n'],repmat('|',1,round(pct*20)),dnum,nnum)

for kk = 1:length(ind_k)
    k = ind_k(kk);
    if (k == meta.uFrame-1) || (k == meta.uFrame+1)
        kprev = meta.uFrame;
    else
        kprev = ind_k(kk-1);
    end
    
    [pxk,pyk] = find_centroids_new(px_full(:,:,kprev),py_full(:,:,kprev),images(:,:,k),celldata,meta);
%     [pxk,pyk] = find_centroids(px_full(:,:,kprev),py_full(:,:,kprev),images(:,:,k),celldata,meta);
    
    px_full(:,:,k) = pxk;
    py_full(:,:,k) = pyk;

    pct = pct + 1/nnum;
    dnum = dnum + 1;
    fprintf([repmat('\b',1,25+2*nlen) '[%-20s] %' num2str(nlen) '.0f/%-' num2str(nlen) '.0f\n'],repmat('|',1,round(pct*20)),dnum,nnum)
end

% for jx = 1:num_X
%     for iy = 1:num_Y
%         if real_points(iy,jx)
% %             set(p1,'XData',px(iy,jx),'YData',py(iy,jx))
%             bad_dot = false;
%             for kk = 1:length(ind_k)
%                 k = ind_k(kk);
%                 img2 = images(:,:,k);
%                 
%                 if k < meta.uFrame
%                     rect = [px_full(iy,jx,k+1), py_full(iy,jx,k+1), 0, 0] + [-0.5, -0.5, 1, 1]*dotspacing_px;
%                     rect = round(rect);
%                 else
%                     rect = [px_full(iy,jx,k-1), py_full(iy,jx,k-1), 0, 0] + [-0.5, -0.5, 1, 1]*dotspacing_px;
%                     rect = round(rect);
%                 end
%                 
%                 radius = dotspacing_px/2;
%                 circ = [px_full(iy,jx,k-1), py_full(iy,jx,k-1)] + radius*[cos(linspace(0,2*pi,10))', sin(linspace(0,2*pi,10))'];
%                 mask = poly2mask(circ(:,1),circ(:,2),celldata.M,celldata.N);
%                 img_crop = imcrop((img2  + 1).*mask - 1,rect);
%                 img_crop(img_crop < 0) = max(img_crop(:));
%                 [M,N] = size(img_crop);
%                 
%                 img_crop_filt = imgaussfilt(img_crop,1);
%                 img_crop_filt_bw = ~imbinarize(img_crop_filt);
%                 D = -bwdist(~img_crop_filt_bw);
%                 mask = imextendedmin(D,2);
%                 D2 = imimposemin(D,mask);
%                 D2(~img_crop_filt_bw) = Inf;
%                 L = watershed(D2,8);
%                 L(~img_crop_filt_bw) = 0;
%                 img_crop_filt_bw_watershed = L > 0;
%                 
%                 img_final = bwselect(img_crop_filt_bw_watershed,N/2,M/2);
%                 temp = regionprops(img_final,'Centroid','Area','Eccentricity');
%                 if isempty(temp)
%                     temp = regionprops(img_crop_filt_bw_watershed,'Centroid','Area','Eccentricity');
% %                     temp = temp([temp.Area] == max([temp.Area]));
%                     if length(temp) > 1
%                         dist = sqrt(sum((vertcat(temp.Centroid) - [N/2,M/2]).^2,2));
%                         temp = temp(dist == min(dist));
%                     end
%                 end
%                 
%                 if isempty(temp)
%                     % dot is bad, do something!
%                     bad_dot = true;
%                     break
%                 end
%                 
%                 if temp.Eccentricity > 0.7
%                     % more processing required to extract dot
%                     temp2 = regionprops(bwulterode(img_final),'Centroid');
%                     if length(temp2) > 1
%                         dist = sqrt(sum((vertcat(temp2.Centroid) - [N/2,M/2]).^2,2));
%                         temp2 = temp2(dist == min(dist));
%                     end
%                     temp.Centroid = temp2.Centroid;
%                 end
%                 
%                 px_full(iy,jx,k) = px_full(iy,jx,k-1) - 0.5*N + temp.Centroid(1) - 1;
%                 py_full(iy,jx,k) = py_full(iy,jx,k-1) - 0.5*M + temp.Centroid(2) - 1;
%                 
% %                 set(i1,'CData',img2)
% %                 set(p1,'XData',px_full(iy,jx,k),'YData',py_full(iy,jx,k))
%             end
%             if bad_dot
%                 real_points(iy,jx) = false;
%                 px_full(iy,jx,:) = px_full(iy,jx,meta.uFrame);
%                 py_full(iy,jx,:) = py_full(iy,jx,meta.uFrame);
%             end
%             
%             
%         end
%     end
% end