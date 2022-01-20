function [px_full, py_full, real_points] = track_dots_across_frames(images,px,py,celldata,meta)
real_points = celldata.real_points;

px_full = repmat(px,1,1,meta.nFrames);
py_full = repmat(py,1,1,meta.nFrames);

% already analyzed uFrame, analyze outwards from that frame
ind_k = [(meta.uFrame-1:-1:1) (meta.uFrame+1:meta.nFrames)];

% figure
% i1 = imagesc(images(:,:,1));
% hold on
% plot(px,py,'.r')
% p1 = plot(px(1,1),py(1,1),'or');
% hold off

pct = 0;
dnum = 1;
nnum = meta.nFrames;
nlen = length(num2str(nnum));
% fprintf(['[%-20s] %' num2str(nlen) '.0f/%-' num2str(nlen) '.0f\n'],repmat('|',1,round(pct*20)),dnum,nnum)
fprintf('[%-20s] %3.0f%%\n',repmat('|',1,round(pct/100*20)),pct)

for kk = 1:length(ind_k)
    k = ind_k(kk);
    if (k == meta.uFrame-1) || (k == meta.uFrame+1)
        kprev = meta.uFrame;
    else
        kprev = ind_k(kk-1);
    end
    
    [pxk,pyk,real_pointsk] = find_centroids_new(px_full(:,:,kprev),py_full(:,:,kprev),images(:,:,k),celldata,meta);
%     [pxk,pyk] = find_centroids(px_full(:,:,kprev),py_full(:,:,kprev),images(:,:,k),celldata,meta);
    
    if real_pointsk ~= real_points
        % something went wrong probably
        % one of the dots might have moved too close to the image boundary
        % should update real_points to exclude that one
    end
    
    px_full(:,:,k) = pxk;
    py_full(:,:,k) = pyk;

    pct = pct + 1/(nnum-1);
    dnum = dnum + 1;
%     fprintf([repmat('\b',1,25+2*nlen) '[%-20s] %' num2str(nlen) '.0f/%-' num2str(nlen) '.0f\n'],repmat('|',1,round(pct*20)),dnum,nnum)
    fprintf([repmat('\b',1,28) '[%-20s] %3.0f%%\n'],repmat('|',1,round(pct*20)),pct*100)
end
fprintf(repmat('\b',1,25+2*nlen))
