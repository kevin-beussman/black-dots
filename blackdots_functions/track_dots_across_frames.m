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
    
    if ~isequal(real_pointsk,real_points)
        % something went wrong probably
        % one of the dots might have moved too close to the image boundary
        % should update real_points to exclude that one
        % it's possible this would reduce #points in a row to less than 4,
        % removing the row

        [Mrk,Nrk] = size(real_pointsk);
        [Mr,Nr] = size(real_points);
        
        remove_row = [];
        remove_col = [];
        if Mrk ~= Mr && Nrk == Nr
            ix = 1;
            while ix <= Mrk && isequal(real_pointsk(ix,:),real_points(ix,:))
                ix = ix + 1;
            end
            remove_row = ix;
        elseif Nrk ~= Nr && Mrk == Mr
            jy = 1;
            while jy <= Nrk && isequal(real_pointsk(:,jy),real_points(:,jy))
                jy = jy + 1;
            end
            remove_col = jy;
        elseif Mrk ~= Mr && Nrk ~= Nr
            for ix = 1:Mr
                for jy = 1:Nr
                    if isequal(real_points([1:ix-1, ix+1:Mr],[1:jy-1,jy+1:Nr]),real_pointsk)
                        remove_row = ix;
                        remove_col = jy;
                    end
                end
            end
        else
            [update_row,update_col] = find(real_pointsk ~= real_points);
            pxk(update_row,update_col) = px_full(update_row,update_col,meta.uFrame);
            pyk(update_row,update_col) = py_full(update_row,update_col,meta.uFrame);
            px_full(update_row,update_col,:) = px_full(update_row,update_col,meta.uFrame);
            py_full(update_row,update_col,:) = py_full(update_row,update_col,meta.uFrame);
        end
        
        if ~isempty(remove_row)
            px_full = px_full([1:remove_row-1,remove_row+1:end],:,:);
            py_full = py_full([1:remove_row-1,remove_row+1:end],:,:);
        end
        if ~isempty(remove_col)
            px_full = px_full(:,[1:remove_col-1,remove_col+1:end],:);
            py_full = py_full(:,[1:remove_col-1,remove_col+1:end],:);
        end

        real_points = real_pointsk;
        celldata.real_points = real_points;
    end
    
    px_full(:,:,k) = pxk;
    py_full(:,:,k) = pyk;

    pct = pct + 1/(nnum-1);
    dnum = dnum + 1;
%     fprintf([repmat('\b',1,25+2*nlen) '[%-20s] %' num2str(nlen) '.0f/%-' num2str(nlen) '.0f\n'],repmat('|',1,round(pct*20)),dnum,nnum)
    fprintf([repmat('\b',1,28) '[%-20s] %3.0f%%\n'],repmat('|',1,round(pct*20)),pct*100)
end
fprintf(repmat('\b',1,25+2*nlen))
