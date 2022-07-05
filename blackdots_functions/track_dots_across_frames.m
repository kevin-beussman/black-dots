function [px_full, py_full, real_points, px0, py0] = track_dots_across_frames(images,px,py,px0,py0,celldata,meta)
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
    
    [pxk,pyk,real_pointsk,~,removed_rows,removed_cols] = find_centroids_new(px_full(:,:,kprev),py_full(:,:,kprev),images(:,:,k),celldata,meta);
%     [pxk,pyk] = find_centroids(px_full(:,:,kprev),py_full(:,:,kprev),images(:,:,k),celldata,meta);
    
    if ~isequal(real_pointsk,real_points)
        % something went wrong probably
        % one of the dots might have moved too close to the image boundary
        % should update real_points to exclude that one
        % it's possible this would reduce #points in a row to less than 4,
        % removing the row
        
%         [Mrk,Nrk] = size(real_pointsk);
%         [Mr,Nr] = size(real_points);
%         
%         removed_rows = [];
%         removed_cols = [];
%         if Mrk ~= Mr && Nrk == Nr
%             ix = 1;
%             while ix <= Mrk && isequal(real_pointsk(ix,:),real_points(ix,:))
%                 ix = ix + 1;
%             end
%             removed_rows = ix;
%         elseif Nrk ~= Nr && Mrk == Mr
%             jy = 1;
%             while jy <= Nrk && isequal(real_pointsk(:,jy),real_points(:,jy))
%                 jy = jy + 1;
%             end
%             removed_cols = jy;
%         elseif Mrk ~= Mr && Nrk ~= Nr
%             for ix = 1:Mr
%                 for jy = 1:Nr
%                     if isequal(real_points([1:ix-1, ix+1:Mr],[1:jy-1,jy+1:Nr]),real_pointsk)
%                         removed_rows = ix;
%                         removed_cols = jy;
%                     end
%                 end
%             end
%         else
%             [update_row,update_col] = find(real_pointsk ~= real_points);
%             pxk(update_row,update_col) = px_full(update_row,update_col,meta.uFrame);
%             pyk(update_row,update_col) = py_full(update_row,update_col,meta.uFrame);
%             for up = 1:length(update_row)
%                 px_full(update_row(up),update_col(up),:) = px_full(update_row(up),update_col(up),meta.uFrame);
%                 py_full(update_row(up),update_col(up),:) = py_full(update_row(up),update_col(up),meta.uFrame);
%             end
%         end

        if ~isempty(removed_rows)
            px_full = px_full([1:removed_rows-1,removed_rows+1:end],:,:);
            py_full = py_full([1:removed_rows-1,removed_rows+1:end],:,:);
            px0 = px0([1:removed_rows-1,removed_rows+1:end],:);
            py0 = py0([1:removed_rows-1,removed_rows+1:end],:);
            real_points = real_points([1:removed_rows-1,removed_rows+1:end],:);
        end
        if ~isempty(removed_cols)
            px_full = px_full(:,[1:removed_cols-1,removed_cols+1:end],:);
            py_full = py_full(:,[1:removed_cols-1,removed_cols+1:end],:);
            px0 = px0(:,[1:removed_cols-1,removed_cols+1:end]);
            py0 = py0(:,[1:removed_cols-1,removed_cols+1:end]);
            real_points = real_points(:,[1:removed_cols-1,removed_cols+1:end]);
        end

        [update_row,update_col] = find(real_pointsk ~= real_points);
        if ~isempty(update_row)
            pxk(update_row,update_col) = px_full(update_row,update_col,meta.uFrame);
            pyk(update_row,update_col) = py_full(update_row,update_col,meta.uFrame);
            for up = 1:length(update_row)
                px_full(update_row(up),update_col(up),:) = px_full(update_row(up),update_col(up),meta.uFrame);
                py_full(update_row(up),update_col(up),:) = py_full(update_row(up),update_col(up),meta.uFrame);
            end
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
