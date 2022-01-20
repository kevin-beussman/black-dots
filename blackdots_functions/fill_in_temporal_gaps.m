function [px,py] = fill_in_temporal_gaps(res,celldata,meta)
    % "full points" are points that appear on every frame
    % "close points" are points that appear on at least 75% of frames
    % for each missing frame, fill it in with the nearest frame
    
%     full_points = NaN(max(res(:,5)),1);
    close_points = NaN(max(res(:,5)),1);
    p1x = NaN(max(res(:,5)),meta.nFrames);
    p1y = NaN(max(res(:,5)),meta.nFrames);
    p2x = NaN(max(res(:,5)),meta.nFrames);
    p2y = NaN(max(res(:,5)),meta.nFrames);

    ifp = 0; icp = 0;
    pct = 0;
    point_markers = [0; find(diff(res(:,5))); size(res,1)];
    max_npt = max(res(:,5));
%     fprintf('\n[%-20s] %3.0f%%\n',repmat('|',1,round(pct*20)),pct*100)
    for npt = 1:max_npt
%         pct = npt/max_npt;
%         fprintf([repmat('\b',1,28) '[%-20s] %3.0f%%\n'],repmat('|',1,round(pct*20)),pct*100)      
        i_1 = point_markers(npt)+1; % this is much faster!
        i_2 = point_markers(npt+1);
        if nnz(res(i_1:i_2,5) == npt) == meta.nFrames
            ifp = ifp + 1;
%             full_points(ifp) =  npt;
            p1x(ifp,:) = res(i_1:i_2,1)';
            p1y(ifp,:) = res(i_1:i_2,2)';
        elseif nnz(res(i_1:i_2,5) == npt) >= round(0.75*meta.nFrames)
            icp = icp + 1;
            close_points(icp) = npt;

            inds = find((res(:,5) == close_points(icp)));
            inds2 = res(inds,4);

            tempx = NaN(1,meta.nFrames);
            tempx(inds2) = res(inds,1)';
            tempy = NaN(1,meta.nFrames);
            tempy(inds2) = res(inds,2)';
            p2x(icp,:) = tempx;
            p2y(icp,:) = tempy;
        end
    end
%     fprintf(repmat('\b',1,28))

%     full_points = full_points(~isnan(full_points));
    close_points = close_points(~isnan(close_points));
    p1x = p1x(1:ifp,:);
    p1y = p1y(1:ifp,:);
    p2x = p2x(1:icp,:);
    p2y = p2y(1:icp,:);
    
    p2x_f = p2x;
    p2y_f = p2y;
    for icp = 1:length(close_points)
        for k = find(isnan(p2x(icp,:)))
            cL = (k - 1)*(k > 1) + 1*(k == 1);
            cR = (k + 1)*(k < meta.nFrames) + meta.nFrames*(k == meta.nFrames);
            while isnan(p2x(icp,cL)) && isnan(p2x(icp,cR))
                cL = (cL - 1)*(cL > 1) + 1*(cL == 1);
                cR = (cR + 1)*(cR < meta.nFrames) + meta.nFrames*(cR == meta.nFrames);
            end
            if ~isnan(p2x(icp,cL))
                p2x_f(icp,k) = p2x(icp,cL);
                p2y_f(icp,k) = p2y(icp,cL);
            elseif ~isnan(p2x(icp,cR))
                p2x_f(icp,k) = p2x(icp,cR);
                p2y_f(icp,k) = p2y(icp,cR);
            end
        end
    end

    px = [p1x; p2x_f]; % these are the final x,y positions
    py = [p1y; p2y_f];