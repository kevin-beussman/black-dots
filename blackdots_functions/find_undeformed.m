function [px,py,px0,py0] = find_undeformed(px,py,celldata,meta)
real_points = celldata.real_points;

[num_Y,num_X] = size(px);

%% find slope of each line of dots
xlines = zeros(num_Y,2);
for iy = 1:num_Y
    num_points = sum(real_points(iy,:));
    if num_points >= 2
        if isempty(meta.uPoints)
            jx = find(real_points(iy,:));
        elseif num_points >= 2*meta.uPoints
            jx = find(real_points(iy,:));
            jx = jx([1:meta.uPoints, (end-meta.uPoints+1):end]);
        else
            jx = find(real_points(iy,:));
        end
        pca_coeff = pca([px(iy,jx)', py(iy,jx)']);
        slope = pca_coeff(2,1)/pca_coeff(1,1);
        intercept = -slope*mean(px(iy,jx)) + mean(py(iy,jx));
        xlines(iy,:) = [slope, intercept];
    end
end

ylines = zeros(num_X,2);
for jx = 1:num_X
    num_points = sum(real_points(:,jx));
    if num_points >= 2
        if isempty(meta.uPoints)
            iy = find(real_points(:,jx))';
        elseif num_points >= 2*meta.uPoints
            iy = find(real_points(:,jx))';
            iy = iy([1:meta.uPoints, (end-meta.uPoints+1):end]);
        else
            iy = find(real_points(:,jx))';
        end
        pca_coeff = pca([px(iy,jx), py(iy,jx)]);
        slope = pca_coeff(2,1)/pca_coeff(1,1);
        intercept = -slope*mean(px(iy,jx)) + mean(py(iy,jx));
        ylines(jx,:) = [slope, intercept];
    end
end

%% correct them with median slope
for iy = 1:num_Y
    num_points = sum(real_points(iy,:));
    if num_points >= 2
        if isempty(meta.uPoints)
            jx = find(real_points(iy,:));
        elseif num_points >= 2*meta.uPoints
            jx = find(real_points(iy,:));
            jx = jx([1:meta.uPoints, (end-meta.uPoints+1):end]);
        else
            jx = find(real_points(iy,:));
        end
%         pca_coeff = pca([px(iy,jx)', py(iy,jx)']);
        slope = median(xlines(:,1));
        intercept = -slope*mean(px(iy,jx)) + mean(py(iy,jx));
        xlines(iy,:) = [slope, intercept];
    end
end

for jx = 1:num_X
    num_points = sum(real_points(:,jx));
    if num_points >= 2
        if isempty(meta.uPoints)
            iy = find(real_points(:,jx))';
        elseif num_points >= 2*meta.uPoints
            iy = find(real_points(:,jx))';
            iy = iy([1:meta.uPoints, (end-meta.uPoints+1):end]);
        else
            iy = find(real_points(:,jx))';
        end
%         pca_coeff = pca([px(iy,jx), py(iy,jx)]);
        slope = median(ylines(:,1));
        intercept = -slope*mean(px(iy,jx)) + mean(py(iy,jx));
        ylines(jx,:) = [slope, intercept];
    end
end

%% find undeformed by intersection of 2 lines

px0 = zeros(size(px));
py0 = zeros(size(py));
for jx = 1:num_X
    for iy = 1:num_Y
        px0(iy,jx) = (xlines(iy,2) - ylines(jx,2))/...
            (ylines(jx,1) - xlines(iy,1));
        py0(iy,jx) = xlines(iy,1)*px0(iy,jx) + xlines(iy,2);
        if ~real_points(iy,jx)
            px(iy,jx) = px0(iy,jx);
            py(iy,jx) = py0(iy,jx);
        end
    end
end
