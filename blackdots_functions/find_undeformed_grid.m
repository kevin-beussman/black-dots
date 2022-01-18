function [px_grid,py_grid,px0_grid,py0_grid,real_points] = find_undeformed_grid(px,py,vM)
% 1) fit a line through grid points using least-squares to find the
% intercept
% 2) repeat this except increase the angle by 90 degrees
% 3) for each point, find the intersection of the two lines that pass
% through it

% px            x-positions
% py            y-positions
% vM            video metadata
% vM.uPoints    # of points on each end to use for least-squares (inf=all)
    
px_uFrame = px(:,vM.uFrame);
py_uFrame = py(:,vM.uFrame);

%% start finding lines of points
px_temp = px_uFrame;
py_temp = py_uFrame;
check_xlines = true;
num_x_lines = 0;
xline_points = {};
max_skip = 3;
while check_xlines
    p_dist_from_center = sqrt((px_temp - vM.N/2).^2 + (py_temp - vM.M/2).^2);
    np = find(p_dist_from_center == min(p_dist_from_center));
    np_orig = np;
        
    np_line = np;
    slope = tan(vM.rot_angle);
    slope_sign = sign(slope);
    check_dir1 = 1;
    while check_dir1
        for dot_skip = 1:max_skip
            for fp = fliplr(1:-0.1:0.51)
                px_guess = (px_temp(np) + fp*dot_skip*vM.DotSpacing/vM.Calibration*cos(atan(slope)));
                py_guess = (py_temp(np) + fp*dot_skip*vM.DotSpacing/vM.Calibration*sin(atan(slope)));
                p_dist_from_guess = sqrt((px_temp - px_guess).^2 + (py_temp - py_guess).^2);
                
                np_next = find(p_dist_from_guess == min(p_dist_from_guess));
                if any(p_dist_from_guess*vM.Calibration < vM.DotSpacing/2)
                    break
                end
            end
%             px_guess = (px_temp(np) + dot_skip*vM.DotSpacing/vM.Calibration*cos(atan(slope)));
%             py_guess = (py_temp(np) + dot_skip*vM.DotSpacing/vM.Calibration*sin(atan(slope)));
%             p_dist_from_guess = sqrt((px_temp - px_guess).^2 + (py_temp - py_guess).^2);
% 
%             np_next = find(p_dist_from_guess == min(p_dist_from_guess));
% %             figure(1)
% %             plot(px,py,'.k',px(np),py(np),'or',px_guess,py_guess,'*g',px(np_next),py(np_next),'og')
% %             set(gca,'ydir','reverse')
% %             drawnow
            if p_dist_from_guess(np_next)*vM.Calibration < vM.DotSpacing/2
                np_line = [np_line; np_next];
% % %                 slope = (py_temp(np_next) - py_temp(np))/(px_temp(np_next) - px_temp(np));
% %                 slope = slope_sign*abs((py_temp(np_next) - py_temp(np_orig))/(px_temp(np_next) - px_temp(np_orig)));
%                 slope = mean([slope_sign*abs((py_temp(np_next) - py_temp(np_orig))/(px_temp(np_next) - px_temp(np_orig))), slope]);
                np = np_next;
                break
            else
%                 figure(1)
%                 plot(px,py,'.k',px(np),py(np),'or',px_guess,py_guess,'*g',px(np_next),py(np_next),'og')
%                 drawnow

                if dot_skip == max_skip
                    check_dir1 = 0;
                end
            end
        end
    end
    
    np = np_orig;
    slope = tan(vM.rot_angle);
    slope_sign = sign(slope);
    check_dir2 = 1;
    while check_dir2
        for dot_skip = 1:max_skip
            for fp = fliplr(1:-0.1:0.51)
                px_guess = (px_temp(np) - fp*dot_skip*vM.DotSpacing/vM.Calibration*cos(atan(slope)));
                py_guess = (py_temp(np) - fp*dot_skip*vM.DotSpacing/vM.Calibration*sin(atan(slope)));
                p_dist_from_guess = sqrt((px_temp - px_guess).^2 + (py_temp - py_guess).^2);
                
                np_next = find(p_dist_from_guess == min(p_dist_from_guess));
                
                if any(p_dist_from_guess*vM.Calibration < vM.DotSpacing/2)
                    break
                end
            end
%             figure(1)
%             plot(px,py,'.k',px(np),py(np),'or',px_guess,py_guess,'*g',px(np_next),py(np_next),'og')
%             set(gca,'ydir','reverse')
%             drawnow
            if p_dist_from_guess(np_next)*vM.Calibration < vM.DotSpacing/2
                np_line = [np_next; np_line];
%                 slope = mean([slope_sign*abs((py_temp(np_next) - py_temp(np_orig))/(px_temp(np_next) - px_temp(np_orig))), slope]);
                np = np_next;
                break
            else
%                 figure(1)
%                 plot(px,py,'.k',px(np),py(np),'or',px_guess,py_guess,'*g',px(np_next),py(np_next),'og')
%                 drawnow

                if dot_skip == max_skip
                    check_dir2 = 0;
                end
            end
        end
    end
    
    num_x_lines = num_x_lines + 1;
    xline_points{num_x_lines} = np_line;
    
    if ismember(2020,np_line)
        1;
    end

    px_temp(np_line) = NaN;
    py_temp(np_line) = NaN;
    
    if all(isnan(py_temp))
        check_xlines = false;
    end
end

px_temp = px_uFrame;
py_temp = py_uFrame;
check_ylines = true;
num_y_lines = 0;
yline_points = {};
while check_ylines
    p_dist_from_center = sqrt((px_temp - vM.N/2).^2 + (py_temp - vM.M/2).^2);
    np = find(p_dist_from_center == min(p_dist_from_center));
    np_orig = np;
        
    np_line = np;
    slope = tan(vM.rot_angle - pi/2);
    slope_sign = sign(slope);
    check_dir1 = 1;
    while check_dir1
        for dot_skip = 1:max_skip
            for fp = fliplr(1:-0.1:0.51)
                px_guess = (px_temp(np) + fp*dot_skip*vM.DotSpacing/vM.Calibration*cos(atan(slope)));
                py_guess = (py_temp(np) + fp*dot_skip*vM.DotSpacing/vM.Calibration*sin(atan(slope)));
                p_dist_from_guess = sqrt((px_temp - px_guess).^2 + (py_temp - py_guess).^2);
                
                np_next = find(p_dist_from_guess == min(p_dist_from_guess));
                
                if any(p_dist_from_guess*vM.Calibration < vM.DotSpacing/2)
                    break
                end
            end
            if p_dist_from_guess(np_next)*vM.Calibration < vM.DotSpacing/2
                np_line = [np_line; np_next];
%                 slope = mean([slope_sign*abs((py_temp(np_next) - py_temp(np_orig))/(px_temp(np_next) - px_temp(np_orig))), slope]);
                np = np_next;
                break
            else
%                 figure(1)
%                 plot(px,py,'.k',px(np),py(np),'or',px_guess,py_guess,'*g',px(np_next),py(np_next),'og')
%                 drawnow
                
                if dot_skip == max_skip
                    check_dir1 = 0;
                end
            end
        end
    end
    
    np = np_orig;
    slope = tan(vM.rot_angle - pi/2);
    slope_sign = sign(slope);
    check_dir2 = 1;
    while check_dir2
        for dot_skip = 1:max_skip
            for fp = fliplr(1:-0.1:0.51)
                px_guess = (px_temp(np) - fp*dot_skip*vM.DotSpacing/vM.Calibration*cos(atan(slope)));
                py_guess = (py_temp(np) - fp*dot_skip*vM.DotSpacing/vM.Calibration*sin(atan(slope)));
                p_dist_from_guess = sqrt((px_temp - px_guess).^2 + (py_temp - py_guess).^2);
                
                np_next = find(p_dist_from_guess == min(p_dist_from_guess));
                
                if any(p_dist_from_guess*vM.Calibration < vM.DotSpacing/2)
                    break
                end
            end
            if p_dist_from_guess(np_next)*vM.Calibration < vM.DotSpacing/2
                
                if ismember(np_next,np_line)
                    1;
                end
                np_line = [np_next; np_line];
%                 slope = mean([slope_sign*abs((py_temp(np_next) - py_temp(np_orig))/(px_temp(np_next) - px_temp(np_orig))), slope]);
                np = np_next;
                break
            else
%                 figure(1)
%                 plot(px,py,'.k',px(np),py(np),'or',px_guess,py_guess,'*g',px(np_next),py(np_next),'og')
%                 drawnow
                
                if dot_skip == max_skip
                    check_dir2 = 0;
                end
            end
        end
    end
    

    num_y_lines = num_y_lines + 1;
    yline_points{num_y_lines} = np_line;

    px_temp(np_line) = NaN;
    py_temp(np_line) = NaN;
    
    if all(isnan(py_temp))
        check_ylines = false;
    end
end

check_short_lines = 1;
while check_short_lines
    p_to_remove2 = [];
    for n_xline = 1:num_x_lines
        if length(xline_points{n_xline}) <= 3
            p_to_remove2 = [p_to_remove2; xline_points{n_xline}];
        end
%         figure(1)
%         plot(px,py,'.k',px(xline_points{n_xline}),py(xline_points{n_xline}),'or')
%         drawnow
%         pause
    end
    for n_yline = 1:num_y_lines
        if length(yline_points{n_yline}) <= 3
            p_to_remove2 = [p_to_remove2; yline_points{n_yline}];
        end
%         figure(1)
%         plot(px,py,'.k',px(yline_points{n_yline}),py(yline_points{n_yline}),'or')
%         drawnow
%         pause
    end
    
    if isempty(p_to_remove2)
        check_short_lines = 0;
        continue
    end
    
    for np_remove = flipud(unique(p_to_remove2))'
        for n_xline = 1:num_x_lines
            xline_points{n_xline}(xline_points{n_xline} == np_remove) = [];
            xline_points{n_xline}(xline_points{n_xline} > np_remove) = xline_points{n_xline}(xline_points{n_xline} > np_remove) - 1;
        end
        for n_yline = 1:num_y_lines
            yline_points{n_yline}(yline_points{n_yline} == np_remove) = [];
            yline_points{n_yline}(yline_points{n_yline} > np_remove) = yline_points{n_yline}(yline_points{n_yline} > np_remove) - 1;
        end
        xline_points = xline_points(~cellfun(@isempty,xline_points));
        yline_points = yline_points(~cellfun(@isempty,yline_points));
        
        num_x_lines = length(xline_points);
        num_y_lines = length(yline_points);
        
        px(np_remove,:) = [];
        py(np_remove,:) = [];
    end
end

px_uFrame = px(:,vM.uFrame);
py_uFrame = py(:,vM.uFrame);

%% get undeformed point locations
px_temp2 = mean(px,2);
py_temp2 = mean(py,2);
point_line_data = NaN(length(px_uFrame),4);
xline_data = NaN(num_x_lines,2);
yline_data = NaN(num_y_lines,2);
px0 = NaN(length(px_uFrame),1);
py0 = NaN(length(px_uFrame),1);
for np = 1:length(px_uFrame)
    n_xline = find(cellfun(@ismember,num2cell(repelem(np,size(xline_points,1),size(xline_points,2))),xline_points));
    n_yline = find(cellfun(@ismember,num2cell(repelem(np,size(yline_points,1),size(yline_points,2))),yline_points));
    
    check_xlines = 0;
    if size(xline_points{n_xline},1) > 3
        if isnan(xline_data(n_xline,1))
            if vM.uPoints > length(xline_points{n_xline})
                points_to_use = 1:length(xline_points{n_xline});
            else
                points_to_use = [1:vM.uPoints, length(xline_points{n_xline})+(-vM.uPoints:0)];
            end
%             px_xline = px_uFrame(xline_points{n_xline}(points_to_use));
%             py_xline = py_uFrame(xline_points{n_xline}(points_to_use));
            px_xline = px_temp2(xline_points{n_xline}(points_to_use));
            py_xline = py_temp2(xline_points{n_xline}(points_to_use));

            pca_coeff = pca([px_xline, py_xline]);
            slope_xline = pca_coeff(2,1)/pca_coeff(1,1);
            intercept_xline = -slope_xline*mean(px_xline) + mean(py_xline);
            
            point_line_data(np,1:2) = [slope_xline, intercept_xline];
            xline_data(n_xline,:) = [slope_xline, intercept_xline];
        else
            slope_xline = xline_data(n_xline,1);
            intercept_xline = xline_data(n_xline,2);
        end
        check_xlines = 1;
    end
    
    check_ylines = 0;
    if size(yline_points{n_yline},1) > 3
        if isnan(yline_data(n_yline,1))
            if vM.uPoints > length(yline_points{n_yline})
                points_to_use = 1:length(yline_points{n_yline});
            else
                points_to_use = [1:vM.uPoints, length(yline_points{n_yline})+(-vM.uPoints:0)];
            end
%             px_yline = px_uFrame(yline_points{n_yline}(points_to_use));
%             py_yline = py_uFrame(yline_points{n_yline}(points_to_use));
            px_yline = px_temp2(yline_points{n_yline}(points_to_use));
            py_yline = py_temp2(yline_points{n_yline}(points_to_use));

            pca_coeff = pca([px_yline, py_yline]);
            slope_yline = pca_coeff(2,1)/pca_coeff(1,1);
            intercept_yline = -slope_yline*mean(px_yline) + mean(py_yline);
            
            point_line_data(np,3:4) = [slope_yline, intercept_yline];
            yline_data(n_yline,:) = [slope_yline, intercept_yline];
        else
            slope_yline = yline_data(n_yline,1);
            intercept_yline = yline_data(n_yline,2);
        end
        check_ylines = 1;
    end
    
    if check_xlines && check_ylines
        px0(np) = (intercept_xline - intercept_yline)/...
            (slope_yline - slope_xline);
        py0(np) = slope_xline*px0(np) + intercept_xline;
    end
end

for np = find(isnan(px0))'
    n_xline = find(cellfun(@ismember,num2cell(repelem(np,size(xline_points,1),size(xline_points,2))),xline_points));
    n_yline = find(cellfun(@ismember,num2cell(repelem(np,size(yline_points,1),size(yline_points,2))),yline_points));
    
    p_dist_from_np = sqrt((px_uFrame - px_uFrame(np)).^2 + (py_uFrame - py_uFrame(np)).^2);
    p_dist_from_np(np) = NaN;
    
    p_closest = find(p_dist_from_np == min(p_dist_from_np(~isnan(point_line_data(:,1)) & ~isnan(point_line_data(:,3)))));
    
    check_xlines = 0;
    if size(xline_points{n_xline},1) <= 3
        px_xline = px_uFrame(xline_points{n_xline});
        py_xline = py_uFrame(xline_points{n_xline});
        
        slope_xline = point_line_data(p_closest,1);
        intercept_xline = -slope_xline*mean(px_xline) + mean(py_xline);
        
        point_line_data(np,1:2) = [slope_xline, intercept_xline];
        xline_data(n_xline,1:2) = [slope_xline, intercept_xline];
    
        check_xlines = 1;
    else
        slope_xline = xline_data(n_xline,1);
        intercept_xline = xline_data(n_xline,2);
    end
    
    check_ylines = 0;
    if size(yline_points{n_yline},1) <= 3
        px_yline = px_uFrame(yline_points{n_yline});
        py_yline = py_uFrame(yline_points{n_yline});
        
        slope_yline = point_line_data(p_closest,3);
        intercept_yline = -slope_yline*mean(px_yline) + mean(py_yline);
        
        point_line_data(np,3:4) = [slope_yline, intercept_yline];
        yline_data(n_yline,1:2) = [slope_yline, intercept_yline];
        
        check_ylines = 1;
    else
        slope_yline = yline_data(n_yline,1);
        intercept_yline = yline_data(n_yline,2);
    end
    
    if check_xlines || check_ylines
        px0(np) = (intercept_xline - intercept_yline)/...
            (slope_yline - slope_xline);
        py0(np) = slope_xline*px0(np) + intercept_xline;
    end
end

%% fill in data points outside the video to get full grid

px_grid = NaN(num_x_lines,num_y_lines,vM.nFrames);
py_grid = NaN(num_x_lines,num_y_lines,vM.nFrames);
px0_grid = NaN(num_x_lines,num_y_lines);
py0_grid = NaN(num_x_lines,num_y_lines);

np_real = 0;
real_points = false(num_x_lines,num_y_lines);

for jx = 1:num_y_lines
    for iy = 1:num_x_lines
        guess_px0 = -(yline_data(jx,2) - xline_data(iy,2))/...
            (yline_data(jx,1) - xline_data(iy,1));
        guess_py0 = xline_data(iy,1)*guess_px0 + xline_data(iy,2);

        p_dist_from_guess = sqrt((px0 - guess_px0).^2 + (py0 - guess_py0).^2);
        
        if all(p_dist_from_guess*vM.Calibration > vM.DotSpacing/2)
            px_grid(iy,jx,:) = guess_px0;
            py_grid(iy,jx,:) = guess_py0;
            px0_grid(iy,jx) = guess_px0;
            py0_grid(iy,jx) = guess_py0;
        else
            real_points(iy,jx) = true;
            np_real = np_real + 1;
            px_grid(iy,jx,:) = px(p_dist_from_guess == min(p_dist_from_guess),:);
            py_grid(iy,jx,:) = py(p_dist_from_guess == min(p_dist_from_guess),:);
            px0_grid(iy,jx) = px0(p_dist_from_guess == min(p_dist_from_guess));
            py0_grid(iy,jx) = py0(p_dist_from_guess == min(p_dist_from_guess));
        end
    end
end

% re-organizing Xgrid and Ygrid so they make sense
[px0_grid,Xsort] = sortrows(px0_grid',1);
px0_grid = px0_grid';
py0_grid = py0_grid(:,Xsort);
[py0_grid,Ysort] = sortrows(py0_grid,1);
px0_grid = px0_grid(Ysort,:);

px_grid = px_grid(Ysort,Xsort,:);
py_grid = py_grid(Ysort,Xsort,:);

real_points = real_points(Ysort,Xsort);

% figure(1)
% plot(px,py,'.k')
% hold on
% for xl = 1:num_x_lines
%     plot(px(xline_points{xl}),py(xline_points{xl}),'-r')
%     pause
% end
% hold off
% plot(px,py,'.k')
% hold on
% for yl = 1:num_y_lines
%     plot(px(yline_points{yl}),py(yline_points{yl}),'-r')
%     pause
% end
