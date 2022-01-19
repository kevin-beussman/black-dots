function [rot_angle] = get_rot_from_gridpts(input)
% input is a grid of points, just get slope of each line to find rotation
if ndims(input) == 3
    px = input(:,:,1);
    py = input(:,:,2);
    
    [num_Y,num_X] = size(px);
    
    rot_angle_x = zeros(num_Y,1);
    for iy = 1:num_Y
        jx = 1:num_X;
        pca_coeff = pca([px(iy,jx)', py(iy,jx)']);
        slope = pca_coeff(2,1)/pca_coeff(1,1);
        rot_angle_x(iy) = atan(slope);
    end

    rot_angle_y = zeros(num_X,1);
    for jx = 1:num_X
        iy = 1:num_Y;
        pca_coeff = pca([px(iy,jx), py(iy,jx)]);
        slope = pca_coeff(2,1)/pca_coeff(1,1);
        rot_angle_y(jx) = atan(slope) - pi/2;
    end
    
    if ~all(sign(rot_angle_x) == sign(rot_angle_x(1)))
        rot_angle_x(rot_angle_x < 0) = rot_angle_x(rot_angle_x < 0) + 2*pi;
    end
    
%     rot_angle = mean([rot_angle_x; rot_angle_y]);
    rot_angle = mean(rot_angle_x);
end