function BD = calc_dot_size_spacing(px,py,real_points,img,meta)
%all units are pixels except BD.DotSizes and BD.DotSpacings

% find dot sizes
img_1 = imbinarize(imcomplement(img));
img_2 = imfill(img_1,'holes');
dchar = regionprops(bwselect(img_2,px(real_points),py(real_points),4),'EquivDiameter','Eccentricity','Area','Circularity','ConvexHull'); %,'MaxFeretProperties'
% 
% img_edge = imfill(edge(imcomplement(img),'Sobel'),'holes');
% dchar_edge = regionprops(bwselect(img_edge,px(real_points),py(real_points),4),'EquivDiameter','Eccentricity','Area');

fields = fieldnames(dchar);
for i = 1:length(fields)
%     if strcmp(fields{i},'EquivDiameter')
%         BD.(fields{i}) = [dchar_edge.(fields{i})];
%     else
%         BD.(fields{i}) = [dchar.(fields{i})];
%     end
    if strcmp(fields{i},'ConvexHull')
        ConvexPerimeter = cellfun(@(x) sum(sqrt(sum(diff(x,1).^2,2))),{dchar.(fields{i})});
    else
        BD.(fields{i}) = [dchar.(fields{i})];
    end
end
BD.CircularityKB2 = 4*BD.Area./(ConvexPerimeter.^2/pi);
BD.DotSizes = BD.EquivDiameter*meta.Calibration;

% find dot spacings
rn = NaN(size(real_points));
rn(real_points) = 1;
dsp_x = sqrt(diff(px.*rn,1,2).^2 + diff(py.*rn,1,2).^2);
dsp_y = sqrt(diff(px.*rn,1,1).^2 + diff(py.*rn,1,1).^2);
dot_spacings = [dsp_x(~isnan(dsp_x)); dsp_y(~isnan(dsp_y))];
BD.DotSpacings = dot_spacings'*meta.Calibration;

