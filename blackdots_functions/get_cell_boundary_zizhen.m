function [stats] = get_cell_boundary_zizhen(image_raw,vM)
% Zizhen method:
cellFrameImage = image_raw;
[x3, y3] = size(cellFrameImage);
image2 = zeros(x3, y3, 1);
threshold = 0;
if vM.Calibration < 0.048
    threshold = vM.Calibration * 10000;
else 
    threshold = 440;
end
%threshold = 210;
for i = 1: x3
    for j = 1: y3
        if cellFrameImage(i, j)>threshold
            image2(i, j ,1) = 1;
        else
            image2(i, j, 1) = 0;
        end
    end
end
f = figure;
f.Position = [10 350 1000 1000];
%filling holes in image
image2 = imfill(image2,'holes');
min_cell_size = 5;  %value in um^2, default is 5um^2
min_cell_threshold = round(min_cell_size/ vM.Calibration / vM.Calibration); 
seD = strel('diamond',1);
image2 = imerode(image2, seD);
BWfinal = imerode(image2, seD);
BWfinal = bwareaopen(BWfinal,min_cell_threshold);
imshow(BWfinal);
drawnow;

fig_0 = uifigure('Position',[10 10 650 275]);
sld = uislider(fig_0,...
               'Position',[20 60 500 4],...
               'ValueChangingFcn',...
               @(sld,event) imageFunct(fig_0,event,image2,x3,y3, cellFrameImage, vM.Calibration),...
               'Limits', [0 500], 'Value', threshold);

% sld = uicontrol('Parent',f,'Style','slider',...
%                  'Units','normalized',...
%                  'Position',[0.3 0.5 0.4 0.1],...
%                  'Tag','slider1',...
%                  'Callback',@(sld,event) this.imageFunct(f,event,image2,x3,y3, X, celllabelframe));

%'Value', true,...'Position',[0.8 0.95 0.2 0.05],...
btn_done_0 = uibutton(fig_0,...
               'Text', 'Done',...
               'ButtonPushedFcn', 'uiresume(f)');

% but_done_0 = uicontrol('parent',f,'units','normalized','Style','pushbutton',...
%     'Position',[0.8 0.95 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
%     'String','DONE','Callback','uiresume(f)');
uiwait(f);
close(f);




n_threshold = getappdata(fig_0, 'threshold');
if n_threshold ~= threshold
    BWfinal = getappdata(fig_0,'image');

else
end
close(fig_0);


labeled = logical(BWfinal);
stats = regionprops(labeled,{'Area','Centroid'});
stats = struct2table(stats);
stats.Area_um = stats.Area *vM.Calibration * vM.Calibration;
stats.index = (1:height(stats)).'
stats = stats(:,[ 4 2 3 1]);

[bwlabeled, num_of_cells] = bwlabel(BWfinal);
stats2 = cell(num_of_cells, 9); 
for i = 1:(num_of_cells * 9)
    stats2{i} = 0;
end
for i = 1 : num_of_cells
    stats2{i,1} =i;
end
for i = 1: x3
    for j = 1 : y3
        if bwlabeled(i,j) ~= 0
            stats2{bwlabeled(i,j), 2} = stats2{bwlabeled(i,j), 2} +1;
        end
    end
end

fig_BW = figure;

imshow(BWfinal)
hold on
for kk = 1:height(stats)
  text(stats.Centroid(kk,1)+10, stats.Centroid(kk,2),...
      num2str(stats.Area_um(kk)), 'Color', 'g')
end

% %     savefig([save_folder filesep 'selected_cell_area.fig']);
% Calculating total area of cells
area_in_pixels = 0;
for i = 1: x3
    for j = 1: y3
        if BWfinal(i,j,1) == 1
            area_in_pixels = area_in_pixels + 1;
        end
    end
end

cell_area = area_in_pixels * vM.Calibration * vM.Calibration;

end

%%
function imageFunct(fig, event, image, x3, y3, cellFrameImage, cali)
    threshold = event.Value;
    for i = 1: x3
        for j = 1: y3
            if cellFrameImage(i, j)>threshold
                image(i, j ,1) = 1;
            else
                image(i, j, 1) = 0;
            end
        end
    end
    
    % filling holes in image
    image = imfill(image,'holes');

    min_cell_size = 5;  %value in um^2, default is 5um^2
    min_cell_threshold = round(min_cell_size/ cali / cali); 


    seD = strel('diamond',1);
    image = imerode(image, seD);
    image = imerode(image, seD);
    BWFinal = bwareaopen(image,min_cell_threshold);
    
    imshow(BWFinal);
    
    setappdata(fig, 'image', BWFinal);
    setappdata(fig, 'threshold', threshold);
    %fig.image = image;
    %fig.threshold = threshold;
end