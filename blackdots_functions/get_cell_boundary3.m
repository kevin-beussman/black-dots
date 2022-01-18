function [cb] = get_cell_boundary3(image_raw,vM)

% image_raw = mat2gray(image_raw);
image_use = mat2gray(medfilt2(image_raw,[20 20],'symmetric'));

[M,N] = size(image_raw);

fig = figure('units','normalized','position',[0.1 0.1 0.8 0.7]);
ax = axes(fig,'units','normalized','position',[0.05 0.15 0.9 0.8]);
img_show = imagesc(image_raw);
colormap(gray)
axis image
set(ax,'YDir','reverse')
set(ax,'CLim',[min(img_show.CData(:)), max(img_show.CData(:))])

threshold = 0.5;
min_cell_size = 30;  %value in um^2, default is 5um^2
min_cell_threshold = round(min_cell_size/ vM.Calibration / vM.Calibration); 
cell_boundaries = get_boundary(image_use,threshold,min_cell_threshold);
set(fig,'UserData',cell_boundaries);
% set(fig,'UserData',cell2mat(cellfun(@(v,n) [v; nan(max(cellfun('length',fig.UserData))-n,2)],fig.UserData,num2cell(cellfun('length',fig.UserData)),'UniformOutput',false)'));

pat_show = patch('XData',fig.UserData(:,1:2:end),'YData',fig.UserData(:,2:2:end),...
    'FaceVertexCData',hsv(size(fig.UserData,2)/2),'FaceAlpha',0.5,'FaceColor','flat','EdgeColor','w','linewidth',1);

sld_threshselector = uicontrol('parent',fig,'units','normalized','Style','slider',...
    'Position',[0.1 0.05 0.68 0.05],'Min',0,'Max',1,'Value',threshold);
txt_threshselector = uicontrol('parent',fig,'units','normalized','Style','edit',...
    'Position',[0.8 0.05 0.05 0.05],'HorizontalAlignment','left','FontUnits','normalized','FontSize',0.5,...
    'String',sprintf('%0.3f',sld_threshselector.Value),...
    'Callback','sld_threshselector.Value = str2double(txt_threshselector.String);');

% callbacks for slider
fcn1 = @(a,b) set(txt_threshselector,'String',sprintf('%0.3f',sld_threshselector.Value));
fcn2 = @(a,b) set(fig,'UserData',get_boundary(image_use,sld_threshselector.Value,min_cell_threshold));
% fcn3 = @(a,b) set(fig,'UserData',cell2mat(cellfun(@(v,n) [v; nan(max(cellfun('length',fig.UserData))-n,2)],fig.UserData,num2cell(cellfun('length',fig.UserData)),'UniformOutput',false)'));
fcn4 = @(a,b) set(pat_show,'XData',fig.UserData(:,1:2:end),'YData',fig.UserData(:,2:2:end));
fcn5 = @(a,b) set(pat_show,'FaceVertexCData',hsv(size(fig.UserData,2)/2),'FaceColor','flat');
fcn = @(a,b) cellfun(@feval,{fcn1 fcn2 fcn4 fcn5});

addlistener(sld_threshselector,'Value','PostSet',fcn); % real-time slider

% set(sld_threshselector,'Callback',fcn)

txt_instructions = uicontrol('parent',fig,'units','normalized','Style','text',...
    'Position',[0 0.95 0.9 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
    'String','Please select threshold to get cell boundary: ');
but_done = uicontrol('parent',fig,'units','normalized','Style','pushbutton',...
    'Position',[0.9 0.95 0.1 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
    'String','DONE','Callback','uiresume');

uiwait

cell_boundary = get_boundary(image_use,sld_threshselector.Value,min_cell_threshold);

nCells = size(cell_boundary,2)/2;
for nc = 1:nCells
    temp_cb = cell_boundary(:,nc*2-1:nc*2);
    while all(temp_cb(end,:) == temp_cb(end-1,:))
        temp_cb(end,:) = [];
    end
    cb{nc} = temp_cb;
end

close(fig)

function bndry = get_boundary(image,threshold,min_size)
%     I1 = adaptthresh(image,threshold,'Statistic','mean');
%     I2 = imbinarize(image,I1);
    I2 = imbinarize(image,threshold);
    I3 = imfill(I2,'holes');
    I4 = bwmorph(I3,'clean');
    I5 = bwmorph(I4,'majority');
    I6 = bwmorph(I5,'clean');
%     I7 = bwselect(I6,round(select_pt(1)),round(select_pt(2)),4);
    I7 = imfill(I6,'holes');
    I8 = medfilt2(I7,[5 5]);
    I9 = bwareaopen(I8,min_size);
    I10 = bwboundaries(I9,8,'noholes');
    if ~isempty(I10)
        [~,bndry_order] = sort(cellfun(@length,I10),'descend');
        I10 = I10(bndry_order);
        max_bndry_length = max(cellfun(@length,I10))+1;
        bndry = NaN(max_bndry_length,2*length(I10));
        for i = 1:length(I10)
            bndry(:,2*i-1:2*i) = repmat(I10{i}(1,[2 1]),max_bndry_length,1);
            bndry(1:length(I10{i}),2*i-1:2*i) = I10{i}(:,[2 1]);
        end
%         bndry = cell2mat(I10);
%         bndry = [bndry; bndry(1,:)];
%         bndry = bndry(:,[2 1]); % reverse order of x and y
    else
        bndry = [];
    end
end

end