function outline = get_cell_boundary(image_raw)
% cell boundary measurement

image_raw = mat2gray(image_raw);

[M,N] = size(image_raw);

prev_CB = cell(0);
nCB = 0;
more_cells = true;

while more_cells
fig = figure('units','normalized','position',[0.1 0.1 0.8 0.7]);
ax = axes;
img_show = imagesc('CData',imadjust(image_raw,stretchlim(image_raw,[0.01 0.999])),'PickableParts','none');
colormap(gray)
axis image
set(ax,'YDir','reverse')
set(ax,'CLim',[min(img_show.CData(:)), max(img_show.CData(:))])

% % r = min([M,N])/100;
r = 0;
% set(ax,'Units','Points')
% rp = min(ax.Position(3:4))/100;
% set(ax,'Units','normalized')

outline = [];
hold(ax,'on')
p0 = patch('XData',[],'YData',[],'FaceColor','none','EdgeColor','w','linewidth',2,'PickableParts','none');
p1 = plot(NaN,NaN,'-r','PickableParts','none');
% p2 = plot(NaN,NaN,'or','PickableParts','none','markersize',2*rp);
p2 = patch('XData',[],'YData',[],'FaceColor','r','FaceAlpha',0.2,'EdgeColor','r','PickableParts','none');
p3 = patch('XData',[],'YData',[],'FaceColor','y','FaceAlpha',0.2,'EdgeColor','y','PickableParts','none');
hold(ax,'off')

if ~isempty(prev_CB)
    [~,bndry_order] = sort(cellfun(@length,prev_CB),'descend');
    prev_CB = prev_CB(bndry_order);
    max_bndry_length = max(cellfun(@length,prev_CB))+1;
    bndry = NaN(max_bndry_length,2*length(prev_CB));
    for i = 1:length(prev_CB)
        bndry(:,2*i-1:2*i) = repmat(prev_CB{i}(1,:),max_bndry_length,1); % (1,[2 1])
        bndry(1:length(prev_CB{i}),2*i-1:2*i) = prev_CB{i}; % (:,[2 1])
    end
    set(p3,'XData',bndry(:,1:2:end),'YData',bndry(:,2:2:end))
end

set(fig,'WindowButtonDownFcn',@fig_buttondownfcn)
set(fig,'WindowButtonUpFcn',@fig_buttonupfcn)
movept = 1;
% moved = false;

but_done = uicontrol('Units','normalized','Position',[0.9 0.95 0.1 0.05],...
    'style','pushbutton','String','done','callback','uiresume');

but_reset = uicontrol('Units','normalized','Position',[0 0.95 0.1 0.05],...
    'style','pushbutton','String','reset','callback',@callback_reset);

uiwait

close(fig)

% close the loop:
if ~isempty(outline)
    outline = [outline; outline(1,:)];

    nCB = nCB + 1;
    prev_CB{nCB} = outline;
else
    more_cells = false;
end

end

outline = prev_CB;

    function fig_buttondownfcn(obj, ~, ~)
        if ~isa(obj.CurrentObject,'matlab.graphics.axis.Axes')
            return
        end
        if ~isempty(outline)
            dist = sqrt(sum((outline - ax.CurrentPoint(1,1:2)).^2,2));
        else
            dist = inf;
        end
        
        r = min([diff(ax.XLim),diff(ax.YLim)])/50;

        if all(dist > r)
            outline = [outline; ax.CurrentPoint(1,1:2)];
        else
            if strcmp(get(obj,'SelectionType'),'open') && (dist(1) < r)
                sgf_cb = bwboundaries(poly2mask(outline(:,1),outline(:,2),M,N),8);
                sgf_cb = sgf_cb{1}(:,[2 1]);

                sgf_width = length(sgf_cb)/50;
                sgf_order = 2;

                sgf_width = ceil(sgf_width*1.1);
                sgf_width = sgf_width + mod(sgf_width,2) + 1;
                sgf_cb2 = sgolayfilt(repmat(sgf_cb,3,1),sgf_order,sgf_width);
                sgf_cb2 = sgf_cb2(size(sgf_cb,1)+1:2*size(sgf_cb,1),:);
                set(p1,'XData',sgf_cb2(:,1),'YData',sgf_cb2(:,2));
                set(p2,'XData',[],'YData',[]);
                outline = sgf_cb2;
                return
            else
                set(fig,'WindowButtonMotionFcn',@fig_buttonmotionfcn)
                movept = find(dist < r);
            end
%             if dist(1) < r
%                 uiresume
%             else
%                 set(fig,'WindowButtonMotionFcn',@fig_buttonmotionfcn)
%                 movept = find(dist < r);
%             end
        end

        set(p1,'XData',outline(:,1),'YData',outline(:,2))
%         set(p2,'XData',outline(:,1),'YData',outline(:,2))
        xcirc = outline(:,1) + r*cos(linspace(0,2*pi,10));
        ycirc = outline(:,2) + r*sin(linspace(0,2*pi,10));
        set(p2,'XData',xcirc','YData',ycirc')
        set(p0,'XData',xcirc(1,:)','YData',ycirc(1,:)')
    end
    function fig_buttonmotionfcn(~, ~, ~)
%         moved = true;
        outline(movept,:) = ax.CurrentPoint(1,1:2);

        set(p1,'XData',outline(:,1),'YData',outline(:,2))
%         set(p2,'XData',outline(:,1),'YData',outline(:,2))
        xcirc = outline(:,1) + r*cos(linspace(0,2*pi,10));
        ycirc = outline(:,2) + r*sin(linspace(0,2*pi,10));
        set(p2,'XData',xcirc','YData',ycirc')
        set(p0,'XData',xcirc(1,:)','YData',ycirc(1,:)')
    end
    function fig_buttonupfcn(~, ~, ~)
        set(fig,'WindowButtonMotionFcn','')
%         if dist(1) < r && ~moved
%             uiresume
%         end
%         moved = false;
    end
    function callback_reset(~, ~, ~)
        outline = [];
        set(p1,'XData',NaN,'YData',NaN)
%         set(p2,'XData',NaN,'YData',NaN)
        set(p2,'XData',[],'YData',[])
        set(p0,'XData',[],'YData',[])
    end
end