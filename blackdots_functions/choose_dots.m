function [px,py] = choose_dots(img,vM)

dotspacing_px = vM.DotSpacing/vM.Calibration;

img2 = mat2gray(img);
[M,N] = size(img2);

fig_1 = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
ax_1 = axes;
im_1 = imagesc(imadjust(img2));
axis image
axis manual
colormap(ax_1,gray*[1 0 0; 0 0.75 0; 0 0 0])
hold on
p1 = plot(NaN,NaN,'.','Color',[0 0.75 1],'PickableParts','none','markersize',12);
p2 = plot(NaN,NaN,'or','markersize',10,'buttondownfcn',@corner_buttondownfcn);
pat3 = patch('XData',NaN,'YData',NaN,'FaceColor','none','EdgeColor',[0 0.75 1],'PickableParts','none');
hold off

% ax_2 = axes;
% im_2 = imagesc([],'pickableparts','none');
% axis image
% axis manual
% set(ax_2,'color','none','pickableparts','none')
% colormap(ax_2,[0 0.5 1; 0 0 0])
% linkaxes([ax_1,ax_2])

init_point = NaN;
down = 0;
down_corner = 0;
px = [];
py = [];
pxc = [];
pyc = [];
pxc_rot = [];
pyc_rot = [];
num_Xc = 0;
num_Yc = 0;
total_xshift = [];
total_yshift = [];
total_xshift_init = [];
total_yshift_init = [];
corner_num = 0;

R = [cos(-vM.rot_angle) -sin(-vM.rot_angle); sin(-vM.rot_angle) cos(-vM.rot_angle)];
R2 = [cos(vM.rot_angle) -sin(vM.rot_angle); sin(vM.rot_angle) cos(vM.rot_angle)];

% set(fig_1,'WindowButtonDownFcn',@fig_buttondownfcn)
set(im_1,'ButtonDownFcn',@im_buttondownfcn)
set(fig_1,'WindowButtonMotionFcn','')
set(fig_1,'WindowButtonUpFcn',@fig_buttonupfcn)

but_done = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
    'Position',[0.85 0.95 0.15 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
    'String','DONE','Callback','uiresume','Enable','off');
but_calc = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
    'Position',[0.85 0.9 0.15 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
    'String','calc','Callback',@fig_calc,'Enable','off');

skip = false;
but_skip = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
    'Position',[0 0.95 0.15 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
    'String','SKIP','Callback','skip = true; uiresume');
but_reset = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
    'Position',[0 0.9 0.15 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
    'String','RESET','Callback',@fig_reset);

uiwait

if skip
    px = [];
    py = [];
else
    px = pxc;
    py = pyc;
end

close(fig_1)

function corner_buttondownfcn(~, ~, ~)
    down = 0;
    down_corner = 1;
    init_point = ax_1.CurrentPoint(1,1:2);
    pxcorners = pxc([1 end],[1 end]);
    pycorners = pyc([1 end],[1 end]);
    [~,corner_num] = min(vecnorm(([pxcorners(:) pycorners(:)] - init_point)'));
    
    total_xshift_init = total_xshift;
    total_yshift_init = total_yshift;
    
    set(fig_1,'WindowButtonMotionFcn',@fig_buttonmotionfcn)
end

function im_buttondownfcn(~, ~, ~)
    if ~down_corner
        init_point = ax_1.CurrentPoint(1,1:2);
        
        down = 1;

        pxc = init_point(1);
        pyc = init_point(2);
        set(p1,'XData',pxc,'YData',pyc)
        set(p2,'XData',NaN,'YData',NaN)

        set(fig_1,'WindowButtonMotionFcn',@fig_buttonmotionfcn)
    end
end

function fig_buttonmotionfcn(~, ~, ~)
    if down
        curr_point = ax_1.CurrentPoint(1,1:2);
        
        curr_point_rot = (R*(curr_point - [N/2,M/2])')' + [vM.N/2, vM.M/2];
        init_point_rot = (R*(init_point - [N/2,M/2])')' + [vM.N/2, vM.M/2];
        
        num_Xc = abs(round((curr_point_rot(1) - init_point_rot(1))/dotspacing_px)) + 1;
        num_Yc = abs(round((curr_point_rot(2) - init_point_rot(2))/dotspacing_px)) + 1;

        if num_Xc > 1
            xgridc_rot = linspace(init_point_rot(1),curr_point_rot(1),num_Xc);
        else
            xgridc_rot = init_point_rot(1);
        end
        if num_Yc > 1
            ygridc_rot = linspace(init_point_rot(2),curr_point_rot(2),num_Yc);
        else
            ygridc_rot = init_point_rot(2);
        end
        
        [pxc_rot,pyc_rot] = meshgrid(xgridc_rot,ygridc_rot);
        
        temp_pc = (R2*([pxc_rot(:), pyc_rot(:)] - [N/2,M/2])')' + [vM.N/2, vM.M/2];
        pxc = temp_pc(:,1);
        pyc = temp_pc(:,2);
        pxc = reshape(pxc,num_Yc,num_Xc);
        pyc = reshape(pyc,num_Yc,num_Xc);
        
        set(p1,'XData',pxc(:),'YData',pyc(:))
        
        total_xshift = zeros(size(pxc));
        total_yshift = zeros(size(pyc));
        pxcorners = pxc([1 end],[1 end]);
        pycorners = pyc([1 end],[1 end]);
        set(p2,'XData',pxcorners(:),'YData',pycorners(:))
    elseif down_corner
        curr_point = ax_1.CurrentPoint(1,1:2);
        
        curr_point_rot = (R*(curr_point - [N/2,M/2])')' + [vM.N/2, vM.M/2];
        init_point_rot = (R*(init_point - [N/2,M/2])')' + [vM.N/2, vM.M/2];
        
        pxcorners_rot = pxc_rot([1 end],[1 end]);
        pycorners_rot = pyc_rot([1 end],[1 end]);

        xshift = zeros(size(pxcorners_rot));
        yshift = zeros(size(pycorners_rot));
        if corner_num == 1
            xshift(1,1) = curr_point_rot(1) - init_point_rot(1) - total_xshift(1,1) + total_xshift_init(1,1);
            yshift(1,1) = curr_point_rot(2) - init_point_rot(2) - total_yshift(1,1) + total_yshift_init(1,1);
        elseif corner_num == 2
            xshift(end,1) = curr_point_rot(1) - init_point_rot(1) - total_xshift(end,1) + total_xshift_init(end,1);
            yshift(end,1) = curr_point_rot(2) - init_point_rot(2) - total_yshift(end,1) + total_yshift_init(end,1);
        elseif corner_num == 3
            xshift(1,end) = curr_point_rot(1) - init_point_rot(1) - total_xshift(1,end) + total_xshift_init(1,end);
            yshift(1,end) = curr_point_rot(2) - init_point_rot(2) - total_yshift(1,end) + total_yshift_init(1,end);
        elseif corner_num == 4
            xshift(end,end) = curr_point_rot(1) - init_point_rot(1) - total_xshift(end,end) + total_xshift_init(end,end);
            yshift(end,end) = curr_point_rot(2) - init_point_rot(2) - total_yshift(end,end) + total_yshift_init(end,end);
        end
        xshift_2 = interp2(pxcorners_rot,pycorners_rot,xshift,pxc_rot,pyc_rot);
        total_xshift = total_xshift + xshift_2;
        yshift_2 = interp2(pxcorners_rot,pycorners_rot,yshift,pxc_rot,pyc_rot);
        total_yshift = total_yshift + yshift_2;
        
        pxc_rot_shift = pxc_rot + total_xshift;
        pyc_rot_shift = pyc_rot + total_yshift;
        
        temp_pc = (R2*([pxc_rot_shift(:), pyc_rot_shift(:)] - [N/2,M/2])')' + [vM.N/2, vM.M/2];
        pxc = temp_pc(:,1);
        pyc = temp_pc(:,2);
        pxc = reshape(pxc,num_Yc,num_Xc);
        pyc = reshape(pyc,num_Yc,num_Xc);
        
        set(p1,'XData',pxc(:),'YData',pyc(:))
        
        pxcorners = pxc([1 end],[1 end]);
        pycorners = pyc([1 end],[1 end]);
        set(p2,'XData',pxcorners(:),'YData',pycorners(:))
    end
end

function fig_buttonupfcn(~, ~, ~)
    if down || down_corner
        set(fig_1,'WindowButtonMotionFcn','')

        set(but_calc,'Enable','on')
        set(but_done,'Enable','off')
        down = 0;
        down_corner = 0;
    end
end

function fig_calc(~, ~, ~)
    [pxc,pyc,~,img_bw] = find_centroids(pxc,pyc,img,vM);
    
    if isempty(pxc)
        skip = true;
        return
    end
    
    set(but_calc,'Enable','off')
    set(but_done,'Enable','on')
    
    set(p1,'XData',pxc(:),'YData',pyc(:))
    
    pxcorner = pxc([1 end],[1 end]);
    pycorner = pyc([1 end],[1 end]);
    set(p2,'XData',pxcorner(:),'YData',pycorner(:))
    
    img_bw_bndry = bwboundaries(img_bw);
%     p_temp = [];
    p_temp = NaN(max(cellfun(@length,img_bw_bndry)),2*length(img_bw_bndry));
    for i = 1:length(img_bw_bndry)
%         p_temp = [p_temp; img_bw_bndry{i}; NaN(1,2)];
        p_temp(1:length(img_bw_bndry{i}),2*i-1:2*i) = img_bw_bndry{i};
    end
    
    set(pat3,'XData',p_temp(:,2:2:end),'YData',p_temp(:,1:2:end));
end

function fig_reset(~, ~, ~)
    temp_pc = (R2*([pxc_rot(:), pyc_rot(:)] - [N/2,M/2])')' + [vM.N/2, vM.M/2];
    pxc = temp_pc(:,1);
    pyc = temp_pc(:,2);
    pxc = reshape(pxc,num_Yc,num_Xc);
    pyc = reshape(pyc,num_Yc,num_Xc);

    set(p1,'XData',pxc(:),'YData',pyc(:))
    
    pxcorners = pxc([1 end],[1 end]);
    pycorners = pyc([1 end],[1 end]);
    set(p2,'XData',pxcorners(:),'YData',pycorners(:))
    
    total_xshift = zeros(size(pxc));
    total_yshift = zeros(size(pyc));
end
end