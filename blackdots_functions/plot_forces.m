%% plot image with forces for publication
k = meta_BD.uFrame;
s = [0, 0];
for ic = 1:nCells
    if ~isempty(celldata(ic))
        pos = [celldata(ic).Yvector + s(2), celldata(ic).Xvector + s(1)];
        loc = [celldata(ic).Yloc_k(:,k) + s(2), celldata(ic).Xloc_k(:,k) + s(1)];
        disp = [celldata(ic).Ydisp_k(:,k), celldata(ic).Xdisp_k(:,k)];
        forc = [celldata(ic).Yforce_k(:,k), celldata(ic).Xforce_k(:,k)];
        trac = [celldata(ic).Ytrac_k(:,k), celldata(ic).Xtrac_k(:,k)];
        
        set(0,'units','pixels')
        screensize = get(0,'screensize');
        figwidth = 1*min(screensize(3:4));
        fig_disp = figure('units','pixels','position',figwidth*[0.3 0.3 0.4 0.4],'Menu','none','ToolBar','none');
        ax_disp = axes(fig_disp,'units','normalized','position',[0 0 1 1]);
        imagesc(celldata(ic).img_ref,[min(celldata(ic).img_ref(:)), max(celldata(ic).img_ref(:))])
        colormap(gray*[1 0 0;0 130/255 0;0 0 0])
        axis image
        axis manual
        axis off
        fig_disp.Position(4) = celldata(ic).M/celldata(ic).N*fig_disp.Position(3);
        hold on
        plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'-w');
        
        color_max = [1 1 0];
        color_min = [1 0 1];
        
        arrow_start = [loc(:,2),loc(:,1)];
        arrow_end = arrow_start + arrowscale*[forc(:,2), forc(:,1)];
        
        arrowlength = sqrt(sum((arrow_end - arrow_start).^2,2))/2/2;
        arrowwidth = arrowlength*5/15;
        arrow_color = mat2gray(arrowlength).*(color_max - color_min) + color_min;
        q_trac = arrow(arrow_start,arrow_end,...
            'BaseAngle',90,'TipAngle',25,'Width',arrowwidth,'Length',arrowlength,...
            'LineWidth',0.25,'EdgeColor','none','FaceVertexCData',arrow_color,'FaceColor','flat'); % 'FaceColor','w','EdgeColor','none'
        hold off

        % make scalebars
        ax_scalebars = axes(fig_disp,'position',get(ax_disp,'position'),'hittest','off','pickableparts','none');

        ylims = ax_disp.YLim;
        xlims = ax_disp.XLim;
        
        axis image
        axis off
        set(ax_scalebars,'color','none',...
            'XLim',get(ax_disp,'XLim'),...
            'YLim',get(ax_disp,'YLim'),...
            'XLimMode','manual','YLimMode','manual',...
            'YDir','Reverse')
        
        linkaxes([ax_disp, ax_scalebars]);
        
        hold on
        patch(ax_scalebars,'XData',[xlims(2) - scalebar_length/meta_BD.Calibration - 0.2*diff(xlims),xlims(2),xlims(2),xlims(2) - scalebar_length/meta_BD.Calibration - 0.2*diff(xlims)],...
            'YData',[ylims(2),ylims(2),ylims(2) - 0.15*diff(ylims),ylims(2) - 0.15*diff(ylims)],...
            'EdgeColor','none','FaceColor','k','FaceAlpha',1)
        
        plot(ax_scalebars,[xlims(2) - scalebar_length/meta_BD.Calibration - 0.15*diff(xlims), xlims(2) - 0.15*diff(xlims)],...
            [ylims(2) - 0.05*diff(ylims), ylims(2) - 0.05*diff(ylims)],'w','linewidth',5)
        text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.05*diff(ylims),...
            sprintf(['%.0f ' char(181) 'm'],scalebar_length),'color','w','VerticalAlignment','Middle')
        
        arrow([xlims(2) - arrowscale*scalebar_force/0.001 - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
            [xlims(2) - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
            'BaseAngle',90,'TipAngle',25,'Width',max(arrowwidth),'Length',max(arrowlength),'LineWidth',0.25,'FaceColor','w','EdgeColor','none')
        text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
            sprintf('%.0f nN',scalebar_force),'color','w','VerticalAlignment','Middle')
        hold off

        set(gcf, 'InvertHardCopy', 'off')
        print(sprintf('plot_force_cell%i', ic),'-dpng')
        print(sprintf('plot_force_cell%i', ic),'-dsvg','-vector')
    end
end

%% plot rotated image with forces for publication
k = meta_BD.uFrame;
s = [0, 0];
for ic = 1:nCells
    if ~isempty(celldata(ic))
        aa = celldata(ic).rot_angle;
        img_rotated = imrotate(celldata(ic).img_ref,aa*180/pi,'crop');
        
        R = [cos(aa), sin(aa); -sin(aa), cos(aa)];
        % R = [1 0; 0 1];
        centerX = floor(celldata(ic).N/2+1) + s(1);
        centerY = floor(celldata(ic).M/2+1) + s(2);
        pos = [celldata(ic).Yvector + s(2), celldata(ic).Xvector + s(1)];
        loc = [celldata(ic).Yloc_k(:,k) + s(2), celldata(ic).Xloc_k(:,k) + s(1)];
        disp = [celldata(ic).Ydisp_k(:,k), celldata(ic).Xdisp_k(:,k)];
        forc = [celldata(ic).Yforce_k(:,k), celldata(ic).Xforce_k(:,k)];
        CB = celldata(ic).CB(:,[2,1]) + s;
        pos_rot = (pos - [centerY, centerX])*R + [centerY, centerX];
        loc_rot = (loc - [centerY, centerX])*R + [centerY, centerX];
        CB_rot = (CB - [centerY, centerX])*R + [centerY, centerX]; CB_rot = CB_rot(:,[2,1]);
        disp_rot = disp*R;
        forc_rot = forc*R;
        
        set(0,'units','pixels')
        screensize = get(0,'screensize');
        figwidth = 1*min(screensize(3:4));
        fig_disp = figure('units','pixels','position',figwidth*[0.3 0.3 0.4 0.4],'Menu','none','ToolBar','none');
        ax_disp = axes(fig_disp,'units','normalized','position',[0 0 1 1]);
        imagesc(img_rotated,[min(celldata(ic).img_ref(:)), max(celldata(ic).img_ref(:))])
        colormap(gray*[1 0 0;0 130/255 0;0 0 0])
        axis image
        axis manual
        axis off
        fig_disp.Position(4) = celldata(ic).M/celldata(ic).N*fig_disp.Position(3);
        hold on
        plot(CB_rot(:,1) + s(1),CB_rot(:,2) + s(2),'-w');
        % plot(pos_rot(:,2),pos_rot(:,1),'+k','markersize',10)
        % plot(loc_rot(:,2),loc_rot(:,1),'.k','markersize',25)
        % quiver(pos_rot(:,2),pos_rot(:,1),disp_rot(:,2),disp_rot(:,1),0,'-w','linewidth',1)
        % arrow([pos_rot(:,2),pos_rot(:,1)],[pos_rot(:,2),pos_rot(:,1)]+[disp_rot(:,2),disp_rot(:,1)],...
        %     'Length',5)
        
        arrow_start = [loc_rot(:,2),loc_rot(:,1)];
        arrow_end = arrow_start + arrowscale*[forc_rot(:,2), forc_rot(:,1)];
        
        arrowlength = sqrt(sum((arrow_end - arrow_start).^2,2))/2/2;
        arrowwidth = arrowlength*5/15;
        arrow_color = mat2gray(arrowlength).*(color_max - color_min) + color_min;
        q_trac = arrow(arrow_start,arrow_end,...
            'BaseAngle',90,'TipAngle',25,'Width',arrowwidth,'Length',arrowlength,...
            'LineWidth',0.25,'EdgeColor','none','FaceVertexCData',arrow_color,'FaceColor','flat'); % 'FaceColor','w','EdgeColor','none'
        hold off
        
        ax_scalebars = axes(fig_disp,'position',get(ax_disp,'position'),'hittest','off','pickableparts','none');
        
        ylims = ax_disp.YLim;
        xlims = ax_disp.XLim;
        
        axis image
        axis off
        set(ax_scalebars,'color','none',...
            'XLim',get(ax_disp,'XLim'),...
            'YLim',get(ax_disp,'YLim'),...
            'XLimMode','manual','YLimMode','manual',...
            'YDir','Reverse')
        
        linkaxes([ax_disp, ax_scalebars]);
        
        hold on
        patch(ax_scalebars,'XData',[xlims(2) - scalebar_length/meta_BD.Calibration - 0.2*diff(xlims),xlims(2),xlims(2),xlims(2) - scalebar_length/meta_BD.Calibration - 0.2*diff(xlims)],...
            'YData',[ylims(2),ylims(2),ylims(2) - 0.15*diff(ylims),ylims(2) - 0.15*diff(ylims)],...
            'EdgeColor','none','FaceColor','k','FaceAlpha',1)
        
        plot(ax_scalebars,[xlims(2) - scalebar_length/meta_BD.Calibration - 0.15*diff(xlims), xlims(2) - 0.15*diff(xlims)],...
            [ylims(2) - 0.05*diff(ylims), ylims(2) - 0.05*diff(ylims)],'w','linewidth',5)
        text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.05*diff(ylims),...
            sprintf(['%.0f ' char(181) 'm'],scalebar_length),'color','w','VerticalAlignment','Middle')
        
        arrow([xlims(2) - arrowscale*scalebar_force/0.001 - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
            [xlims(2) - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
            'BaseAngle',90,'TipAngle',25,'Width',max(arrowwidth),'Length',max(arrowlength),'LineWidth',0.25,'FaceColor','w','EdgeColor','none')
        text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
            sprintf('%.0f nN',scalebar_force),'color','w','VerticalAlignment','Middle')
        hold off

        set(gcf, 'InvertHardCopy', 'off')
        print(sprintf('plot_force_rot_cell%i', ic),'-dpng')
        print(sprintf('plot_force_rot_cell%i', ic),'-dsvg','-vector')
    end
end

% EVERYTHING BELOW HERE MIGHT NOT WORK YET

%% plot all cell forces on whole image
% figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% imagesc(img_REFBD)
% colormap(gray*[1 0 0;0 130/255 0;0 0 0])
% axis image
% axis manual
% axis off
% hold on
% nCells = length(celldata);
% k = meta_BD.uFrame;
% for ic = 1:nCells
%     if ~isempty(celldata(ic))
%         s = celldata(ic).crop([1 2]) - [1, 1];
%         p_bd = plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'-c','linewidth',2);
%         plot(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),'.w','markersize',8)
%         plot(celldata(ic).px_k(:,:,k) + s(1),celldata(ic).py_k(:,:,k) + s(2),'.w','markersize',8)
%         plot(celldata(ic).Xvector(celldata(ic).celldots) + s(1),celldata(ic).Yvector(celldata(ic).celldots) + s(2),'ow')
%         quiver(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),celldata(ic).Xtrac_k(:,k),celldata(ic).Ytrac_k(:,k),0.5,'-c','linewidth',1)
%     end
% end
% hold off


%% plot each cell's forces separately
% nCells = length(celldata);
% k = meta_BD.uFrame;
% for ic = 1:nCells
%     if ~isempty(celldata(ic))
%         figure('units','normalized','position',[0.3 0.3 0.4 0.4]);
%         imagesc(celldata(ic).img_ref)
%         colormap(gray*[1 0 0;0 130/255 0;0 0 0])
%         axis image
%         axis manual
%         axis off
%         hold on
%         s = [0, 0];
%         p_bd = plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'-c','linewidth',2);
%         plot(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),'.w','markersize',8)
%         plot(celldata(ic).px_k(:,:,k) + s(1),celldata(ic).py_k(:,:,k) + s(2),'.w','markersize',8)
%         plot(celldata(ic).Xvector(celldata(ic).celldots) + s(1),celldata(ic).Yvector(celldata(ic).celldots) + s(2),'ow')
%         quiver(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),celldata(ic).Xtrac_k(:,k),celldata(ic).Ytrac_k(:,k),0.5,'-c','linewidth',1)
%         hold off
%     end
% end

%% traction heatmap
% ic = 1;
% % k = 123;
% k = 1;
% ninterp = 1;
% 
% trac_mag = sqrt(celldata(ic).Xtrac_k(:,k).^2 + celldata(ic).Ytrac_k(:,k).^2);
% [tx,ty] = meshgrid(linspace(0,celldata(ic).N,size(celldata(ic).Xgrid,2)*ninterp),linspace(0,celldata(ic).M,size(celldata(ic).Xgrid,1)*ninterp));
% trac_f = scatteredInterpolant(celldata(ic).Xvector,celldata(ic).Yvector,trac_mag,'natural');
% trac_2 = trac_f(tx,ty);
% 
% fig_ex = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% imagesc(tx(:),ty(:),trac_2)
% colormap(jet)
% axis image
% axis manual
% axis off
% hold on
% % if ~isempty(celldata(ic).CB)
% %     nCells = length(celldata(ic).CB);
% %     for ic = 1:nCells
% %         p_bd = plot(celldata(ic).CB(:,1),celldata(ic).CB(:,2),'-w','linewidth',2);
% %     end
% % end
% p_bd = plot(celldata(ic).CB(:,1),celldata(ic).CB(:,2),'-w','linewidth',2);

%% interactive plot -- calculate specific forces (not working yet)
% fig_1 = figure;
% ax_1 = axes;
% im_1 = imagesc(MergedImage); 
% title('Select Points That the Cell Touches')
% hold on
% p_dot = plot(px(:,meta_BD.uFrame),py(:,meta_BD.uFrame),'.r','markersize',10); % plot x y positions found above
% p2_dot = plot(px(1,meta_BD.uFrame),py(1,meta_BD.uFrame),'or'); 
% p_dot_calc = plot(0,0,'.k','markersize',10);
% hold off
% 
% set(im_1,'hittest','off')
% set(p_dot,'hittest','off')
% set(p2_dot,'hittest','off')
% 
% fcn2_1 = @(a,b) set(fig_1,'UserData',find(sqrt((px(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,2)).^2) == min(sqrt((px(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,2)).^2))));
% fcn2_5 = @(a,b) set(p2_dot,'xdata',px(fig_1.UserData,meta_BD.uFrame),'ydata',py(fig_1.UserData,meta_BD.uFrame));
% fcn2 = @(a,b) cellfun(@feval,{fcn2_1 fcn2_5});
% set(ax_1,'buttondownfcn',fcn2)
% 
% % auto-point selection 
% dist_between_poles_pixels = round(1/(2*Calibration)); 
% th_val = 0; % increase this value to include more nearby cells
% p_calc = false(size(px));
% points_cnt = size(px);
% for i = 1: points_cnt
%     x = round(px(i));
%     y = round(py(i));
%     if image2(y, x, 1) == 1
%         p_calc(i) = true;
%     else % checking for if post is close to cell but not touched directly
%         for a = 1 : (dist_between_poles_pixels * 2+ th_val)
%             for b = 1 : (dist_between_poles_pixels * 2+ th_val)
%                 try
%                     x_mod = x - dist_between_poles_pixels + a;
%                     y_mod = y - dist_between_poles_pixels + b;
%                     if image2(x_mod, y_mod,1) == 1
%                         if ((x_mod - x)*(x_mod - x)+(y_mod - y)*(y_mod - y)) < ((dist_between_poles_pixels + th_val/2)  * (dist_between_poles_pixels + th_val/2))
%                             p_calc(i) = true;
%                         end
%                     end
%                 catch
%                 end
%             end
%         end
%     end 
% end
% 
% set(p_dot_calc,'xdata',px(p_calc,meta_BD.uFrame),'ydata',py(p_calc,meta_BD.uFrame));
% 
% but_calc = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
%     'Position',[0.8 0 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
%     'String','CALC','Callback','p_calc(fig_1.UserData) = true; set(p_dot_calc,''xdata'',px(p_calc,meta_BD.uFrame),''ydata'',py(p_calc,meta_BD.uFrame));');
% 
% but_done = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
%     'Position',[0.8 0.95 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
%     'String','DONE','Callback','uiresume');
% 
% uiwait
% 
% if sum(p_calc)>0 % if any points are selected
%     def_calc = dot_displacements(p_calc,:);
%     def_noise = dot_displacements(~p_calc,:);
%     def_calc_mean = mean(dot_displacements(p_calc,:))
%     def_noise_mean = mean(dot_displacements(~p_calc,:));
% end

%% Plot results over time



% set(0,'units','pixels')
% screensize = get(0,'screensize');
% figwidth = 1*min(screensize(3:4));
% fig_disp = figure('units','pixels','position',figwidth*[0.3 0.3 0.4 0.4],'Menu','none','ToolBar','none');
% ax_disp = axes(fig_disp,'units','normalized','position',[0 0 1 1]);
% imagesc(celldata(ic).img_ref,[min(celldata(ic).img_ref(:)), max(celldata(ic).img_ref(:))])
% colormap(gray*[1 0 0;0 130/255 0;0 0 0])
% axis image
% axis manual
% axis off
% fig_disp.Position(4) = celldata(ic).M/celldata(ic).N*fig_disp.Position(3);
% hold on
% plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'-w');
% 
% color_max = [1 1 0];
% color_min = [1 0 1];
% cellboundary = [];

if meta_BD.nFrames >= 2
    s = [0, 0];
    ic = cell_to_plot;
    k = 1;
    
    xrange = celldata(ic).crop([2,4]); xrange(2) = round(xrange(2) + xrange(1) - 1); xrange(1) = round(xrange(1));
    yrange = celldata(ic).crop([1,3]); yrange(2) = round(yrange(2) + yrange(1) - 1); yrange(1) = round(yrange(1));
    
    img_BD_crop = img_BD(xrange(1):xrange(2),yrange(1):yrange(2),:);
    % autocontrast = stretchlim(img_BD_crop(:,:,meta_BD.uFrame));
%     autocontrast = stretchlim(img_BD_crop(:));
    autocontrast = [min(min(img_BD_crop(:,:,1))), max(max(img_BD_crop(:,:,1)))];
    
    fig_video = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
    img_video = imagesc(img_BD_crop(:,:,1), autocontrast);
    axis image
    axis manual
    colormap(gray)
    hold on
    % s = celldata(ic).crop([1 2]);
    
    pos = [celldata(ic).Yvector + s(2), celldata(ic).Xvector + s(1)];
    loc = [celldata(ic).Yloc_k(:,k) + s(2), celldata(ic).Xloc_k(:,k) + s(1)];
    loc_filt = [celldata(ic).Yloc_k_filt(:,k) + s(2), celldata(ic).Xloc_k_filt(:,k) + s(1)];
    disp = [celldata(ic).Ydisp_k(:,k), celldata(ic).Xdisp_k(:,k)];
    disp_filt = [celldata(ic).Ydisp_k_filt(:,k), celldata(ic).Xdisp_k_filt(:,k)];
    forc = [celldata(ic).Yforce_k(:,k), celldata(ic).Xforce_k(:,k)];
    trac = [celldata(ic).Ytrac_k(:,k), celldata(ic).Xtrac_k(:,k)];
    
    % plot(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),'.w','markersize',8)
    % plot(celldata(ic).Xvector(celldata(ic).celldots) + s(1),celldata(ic).Yvector(celldata(ic).celldots) + s(2),'ow')
    % plot(celldata(ic).Xvector(celldata(ic).real_points) + s(1),celldata(ic).Yvector(celldata(ic).real_points) + s(2),'.w')
    p_pos = plot(loc_filt(:,2),loc_filt(:,1),'.w');
    p_bd = plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'--g');
%     q_disp = quiver(loc_filt(:,2),loc_filt(:,1),disp_filt(:,2),disp_filt(:,1),1,'-w','linewidth',1);
    q_forc = quiver(loc_filt(:,2),loc_filt(:,1),arrowscale*forc(:,2),arrowscale*forc(:,1),0,'-c','linewidth',1);
    tl = title('0');
    hold off
    set(gca,'YDir','reverse')
    
    k = 1;
    while true
        
        if k > meta_BD.nFrames
            k = 1;
        end
    
        pos = [celldata(ic).Yvector + s(2), celldata(ic).Xvector + s(1)];
        loc = [celldata(ic).Yloc_k(:,k) + s(2), celldata(ic).Xloc_k(:,k) + s(1)];
        loc_filt = [celldata(ic).Yloc_k_filt(:,k) + s(2), celldata(ic).Xloc_k_filt(:,k) + s(1)];
        disp = [celldata(ic).Ydisp_k(:,k), celldata(ic).Xdisp_k(:,k)];
        disp_filt = [celldata(ic).Ydisp_k_filt(:,k), celldata(ic).Xdisp_k_filt(:,k)];
        forc = [celldata(ic).Yforce_k(:,k), celldata(ic).Xforce_k(:,k)];
        trac = [celldata(ic).Ytrac_k(:,k), celldata(ic).Xtrac_k(:,k)];
    
        if exist('img_BD','var')
            set(img_video,'CData',img_BD_crop(:,:,k))
        end
        set(p_pos,'XData',loc_filt(:,2),'YData',loc_filt(:,1))
%         set(q_disp,'XData',loc_filt(:,2),'YData',loc_filt(:,1),'UData',disp_filt(:,2),'VData',disp_filt(:,1))
        set(q_forc,'XData',loc_filt(:,2),'YData',loc_filt(:,1),'UData',arrowscale*forc(:,2),'VData',arrowscale*forc(:,1))
        set(tl,'String',sprintf('%d',k))
        drawnow
    %     pause(0.1)
        k = k + 1;
    end
end

%% polar histogram
% ic = 1;
% bin_edges = linspace(0,2*pi,9);
% force_rose_sum = zeros(length(bin_edges)-1,meta_BD.nFrames);
% [force_angles,force_mag] = cart2pol(celldata(ic).Xforce_k,celldata(ic).Yforce_k);
% fangle_in_cell = force_angles(celldata(ic).celldots(:),:) + pi;
% fmag_in_cell = force_mag(celldata(ic).celldots(:),:);
% for k = 1:meta_BD.nFrames
%     for ibin = 1:length(bin_edges)-1
%         fmag_in_bin = fmag_in_cell(fangle_in_cell(:,k) > bin_edges(ibin) & fangle_in_cell(:,k) < bin_edges(ibin+1),k);
%         fangle_in_bin = fangle_in_cell(fangle_in_cell(:,k) > bin_edges(ibin) & fangle_in_cell(:,k) < bin_edges(ibin+1),k);
%     %     [bin_edges(ibin)*180/pi bin_edges(ibin+1)*180/pi length(fmag_in_bin)]
%         force_rose_sum(ibin,k) = sum(fmag_in_bin,1);
%     end
% end
% 
% sc = 1e-3; % scale factor for arrows
% figure
% pax = polaraxes;
% for k = 1:meta_BD.nFrames
% % for k = meta_BD.uFrame
%     polarhistogram(pax,'BinEdges',bin_edges,'BinCounts',1e-3*force_rose_sum(end:-1:1,k))
%     pax.RLim = [0 1e-3*max(force_rose_sum(:))];
%     pax.RLimMode = 'manual';
%     hold on
% %     [x,y] = pol2cart(fangle_in_bin,fmag_in_bin);
% %     quiver(zeros(size(x)),zeros(size(y)),x,y,0,'-k')
%     polarplot(-[fangle_in_cell(:,k),fangle_in_cell(:,k)]', sc*[zeros(size(fmag_in_cell(:,k))),fmag_in_cell(:,k)]','-k')
%     hold off
%     title(sprintf('Frame #: %i',k));
%     drawnow
% end

%% force transient over time
for ic = 1:nCells
    if ~isempty(celldata(ic))
        figure('Position',[700,100,600,200])
        plot(meta_BD.Time,10^-3*celldata(ic).total_force,'-k','linewidth',2)
        xlabel('Time [s]')
        ylabel('Total Force [nN]')
        box off
        set(gca,'linewidth',1.5,'tickdir','out','XColor','k','YColor','k')
        
        print(sprintf('plot_force_transient_cell%i', ic),'-dpng')
        print(sprintf('plot_force_transient_cell%i', ic),'-dsvg','-vector')
    end
end

%% plot image with forces for publication
% ic = 1;
% % k = meta_BD.uFrame;
% % k = 100;
% k = 1;
% 
% s = celldata(ic).crop([1 2]);
% 
% pos = [celldata(ic).Yvector + s(2), celldata(ic).Xvector + s(1)];
% loc = [celldata(ic).Yloc_k(:,k) + s(2), celldata(ic).Xloc_k(:,k) + s(1)];
% disp = [celldata(ic).Ydisp_k(:,k), celldata(ic).Xdisp_k(:,k)];
% forc = [celldata(ic).Yforce_k(:,k), celldata(ic).Xforce_k(:,k)];
% trac = [celldata(ic).Ytrac_k(:,k), celldata(ic).Xtrac_k(:,k)];
% 
% fig_disp = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% ax_disp = axes(fig_disp,'units','normalized','position',[0 0 1 1]);
% % all_cell_dots = any(cat(3,celldata(ic).celldots),3);
% imagesc(img_REFBD,[min(img_REFBD(:)), max(img_REFBD(:))])
% % imagesc(img_filt)
% colormap(gray*[1 0 0;0 130/255 0;0 0 0])
% axis image
% axis manual
% axis off
% % axis([435 795 207 575])
% hold on
% p_bd = plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'-w');
% 
% color_max = [1 1 0];
% color_min = [1 0 1];
% 
% arrowscale = 0.002; % 0.002
% arrow_start = [loc(:,2),loc(:,1)];
% arrow_end = arrow_start + arrowscale*[forc(:,2), forc(:,1)];
% 
% arrowlength = sqrt(sum((arrow_end - arrow_start).^2,2))/2/2;
% arrowwidth = arrowlength*5/15;
% arrow_color = mat2gray(arrowlength).*(color_max - color_min) + color_min;
% q_trac = arrow(arrow_start,arrow_end,...
%     'BaseAngle',90,'TipAngle',25,'Width',arrowwidth,'Length',arrowlength,...
%     'LineWidth',0.25,'EdgeColor','none','FaceVertexCData',arrow_color,'FaceColor','flat'); % 'FaceColor','w','EdgeColor','none'
% 
% ax_scalebars = axes(fig_disp,'position',get(ax_disp,'position'),'hittest','off','pickableparts','none');
% 
% scalebar_length = 20/meta_BD.Calibration;
% scalebar_disp = floor(20/arrowscale)*arrowscale/meta_BD.Calibration;
% % scalebar_force = floor(20/meta_BD.Calibration/(arrowscale/0.001)/10)*(arrowscale/0.001)*10;
% % scalebar_force = scalebar_force*2/10;
% scalebar_force = 50;
% ylims = ax_disp.YLim;
% xlims = ax_disp.XLim;
% 
% axis image
% axis off
% set(ax_scalebars,'color','none',...
%     'XLim',get(ax_disp,'XLim'),...
%     'YLim',get(ax_disp,'YLim'),...
%     'XLimMode','manual','YLimMode','manual',...
%     'YDir','Reverse')
% 
% linkaxes([ax_disp, ax_scalebars]);
% 
% hold on
% patch(ax_scalebars,'XData',[xlims(2) - scalebar_length - 0.2*diff(xlims),xlims(2),xlims(2),xlims(2) - scalebar_length - 0.2*diff(xlims)],...
%     'YData',[ylims(2),ylims(2),ylims(2) - 0.15*diff(ylims),ylims(2) - 0.15*diff(ylims)],...
%     'EdgeColor','none','FaceColor','k','FaceAlpha',1)
% 
% plot(ax_scalebars,[xlims(2) - scalebar_length - 0.15*diff(xlims), xlims(2) - 0.15*diff(xlims)],...
%     [ylims(2) - 0.05*diff(ylims), ylims(2) - 0.05*diff(ylims)],'w','linewidth',5)
% text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.05*diff(ylims),...
%     sprintf('%.0f um',scalebar_length*meta_BD.Calibration),'color','w','VerticalAlignment','Middle')
% 
% arrow([xlims(2) - arrowscale*scalebar_force/0.001 - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
%     [xlims(2) - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
%     'BaseAngle',90,'TipAngle',25,'Width',max(arrowwidth),'Length',max(arrowlength),'LineWidth',0.25,'FaceColor','w','EdgeColor','none')
% % text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
% %     sprintf('%.0f nN',1/arrowscale*scalebar_force*0.001),'color','w','VerticalAlignment','Middle')
% text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
%     sprintf('%.0f nN',scalebar_force),'color','w','VerticalAlignment','Middle')
% 
% % quiver(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),celldata(ic).Xtrac_k(:,k),celldata(ic).Ytrac_k(:,k),0.5,'-c','linewidth',1)
% 
% % print('arrow_forc','-dsvg','-painters')

%% plot rotated image with forces for publication
% ic = 1;
% % k = meta_BD.uFrame;
% % k = 100;
% k = 1;
% aa = celldata(ic).rot_angle;
% s = celldata(ic).crop([1 2]);
% img_rotated = imrotate(img_REFBD,aa*180/pi,'crop');
% 
% R = [cos(aa), sin(aa); -sin(aa), cos(aa)];
% % R = [1 0; 0 1];
% centerX = floor(celldata(ic).N/2+1) + s(1);
% centerY = floor(celldata(ic).M/2+1) + s(2);
% pos = [celldata(ic).Yvector + s(2), celldata(ic).Xvector + s(1)];
% loc = [celldata(ic).Yloc_k(:,k) + s(2), celldata(ic).Xloc_k(:,k) + s(1)];
% disp = [celldata(ic).Ydisp_k(:,k), celldata(ic).Xdisp_k(:,k)];
% forc = [celldata(ic).Yforce_k(:,k), celldata(ic).Xforce_k(:,k)];
% pos_rot = (pos - [centerY, centerX])*R + [centerY, centerX];
% loc_rot = (loc - [centerY, centerX])*R + [centerY, centerX];
% disp_rot = disp*R;
% forc_rot = forc*R;
% 
% fig_disp = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% ax_disp = axes(fig_disp,'units','normalized','position',[0 0 1 1]);
% % all_cell_dots = any(cat(3,celldata(ic).celldots),3);
% imagesc(img_rotated,[min(img_REFBD(:)), max(img_REFBD(:))])
% % imagesc(img_filt)
% colormap(gray*[1 0 0;0 130/255 0;0 0 0])
% axis image
% axis manual
% axis off
% % axis([435 795 207 575])
% hold on
% p_bd = plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'--g');
% % plot(pos_rot(:,2),pos_rot(:,1),'+k','markersize',10)
% % plot(loc_rot(:,2),loc_rot(:,1),'.k','markersize',25)
% % quiver(pos_rot(:,2),pos_rot(:,1),disp_rot(:,2),disp_rot(:,1),0,'-w','linewidth',1)
% % arrow([pos_rot(:,2),pos_rot(:,1)],[pos_rot(:,2),pos_rot(:,1)]+[disp_rot(:,2),disp_rot(:,1)],...
% %     'Length',5)
% 
% arrowscale = 0.002;
% arrow_start = [loc_rot(:,2),loc_rot(:,1)];
% arrow_end = arrow_start + arrowscale*[forc_rot(:,2), forc_rot(:,1)];
% 
% arrowlength = sqrt(sum((arrow_end - arrow_start).^2,2))/2;
% arrowwidth = arrowlength*5/15;
% q_trac = arrow(arrow_start,arrow_end,...
%     'BaseAngle',90,'TipAngle',25,'Width',arrowwidth,'Length',arrowlength,...
%     'LineWidth',0.25,'FaceColor','w','EdgeColor','none');
% 
% ax_scalebars = axes(fig_disp,'position',get(ax_disp,'position'),'hittest','off','pickableparts','none');
% 
% scalebar_length = 20/meta_BD.Calibration;
% scalebar_disp = floor(20/arrowscale)*arrowscale/meta_BD.Calibration;
% % scalebar_force = floor(20/meta_BD.Calibration/(arrowscale/0.001)/10)*(arrowscale/0.001)*10;
% % scalebar_force = scalebar_force*2/10;
% scalebar_force = 20;
% ylims = ax_disp.YLim;
% xlims = ax_disp.XLim;
% 
% axis image
% axis off
% set(ax_scalebars,'color','none',...
%     'XLim',get(ax_disp,'XLim'),...
%     'YLim',get(ax_disp,'YLim'),...
%     'XLimMode','manual','YLimMode','manual',...
%     'YDir','Reverse')
% 
% linkaxes([ax_disp, ax_scalebars]);
% 
% hold on
% patch(ax_scalebars,'XData',[xlims(2) - scalebar_length - 0.2*diff(xlims),xlims(2),xlims(2),xlims(2) - scalebar_length - 0.2*diff(xlims)],...
%     'YData',[ylims(2),ylims(2),ylims(2) - 0.15*diff(ylims),ylims(2) - 0.15*diff(ylims)],...
%     'EdgeColor','none','FaceColor','k','FaceAlpha',1)
% 
% plot(ax_scalebars,[xlims(2) - scalebar_length - 0.15*diff(xlims), xlims(2) - 0.15*diff(xlims)],...
%     [ylims(2) - 0.05*diff(ylims), ylims(2) - 0.05*diff(ylims)],'w','linewidth',5)
% text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.05*diff(ylims),...
%     sprintf('%.0f um',scalebar_length*meta_BD.Calibration),'color','w','VerticalAlignment','Middle')
% 
% arrow([xlims(2) - arrowscale*scalebar_force/0.001 - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
%     [xlims(2) - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
%     'BaseAngle',90,'TipAngle',25,'Width',max(arrowwidth),'Length',max(arrowlength),'LineWidth',0.25,'FaceColor','w','EdgeColor','none')
% % text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
% %     sprintf('%.0f nN',1/arrowscale*scalebar_force*0.001),'color','w','VerticalAlignment','Middle')
% text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
%     sprintf('%.0f nN',scalebar_force),'color','w','VerticalAlignment','Middle')
% 
% % quiver(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),celldata(ic).Xtrac_k(:,k),celldata(ic).Ytrac_k(:,k),0.5,'-c','linewidth',1)
% 
% % print('arrow_forc','-dsvg','-painters')