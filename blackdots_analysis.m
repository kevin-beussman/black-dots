%% gridded TFM (black dots) analysis
% Kevin Beussman / Andrea Leonard / Molly Mollica / Ziz
% % beussk@uw.edu
clearvars; clc; close all;

addpath(genpath('G:\Shared drives\CBL\Kevin_Beussman\Analysis Codes\tools'))

uM.force_user_vals = false; % set this to 'true' to if you want to use the below values instead of the ND2 values
uM.track_method = 1; % method 1 = more manual method, method 2 = more automated method

% these following settings are used only if they cannot be pulled from the video data
uM.CameraPixelSize = 3.63; % physical pixel size [um], from camera manufacturer V3=6.5, Flash2.8=3.63, Andor=6.45
uM.Objective = 60; % e.g. 60x objective
uM.Binning = 1; % e.g. 1x1 binning
uM.FrameRate = 1; % [Hz]
uM.CouplerRatio = 1; % e.g. 0.7x coupler
uM.MagMultiplier = 1; % e.g. 1.5x magnification toggle

% DotSize/DotSpacing
% N1: 1.2/2.0 or 1.0/1.95 or 0.8/1.95
% 2E: 3.0/6.0
% 2A: 3.5/8.0

% PDMS properties:
% 0% = 5 kPa
% 5% = 13.5 kPa
% 10% = 47.8 kPa

uM.DotSize = 0.8; % [um] microns
uM.DotSpacing = 1.95; % [um] microns
uM.YoungsModulus = 13500; % [N*m^-2] or equivalently [pN*um^-2]
uM.Poisson = 0.5; % [] unitless

% uM.Frames = 1; % frames to analyze. comment this out to use all frames
uM.BDchannel = 4; % channel to use for black dots. comment this out if there is only one channel
uM.CBchannel = 1; % channel to use for cell boundary. comment this out if there is only one channel
% uM.FAchannel = 2; % channel to use for focal adhesions. comment this out if there is only one channel
uM.uFrame = 1; % relaxed frame -- this is relative to uM.Frames
uM.Crop = false; % true or [xmin ymin, width height]
uM.crop_around_cells = true;
uM.manual_boundary = true;
uM.uPoints = 2; % number of outer undeformed points in each direction to use for grid estimation. comment out to use all points
uM.useLcurve = false; % calculate regularization parameter from L curve (takes a while)
uM.regParam = 5e-8; % regularization parameter guess

%% Choose file to analyze
[file_raw, path_raw] = uigetfile({'*.nd2','NIS Elements';'*.cxd','HCImage Live';'*.avi','(AVI) Audio Video Interleave';'*.tif;*.tiff','TIFF Image Stacks'},'Select file to analyze. NOTE: .nd2 files require bioformats plugin.','MultiSelect','off');
% path_raw = 'E:\k\2019-05-30 BDs CMs Live\10pct_20k_1\';
% files_raw = dir;
% files_raw = {files_raw([files_raw.isdir]).name};
% files_raw = files_raw(contains(files_raw,'_blackdots_analysis'));
% name_ind = regexp(files_raw,'_blackdots_analysis');
% for f = 1:length(files_raw)
%     files_raw{f} = [files_raw{f}(1:name_ind{f}-1) '.nd2'];
% end

cd(path_raw)

[vidFrames,vM] = cbl_vidread(file_raw,uM);
% vidFrames = imcomplement(vidFrames); % use this if dots are bright
vidFrames_raw = vidFrames;
vMraw = vM;

% %%%%%%%
% pre_rotate = 0; % -23.448;
% % pre_crop = [502 451 1089 1089];
% pre_crop = [638 544 1135 1039];
% vidFrames = [];
% for k = 1:vM.nFrames
%     vidFrames(:,:,k) = imcrop(imrotate(vidFrames_raw(:,:,k),pre_rotate,'bilinear'),pre_crop);
% end
% vidFrames_raw = vidFrames;
% [vM.M, vM.N] = size(vidFrames(:,:,1));
% vM.Crop = pre_crop;
% %%%%%%% KB2 2020/05/27

[~,filename,~] = fileparts(file_raw);

% create save folder(s)
cur_time = char(datetime('now','format','yyyyMMdd_HHmmss'));

path_save = [path_raw filename '_blackdots_analysis' filesep];
if ~exist(path_save,'dir')
    mkdir(path_save)
end

path_save = [path_save cur_time filesep];
if ~exist(path_save,'dir')
    mkdir(path_save)
end

file_save = ['blackdots_data_' cur_time '.mat'];

%% Find cell boundary
[file_cb, path_cb] = uigetfile({'*.mat'},'Select file containing cell boundary coordinates.');
if file_cb ~= 0
    load([path_cb file_cb],'CB_uncrop','cbFrame');
    
    if vM.Crop ~= false
%         R = [cosd(-pre_rotate), sind(-pre_rotate); -sind(-pre_rotate), cosd(-pre_rotate)];
%         centerX = floor(vMraw.N/2+1);
%         centerY = floor(vMraw.M/2+1);
%         [newM, newN] = size(imrotate(cbFrame,pre_rotate));
%         centerX2 = floor(newN/2+1);
%         centerY2 = floor(newM/2+1);
%         CB_uncrop{1} = (CB_uncrop{1} - [centerX, centerY])*R + [centerX2, centerY2];
        CB = cellfun(@plus,CB_uncrop,{-vM.Crop([1 2])},'UniformOutput',false);
    else
        CB = CB_uncrop;
    end
    
     save([path_save file_save],'CB_uncrop','cbFrame')
else
    if ~isfield(vM,'CBchannel')
        CB{1} = [1,1; vM.N,1; vM.N,vM.M; 1,vM.M];
        cbFrame = vidFrames_raw(:,:,1);
        cbFrame(:,:) = 0;
    else
        [cbFrame,~] = cbl_vidread(file_raw,vM,vM.CBchannel);
        if vM.manual_boundary
            CB = get_cell_boundary(cbFrame); % manual clicking boundary
        else
            CB = get_cell_boundary2(cbFrame,vM); % semi-automatic segmentation
    %         stats = get_cell_boundary_zizhen(cbFrame,vM); % doesnt work on kb2 computer
        end
    end

    if vM.Crop ~= false
        CB_uncrop = cellfun(@plus,CB,{vM.Crop([1 2])},'UniformOutput',false);
    else
        CB_uncrop = CB;
    end
    
    save([path_save file_save],'CB_uncrop','cbFrame')
end
CBraw = CB;

%% cell selection
nCells = length(CB);
fprintf('\t%i objects (cells) detected\n',nCells)
skipped = false(nCells,1);
for ic = 1:nCells
celldata(ic).crop_tight = [min(CB{ic}), max(CB{ic})-min(CB{ic})];

expand = vM.DotSize/vM.Calibration*[-1 -1 2 2];
celldata(ic).crop = celldata(ic).crop_tight;

last_crop = celldata(ic).crop;

fig_cellselect = figure('Units','Normalized','Position',[0.2 0.1 0.6 0.6]);
ax2 = subplot(1,2,2);
im1 = imagesc(vidFrames(:,:,vM.uFrame));
axis image
hold on
plot(CB{ic}(:,1),CB{ic}(:,2),'-r')
hold off
colormap(ax2,gray*[1 0 0; 0 0.75 0; 0 0 0])

ax1 = subplot(1,2,1);
imagesc(cbFrame);
% % imagesc(imcrop(imrotate(cbFrame,pre_rotate),pre_crop))
axis image
hold on
plot(CB{ic}(:,1),CB{ic}(:,2),'-r')
hold off
colormap(ax1,gray*[0 0 0; 0 0.75 0; 0 0 1])
linkaxes([ax1 ax2])
axis(celldata(ic).crop_tight([1 1 2 2]) + [0, celldata(ic).crop_tight(3), 0, celldata(ic).crop_tight(4)] - 0.5)

but_expand = uicontrol(fig_cellselect,'Units','Normalized','Position',[0 0 0.2 0.1],'Style','PushButton','String','expand',...
    'Callback',['celldata(ic).crop = celldata(ic).crop + expand;', ...
                'if celldata(ic).crop(1) < 1, celldata(ic).crop(1) = 1; end;', ...
                'if celldata(ic).crop(2) < 1, celldata(ic).crop(2) = 1; end;', ...
                'if sum(celldata(ic).crop([2 4])) > vM.M, celldata(ic).crop(4) = vM.M - celldata(ic).crop(2) + 1; end;', ...
                'if sum(celldata(ic).crop([1 3])) > vM.N, celldata(ic).crop(3) = vM.N - celldata(ic).crop(1) + 1; end;', ...
                'axis(celldata(ic).crop([1 1 2 2]) + [0, celldata(ic).crop(3), 0, celldata(ic).crop(4)] - 0.5);', ...
                'last_crop = [last_crop; celldata(ic).crop];']);

but_undo = uicontrol(fig_cellselect,'Units','Normalized','Position',[0.2 0 0.2 0.1],'Style','PushButton','String','undo',...
    'Callback',['if size(last_crop,1) > 1, last_crop(end,:) = []; end;', ...
                'celldata(ic).crop = last_crop(end,:);', ...
                'axis(celldata(ic).crop([1 1 2 2]) + [0, celldata(ic).crop(3), 0, celldata(ic).crop(4)] - 0.5);']);

but_ok = uicontrol(fig_cellselect,'Units','Normalized','Position',[0 0.1 0.2 0.1],'Style','PushButton','String','OK','Callback','uiresume;');

but_skip = uicontrol(fig_cellselect,'Units','Normalized','Position',[0.2 0.1 0.2 0.1],'Style','PushButton','String','skip',...
    'Callback','uiresume; skipped(ic) = true;');

uiwait

close(fig_cellselect)

% if skipped(ic)
%     continue
% end

%% Bandpass filter all video frames
fprintf('Filtering image frames...')

bpass_noise = 1;
dotsize_px = 2*round((vM.DotSize/vM.Calibration + 1)/2) - 1;
vidFrames_filt = zeros(size(vidFrames));
for k = 1:vM.nFrames
    vidFrames_filt(:,:,k) = mat2gray(imcomplement(bpass_kb2(imcomplement(vidFrames(:,:,k)),bpass_noise,dotsize_px)));
end

img_raw = vidFrames(:,:,vM.uFrame);
img_filt = vidFrames_filt(:,:,vM.uFrame);

fprintf('DONE\n')

%% crop to selected cell
fprintf('Cropping around cell...')

vidFrames_crop_raw = imcrop(vidFrames_raw(:,:,1),celldata(ic).crop);
vidFrames_crop_filt = imcrop(vidFrames_filt(:,:,1),celldata(ic).crop);
for k = 2:vM.nFrames
    vidFrames_crop_raw(:,:,k) = imcrop(vidFrames_raw(:,:,k),celldata(ic).crop);
    vidFrames_crop_filt(:,:,k) = imcrop(vidFrames_filt(:,:,k),celldata(ic).crop);
end

vMcrop = vM;
vMcrop.Crop = celldata(ic).crop;
vMcrop.M = size(vidFrames_crop_raw(:,:,1),1);
vMcrop.N = size(vidFrames_crop_raw(:,:,1),2);

celldata(ic).CB = CB{ic} -vMcrop.Crop([1 2]);

fprintf('DONE\n')

%% Start analysis
if vM.track_method == 1
    %% Measure rotation angle from image
    fprintf('Measuring rotation angle...')

    vMcrop.rot_angle = get_rot_angle(vidFrames_crop_filt(:,:,vMcrop.uFrame),vMcrop);

    fprintf('DONE\n')
    
    %% Choose the sampling dots
    % allows user to estimate dot grid positions
    fprintf('Choosing dots...')

    [px_sample,py_sample] = choose_dots(vidFrames_crop_filt(:,:,vMcrop.uFrame),vMcrop);
    
    fprintf('DONE\n')
    
    %%
    if isempty(px_sample)
        skipped(ic) =  true;
        continue
    end
    
    %% Update rotation angle
    fprintf('Updating rotation angle...')

    vMcrop.rot_angle = get_rot_angle(cat(3,px_sample,py_sample),vMcrop);

    fprintf('DONE\n')
    
    %% Extend the sampling dots
    % fits a grid to the sampling dots, and extends the grid of dots
    % to fill the full image
    fprintf('Extending dots to image bounds...')

    [px_ext,py_ext] = extend_dots(px_sample,py_sample,vMcrop);

    fprintf('DONE\n')
    
    %% Find centroids
    fprintf('Finding dot centroids...')

    [px,py,real_points] = find_centroids(px_ext,py_ext,vidFrames_crop_filt(:,:,vMcrop.uFrame),vMcrop);

    fprintf('DONE\n')
    
    %STOP HERE?
    
    %% Characterize dots
    fprintf('Characterizing dots...')
    
    BD = calc_dot_size_spacing(px,py,real_points,vidFrames_crop_filt(:,:,vMcrop.uFrame),vMcrop);
    
    celldata(ic).BD = BD;
    
    fprintf('DONE\n')
    
    %% Find undeformed dot positions
    fprintf('Finding undeformed dot centroids...')

    [px,py,px0,py0] = find_undeformed(px,py,real_points,vMcrop);

%     L0 = median(BD.DotSpacings)/vM.Calibration;
%     [px2,py2,px02,py02] = find_undeformed_springs(px,py,real_points,L0,vM);
    
%     % plotting
% %     plot(px,py,'.k',px0,py0,'ok',px02,py02,'or')
% %     hold on
% %     quiver(px0,py0,(px - px0),(py - py0),'-k')
% %     quiver(px02,py02,(px2 - px02),(py2 - py02),'-r')
% %     axis image
% %     set(gca,'ydir','reverse')
    
    fprintf('DONE\n')
    
    %% Track dots across all frames
    % just basically repeats find_centroids for each frame
    fprintf('Finding dot centroids in each frame...')

    if vM.nFrames > 1
        [px_k,py_k,real_points] = track_dots_across_frames(vidFrames_crop_filt,px,py,real_points,vM);
    else
        px_k = px(:);
        py_k = py(:);
    end

    fprintf('DONE\n')
    
elseif vM.track_method == 2
    %% Find object centers
    fprintf('Finding centers in each frame...')

    dotsize = 2*round((vMcrop.DotSize/vMcrop.Calibration + 1)/2) - 1; % pixels, nearest odd integer
    dotspacing = round(vMcrop.DotSpacing/vMcrop.Calibration);

    bpass_noise = 1; % bandpass noise: 1 pixel is typical      % 1*vMcrop.Calibration
    pkfnd_threshold = 0.15; % object peak threshold: find by trial and error, 0.15 is good

    centers = []; % centers has 3 columns: [x, y, frame#]. we don't know how many objects there will be
    pct = 0;
    fprintf('\n[%-20s] %3.0f%%\n',repmat('|',1,round(pct/100*20)),pct)
    for k = 1:vMcrop.nFrames
        pct = k/vMcrop.nFrames;
        fprintf([repmat('\b',1,28) '[%-20s] %3.0f%%\n'],repmat('|',1,round(pct*20)),pct*100)
        I1 = vidFrames_crop_raw(:,:,k);
        I2 = bpass_kb2(imcomplement(I1),bpass_noise,dotsize); % spatial bandpass filter
        pk = pkfnd_kb2(mat2gray(I2),pkfnd_threshold,dotsize); % find objects whose intensity is above threshold
        dsz = ceil(dotsize*1.2); if ~mod(dsz,2), dsz = dsz + 1; end % KB2 edit 9/3/2020
        cnt = cntrd_kb2(I2,pk,dsz); % KB2 edit 9/3/2020
        
%         if k == 125;
%             1;
%         elseif k == 150;
%             1;
%         end

        centers = [centers; cnt(:,1:2), ones(size(cnt,1),2)*k];
    end

%     figure(1)
%     imagesc(vidFrames(:,:,1))
%     hold on
%     plot(centers(:,1),centers(:,2),'.r','markersize',10)
%     hold off
    
    fprintf('DONE\n')
    
    %% Track objects across video frames
    fprintf('Tracking centers...')

    if ~exist('track','file')
        error('''track'' not found, visit http://site.physics.georgetown.edu/matlab/')
    end
    if vMcrop.nFrames >= 2

        param.mem = 5; % mem = # of "dropped" frames before calling it a new object
        param.good = 0; % good = ?
        param.dim = 2; % dim = dimensions, (x, y)
        param.quiet = 0; % quiet = no command line feedback

        res = track(centers,floor(dotspacing/4),param);
        res = sortrows(res,5);% res has 5 columns: [x, y, frame#, frame#, object#]
    else
        res = [centers(:,1:2), ones(size(centers,1),2), (1:size(centers,1))'];
    end

    fprintf('DONE\n')

    %% Find and analyze tracked points
    % "close points" are points that appear on at least 75% of frames
    fprintf('Filling in temporal gaps in position data...')

    if vMcrop.nFrames >= 2
        full_points = NaN(max(res(:,5)),1);
        close_points = NaN(max(res(:,5)),1);
        p1x = NaN(max(res(:,5)),vMcrop.nFrames);
        p1y = NaN(max(res(:,5)),vMcrop.nFrames);
        p2x = NaN(max(res(:,5)),vMcrop.nFrames);
        p2y = NaN(max(res(:,5)),vMcrop.nFrames);

        ifp = 0; icp = 0;
        pct = 0;
        point_markers = [0; find(diff(res(:,5))); size(res,1)];
        max_npt = max(res(:,5));
        fprintf('\n[%-20s] %3.0f%%\n',repmat('|',1,round(pct*20)),pct*100)
    for npt = 1:max_npt
            pct = npt/max_npt;
            fprintf([repmat('\b',1,28) '[%-20s] %3.0f%%\n'],repmat('|',1,round(pct*20)),pct*100)      
            i_1 = point_markers(npt)+1; % this is much faster!
            i_2 = point_markers(npt+1);
            if nnz(res(i_1:i_2,5) == npt) == vMcrop.nFrames
                ifp = ifp + 1;
                full_points(ifp) =  npt;
                p1x(ifp,:) = res(i_1:i_2,1)';
                p1y(ifp,:) = res(i_1:i_2,2)';
            elseif nnz(res(i_1:i_2,5) == npt) >= round(0.75*vMcrop.nFrames)
                icp = icp + 1;
                close_points(icp) = npt;

                inds = find((res(:,5) == close_points(icp)));
                inds2 = res(inds,4);

                tempx = NaN(1,vMcrop.nFrames);
                tempx(inds2) = res(inds,1)';
                tempy = NaN(1,vMcrop.nFrames);
                tempy(inds2) = res(inds,2)';
                p2x(icp,:) = tempx;
                p2y(icp,:) = tempy;
            end
    end 

        full_points = full_points(~isnan(full_points));
        close_points = close_points(~isnan(close_points));
        p1x = p1x(1:ifp,:);
        p1y = p1y(1:ifp,:);
        p2x = p2x(1:icp,:);
        p2y = p2y(1:icp,:);
    else
        p1x = res(:,1);
        p1y = res(:,2);
        p2x = [];
        p2y = [];
    end

    if vMcrop.nFrames >= 2
        p2x_f = p2x;
        p2y_f = p2y;
        for icp = 1:length(close_points)
            for k = find(isnan(p2x(icp,:)))
                cL = (k - 1)*(k > 1) + 1*(k == 1);
                cR = (k + 1)*(k < vMcrop.nFrames) + vMcrop.nFrames*(k == vMcrop.nFrames);
                while isnan(p2x(icp,cL)) && isnan(p2x(icp,cR))
                    cL = (cL - 1)*(cL > 1) + 1*(cL == 1);
                    cR = (cR + 1)*(cR < vMcrop.nFrames) + vMcrop.nFrames*(cR == vMcrop.nFrames);
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
    else
        px = p1x;
        py = p1y;
    end

    fprintf('DONE\n')
    
    %% manual point removal (and addition, addition not functional yet)
%     fig_1 = figure;
%     ax_1 = axes;
%     im_1 = imagesc(vidFrames(:,:,vMcrop.uFrame));
%     hold on
%     p_dot = plot(px(:,vMcrop.uFrame),py(:,vMcrop.uFrame),'.r','markersize',10); % plot x y positions found above
%     p2_dot = plot(px(1,vMcrop.uFrame),py(1,vMcrop.uFrame),'or'); 
%     p_dot_remove = plot(0,0,'.k','markersize',10);
%     hold off
% 
%     set(im_1,'hittest','off')
%     set(p_dot,'hittest','off')
%     set(p2_dot,'hittest','off')
% 
%     fcn2_1 = @(a,b) set(fig_1,'UserData',find(sqrt((px(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,2)).^2) == min(sqrt((px(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,2)).^2))));
%     fcn2_5 = @(a,b) set(p2_dot,'xdata',px(fig_1.UserData,vMcrop.uFrame),'ydata',py(fig_1.UserData,vMcrop.uFrame));
%     fcn2 = @(a,b) cellfun(@feval,{fcn2_1 fcn2_5});
%     set(ax_1,'buttondownfcn',fcn2)
% 
%     p_remove = false(size(px));
%     but_remove = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
%         'Position',[0.8 0 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
%         'String','REMOVE','Callback','p_remove(fig_1.UserData) = true; set(p_dot_remove,''xdata'',px(p_remove,vMcrop.uFrame),''ydata'',py(p_remove,vMcrop.uFrame));');
% 
%     % but_add = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
%     %     'Position',[0.1 0 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
%     %     'String','ADD','Callback','p_remove(fig_1.UserData) = true; set(p_dot_remove,''xdata'',px(p_remove,vMcrop.uFrame),''ydata'',py(p_remove,vMcrop.uFrame));');
% 
%     but_done = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
%         'Position',[0.8 0.95 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
%         'String','DONE','Callback','uiresume');
% 
%     uiwait
%     
%     close(fig_1)
% 
%     px(p_remove,:) = [];
%     py(p_remove,:) = [];

    fig_1 = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
    ax_1 = axes;
    im_1 = imagesc(imadjust(mat2gray(vidFrames_crop_raw(:,:,vMcrop.uFrame))));
    axis image
%     colormap(gray)
    hold on
    p_dot = plot(px(:,vMcrop.uFrame),py(:,vMcrop.uFrame),'.r','markersize',10);
    p2_dot = plot(px(1,vMcrop.uFrame),py(1,vMcrop.uFrame),'or');
    p_dot_remove = plot(0,0,'.k','markersize',10);
    hold off

    set(im_1,'hittest','off')
    set(p_dot,'hittest','off')
    set(p2_dot,'hittest','off')

    fcn1_1 = @(a,b) set(fig_1,'UserData',find(sqrt((px(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,2)).^2) == min(sqrt((px(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,2)).^2))));
    fcn1_2 = @(a,b) set(p2_dot,'xdata',px(fig_1.UserData,vMcrop.uFrame),'ydata',py(fig_1.UserData,vMcrop.uFrame));
    fcn1 = @(a,b) cellfun(@feval,{fcn1_1 fcn1_2});

    fcn2_1 = @(a,b) set(fig_1,'UserData',[fig_1.UserData; find(sqrt((px(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,2)).^2) == min(sqrt((px(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,vMcrop.uFrame) - ax_1.CurrentPoint(1,2)).^2)))]);
    fcn2_2 = @(a,b) set(p2_dot,'xdata',px(fig_1.UserData,vMcrop.uFrame),'ydata',py(fig_1.UserData,vMcrop.uFrame));
    fcn2 = @(a,b) cellfun(@feval,{fcn2_1 fcn2_2});

    fcn3 = @(a,b) set(ax_1,'buttondownfcn',fcn2);
    fcn4 = @(a,b) set(ax_1,'buttondownfcn',fcn1);

    set(ax_1,'buttondownfcn',fcn1)
    set(fig_1,'KeyPressFcn',fcn3,'KeyReleaseFcn',fcn4)

    p_remove = false(size(px));
    but_remove = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
        'Position',[0.8 0 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
        'String','REMOVE','Callback','p_remove(fig_1.UserData) = true; set(p_dot_remove,''xdata'',px(p_remove,vMcrop.uFrame),''ydata'',py(p_remove,vMcrop.uFrame));');

    but_done = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
        'Position',[0.8 0.95 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
        'String','DONE','Callback','uiresume');

    uiwait

    px(p_remove,:) = [];
    py(p_remove,:) = [];

    close(fig_1)
    
    %% get rotation angle based on dot locations
    fprintf('Finding rotation angle...')

    vMcrop.rot_angle = get_rot_angle([px(:,vMcrop.uFrame),py(:,vMcrop.uFrame)],vMcrop);

    fprintf('DONE\n')
    
    %% find undeformed state of grid
    fprintf('Finding undeformed state of grid...')

    [px_k,py_k,px0,py0,real_points] = find_undeformed_grid(px,py,vMcrop); % IF YOU HAVE A LOT OF DISPLACEMENT, CHANGE INF TO 2 to find fewer points
    px = px_k(:,:,vMcrop.uFrame);
    py = py_k(:,:,vMcrop.uFrame);
    
    fprintf('DONE\n')
    
    %% calculate dot size and spacing from image
    fprintf('Characterizing dots...')
    
    BD = calc_dot_size_spacing(px,py,real_points,vidFrames_crop_filt(:,:,vMcrop.uFrame),vMcrop);
    
    celldata(ic).BD = BD;
    
    fprintf('DONE\n')
    
end

%% Filter locations
fprintf('Filtering locations...')

Xloc_k = reshape(px_k,[],vMcrop.nFrames);
Xvector = reshape(px0,[],1);

Yloc_k = reshape(py_k,[],vMcrop.nFrames);
Yvector = reshape(py0,[],1);

% filter
uHz = 4; % lowpass cutoff
% lHz = 0.5;
if vMcrop.nFrames >= 2
%     [B,A] = butter(3,[lHz/(vMcrop.FrameRate/2), uHz/(vMcrop.FrameRate/2)],'bandpass');
    [B,A] = butter(3,uHz/(vMcrop.FrameRate/2),'low');

    Xloc_k_filt = NaN(size(Xloc_k));
    Yloc_k_filt = Xloc_k_filt;
    for pt = 1:size(Xloc_k,1)
        Xloc_k_temp = [Xloc_k(pt,end:-1:1), Xloc_k(pt,:), Xloc_k(pt,end:-1:1)];
        Yloc_k_temp = [Yloc_k(pt,end:-1:1), Yloc_k(pt,:), Yloc_k(pt,end:-1:1)];
        Xloc_k_temp_filt = filtfilt(B,A,Xloc_k_temp);
        Yloc_k_temp_filt = filtfilt(B,A,Yloc_k_temp);

        Xloc_k_filt(pt,:) = Xloc_k_temp_filt(size(Xloc_k,2)+1:2*size(Xloc_k,2));
        Yloc_k_filt(pt,:) = Yloc_k_temp_filt(size(Xloc_k,2)+1:2*size(Xloc_k,2));
    end
else
    Xloc_k_filt = Xloc_k;
    Yloc_k_filt = Yloc_k;
end
fprintf('DONE\n')

%% Calculate displacements
fprintf('Calculating displacements...')

Xdisp_k = Xloc_k - Xvector; % pixels
Xdisp_k_filt = Xloc_k_filt - Xvector;
Ydisp_k = Yloc_k - Yvector;
Ydisp_k_filt = Yloc_k_filt - Yvector;

Xgrid = px0;
Ygrid = py0;

% XYdisp = [Xdisp_k_filt(real_points(:),1), Ydisp_k_filt(real_points(:),1)]';
% for f = 1:size(XYdisp,2)
%     totaldisp(f) = norm(XYdisp(:,f))*vMcrop.Calibration;
% end

fprintf('DONE\n')

%% find cell points

% CB = CBraw;
% nCells = length(CBraw);
% fprintf('\t%i cells detected\n',nCells)
% edge_tooclose = 2*mean(BD.DotSpacings)/vM.Calibration;
dot_consider = 0.75*mean(BD.DotSpacings)/vMcrop.Calibration;

tb = celldata(ic).CB; % temporary cell boundary data
%     fprintf('\tCell %i...',nc2)

%     %check that cell is not too close to outside
%     if any([(min(tb(:,1)) - 1), (min(tb(:,2)) - 1), (vM.N - max(tb(:,1))), (vM.M - max(tb(:,2)))] < edge_tooclose)
%        fprintf('TOO CLOSE TO EDGE\n')
%        CB(ic+1) = [];
%        continue
%     else
%         fprintf('OK\n')
%         ic = ic + 1;
%     end
    
dot_dist_from_CB = zeros(size(px));
for iy = 1:size(px,1)
    for jx = 1:size(px,2)
        dot_dist_from_CB(iy,jx) = min(sqrt((px(iy,jx) - tb(:,1)).^2 + (py(iy,jx) - tb(:,2)).^2));
    end
end

celldots = inpolygon(px,py,tb(:,1),tb(:,2)) | (dot_dist_from_CB < dot_consider);

%% threshold displacements
% (if displacement < threshold, set = 0)
fprintf('Thresholding displacements...')

Xdisp_all = Xdisp_k_filt;
Ydisp_all = Ydisp_k_filt;
% Xdisp_out = Xdisp_k_filt(real_points(:) & ~celldots(:));
% Ydisp_out = Ydisp_k_filt(real_points(:) & ~celldots(:));

disp_all = sqrt(Xdisp_all.^2 + Ydisp_all.^2)*vMcrop.Calibration;
% disp_out = sqrt(Xdisp_out.^2 + Ydisp_out.^2)*vMcrop.Calibration;

% dist_out_boot = bootstrp(2000,@median,disp_out);
% noise_limit = quantile(disp_out,0.95); % this is the level that includes the lower 95% of displacements outside the cell;
% 
% noise_limit = 0.16867; % [um] for 2E substrate 0% stiffness, this is peak difference
% noise_limit = 0.22054; % [um] for 2E substrate 0% stiffness, this is 95% quantile
% noise_limit = 0.21856; % [um] for 2E substrate 5% stiffness, this is 95% quantile
% noise_limit = 0.22054; % [um] for 2E substrate 0% stiffness, this is 95% quantile
noise_limit = 0;
noise_points = reshape(all(disp_all < noise_limit,2),size(real_points));

Xdisp_k_dn = Xdisp_k_filt;
Ydisp_k_dn = Ydisp_k_filt;
Xdisp_k_dn(noise_points(:),:) = 0;
Ydisp_k_dn(noise_points(:),:) = 0;

fprintf('DONE\n')

%% calculate traction forces
fprintf('Calculating tractions...')

[n_row,n_col] = size(Xgrid);

vMcrop.regParam = vM.regParam;

E = vMcrop.YoungsModulus;
nu = vMcrop.Poisson;
d = mean(BD.DotSpacings);

Xtrac_k = [];
Ytrac_k = [];
Xforce_k = [];
Yforce_k = [];

ind_k = [(vMcrop.uFrame:-1:1) (vMcrop.uFrame+1:vMcrop.nFrames)];

% FTTC with Regularization
% for k = 1:vMcrop.nFrames
for kk = 1:length(ind_k)
    k = ind_k(kk);
    u_x = reshape(Xdisp_k_dn(:,k),n_row,n_col)*vMcrop.Calibration;
    u_y = reshape(Ydisp_k_dn(:,k),n_row,n_col)*vMcrop.Calibration;
    u = cat(3,u_x,u_y);
    if kk == 1
        if vMcrop.useLcurve
            L = logspace(-3,-8,100);
%             [~,resid_norm(k,:),soln_norm(k,:)] = calcforce_regFTTC(u,E,nu,d,L,real_points);
%             [reg_optimal(k),i_reg_optimal(k)] = find_L(resid_norm(k,:),soln_norm(k,:),L,vMcrop.regParam,'optimal');
%             [reg_corner(k),i_reg_corner(k)] = find_L(resid_norm(k,:),soln_norm(k,:),L,vMcrop.regParam,'corner');
            [~,resid_norm,soln_norm] = calcforce_regFTTC(u,E,nu,d,L,real_points);
            [reg_optimal,i_reg_optimal] = find_L(resid_norm,soln_norm,L,vMcrop.regParam,'optimal');
            [reg_corner,i_reg_corner] = find_L(resid_norm,soln_norm,L,vMcrop.regParam,'corner');
%             L = reg_optimal(k);
            L = vMcrop.regParam;
            
%             loglog(resid_norm,soln_norm,'-k',resid_norm(i_reg_optimal),soln_norm(i_reg_optimal),'or',resid_norm(i_reg_corner),soln_norm(i_reg_corner),'ob')
        else
            L = vMcrop.regParam;
        end
    end

    trac = calcforce_regFTTC(u,E,nu,d,L,real_points);
    % trac has same units as E, young's modulus
    % [N*m^-2] is equivalent to [pN*um^-2]

    Xtrac_k(:,k) = reshape(trac(:,:,1),[],1); % [pN*um^-2] = [uN*mm^-2] = [N*m^-2]
    Ytrac_k(:,k) = reshape(trac(:,:,2),[],1); % [pN*um^-2]
    Xforce_k(:,k) = Xtrac_k(:,k)*d^2; % pN
    Yforce_k(:,k) = Ytrac_k(:,k)*d^2; % pN
end

% force_mag_k = sqrt(Xforce_k.^2 + Yforce_k.^2);
% total_force_k = sum(force_mag_k,1);
    
%     F_mag_cell = sqrt(Xforce_k(cell_points(:),:).^2 + Yforce_k(cell_points(:),:).^2);
%     D_mag_cell = sqrt(Xdisp_k_filt(cell_points(:),:).^2 + Ydisp_k_filt(cell_points(:),:).^2);

% TRPF
% % % if ~isempty(FA)
% % %     pos_u = [Xvector, Yvector]*vMcrop.Calibration;
% % %     pos_f = FA.Centroids*vMcrop.Calibration;
% % % 
% % %     if vMcrop.useLcurve
% % %         if ~exist('regParamSelecetionLcurve','file')
% % %             addpath('G:\Shared Drives\CBL\Kevin_Beussman\Analysis Codes\External Codes\TFM_Version_1.11\TFM Version 1.11')
% % %         end
% % % 
% % %         u = cat(3,Xdisp_k_filt(:,vMcrop.uFrame),Ydisp_k_filt(:,vMcrop.uFrame))*vMcrop.Calibration;
% % %         L = logspace(-3,-8,100);
% % %         [~,resid_norm,soln_norm] = calcforce_regTRPF(pos_u,pos_f,u,E,nu,d,L);
% % %         [reg_corner,i_reg] = find_L(resid_norm',soln_norm',L,5e-6);
% % %         L = reg_corner;
% % %     else
% % %         reg_corner = [];
% % %         L = vMcrop.regParam;
% % %     end
% % % 
% % %     % TRPF with Regularization
% % %     u = cat(3,Xdisp_k_filt,Ydisp_k_filt)*vMcrop.Calibration;
% % %     force_TRPF = calcforce_regTRPF(pos_u,pos_f,u,vMcrop.YoungsModulus,vMcrop.Poisson,mean(BD.DotSpacings),L);
% % %     XforceTRPF_k = force_TRPF(:,:,1); % piconewtons
% % %     YforceTRPF_k = force_TRPF(:,:,2); % piconewtons
% % % else
% % %     reg_corner = [];
% % %     L = vMcrop.regParam;
% % % end

fprintf('DONE\n')

%% calculate data for each cell
fprintf('Calculating data for each cell...\n')

celldata(ic).vM = vMcrop;

celldata(ic).px0 = px0;
celldata(ic).py0 = py0;
celldata(ic).px = px;
celldata(ic).py = py;
celldata(ic).px_k = px_k;
celldata(ic).py_k = py_k;

celldata(ic).Xloc_k = Xloc_k;
celldata(ic).Yloc_k = Yloc_k;
celldata(ic).Xloc_k_filt = Xloc_k_filt;
celldata(ic).Yloc_k_filt = Yloc_k_filt;
celldata(ic).Xvector = Xvector;
celldata(ic).Yvector = Yvector;
celldata(ic).Xgrid = Xgrid;
celldata(ic).Ygrid = Ygrid;

celldata(ic).Xdisp_k = Xdisp_k;
celldata(ic).Ydisp_k = Ydisp_k;
celldata(ic).Xdisp_k_filt = Xdisp_k_filt;
celldata(ic).Ydisp_k_filt = Ydisp_k_filt;

celldata(ic).Xdisp_k_dn = Xdisp_k_dn;
celldata(ic).Ydisp_k_dn = Ydisp_k_dn;

celldata(ic).noise_limit = noise_limit;
celldata(ic).noise_points = noise_points;
celldata(ic).real_points = real_points;

celldata(ic).Xtrac_k = Xtrac_k;% pN*um^-2
celldata(ic).Ytrac_k = Ytrac_k;% pN*um^-2
celldata(ic).Xforce_k = Xforce_k; % pN
celldata(ic).Yforce_k = Yforce_k; % pN

celldata(ic).reg_optimal = reg_optimal;
celldata(ic).reg_corner = reg_corner;
celldata(ic).resid_norm = resid_norm;
celldata(ic).soln_norm = soln_norm;

celldata(ic).celldots = celldots;

% celldata(ic).total_trac = sum(sqrt(Xtrac_k(celldots(:),:).^2 + Ytrac_k(celldots(:),:).^2)); % pN*um^-2

celldata(ic).total_force = sum(sqrt(Xforce_k(celldots(:),:).^2 + Yforce_k(celldots(:),:).^2)); % pN
celldata(ic).net_force = sqrt(sum(Xforce_k(celldots(:),:)).^2 + sum(Yforce_k(celldots(:),:)).^2); % pN

celldata(ic).total_disp = sum(sqrt(Xdisp_k_dn(celldots(:),:).^2 + Ydisp_k_dn(celldots(:),:).^2)); % um
celldata(ic).net_disp = sqrt(sum(Xdisp_k_dn(celldots(:),:)).^2 + sum(Ydisp_k_dn(celldots(:),:)).^2); % um

celldata(ic).nDots = nnz(celldots);
celldata(ic).area = polyarea(tb(:,1),tb(:,2))*vMcrop.Calibration.^2; % [um^2]

% figure
% plot(px,py,'.k')
% set(gca,'ydir','reverse')
% hold on
% plot(CB{1}(:,1),CB{1}(:,2),'-r')
% plot(px(cell_dots{nc}),py(cell_dots{nc}),'ok')
% hold off
fprintf('DONE\n')

end

%% save the data
fprintf('Saving...')

% cur_time = char(datetime('now','format','yyyyMMdd_HHmmss'));

% file_save = 'blackdots_data.mat';
% file_save = ['blackdots_data_' cur_time '.mat'];

save([path_save file_save],...
    'file_raw','path_raw','uM','vM',...
    'img_raw','img_filt',...
    'celldata',...
    '-append')
% 'Xloc_k','Yloc_k','Xloc_k_filt','Yloc_k_filt',...
%     'Xvector','Yvector','Xgrid','Ygrid',...
%     'Xdisp_k','Ydisp_k','Xdisp_k_filt','Ydisp_k_filt',...
%     'Xtrac_k','Ytrac_k','Xforce_k','Yforce_k',...
%     'px0','py0','px','py','px_k','py_k',...
%     'real_points','reg_optimal','reg_corner',...
%     'BD','CB',

%     'BD','FA','px_sample','py_sample',...
%     'cellboundary','file_cb','path_cb','vM_cb',...
%     'tform_data','file_stain','path_stain','vM_stain')

% path_backup = [path_save 'backup\'];
% if ~exist(path_backup,'dir')
%     mkdir(path_backup)
% end
% 
% file_backup = ['blackdots_data_' cur_time '.mat'];
% 
% copyfile([path_save file_save],[path_backup file_backup])

fprintf('DONE\n')

error('We''re done here! Continue code to plot')

%% example plot
fig_ex = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
% all_cell_dots = any(cat(3,celldata(ic).celldots),3);
imagesc(img_raw)
% imagesc(img_filt)
colormap(gray*[1 0 0;0 130/255 0;0 0 0])
axis image
axis manual
axis off
hold on
nCells = length(celldata);
for k = 1%:vM.uFrame
    for ic = 1%:nCells
        if ~isempty(celldata(ic))
            s = celldata(ic).vM.Crop([1 2]) - [1, 1];
            p_bd = plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'-c','linewidth',2);
            plot(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),'.w','markersize',8)
            plot(celldata(ic).px_k(:,:,k) + s(1),celldata(ic).py_k(:,:,k) + s(2),'.w','markersize',8)
            plot(celldata(ic).Xvector(celldata(ic).celldots) + s(1),celldata(ic).Yvector(celldata(ic).celldots) + s(2),'ow')
%             plot(celldata(ic).Xloc_k(:,k) + s(1),celldata(ic).Yloc_k(:,k) + s(2),'.w','markersize',8)
%             quiver(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),celldata(ic).Xdisp_k(:,k),celldata(ic).Ydisp_k(:,k),1,'-w','linewidth',1)
            quiver(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),celldata(ic).Xtrac_k(:,k),celldata(ic).Ytrac_k(:,k),0.5,'-c','linewidth',1)

%             for nx = 1:size(celldata(ic).Xgrid,1)
%                 plot(celldata(ic).Xgrid(nx,:) + s(1),celldata(ic).Ygrid(nx,:) + s(2),'-w')
%             end
%             for ny = 1:size(celldata(ic).Xgrid,2)
%                 plot(celldata(ic).Xgrid(:,ny) + s(1),celldata(ic).Ygrid(:,ny) + s(2),'-w')
%             end
        end
    end
end

% plot(px0,py0,'oy')
% % plot(Xloc_k(:,1),Yloc_k(:,1),'.w','markersize',8)
% % quiver(Xloc_k(:,1),Yloc_k(:,1),Xdisp_k(:,1),Ydisp_k(:,1),1,'-w','linewidth',1)
% % quiver(Xloc_k(:,1),Yloc_k(:,1),Xtrac_k(:,1),Ytrac_k(:,1),0.5,'-c','linewidth',1)

%% traction heatmap
ic = 1;
% k = 123;
k = 1;
ninterp = 1;

trac_mag = sqrt(celldata(ic).Xtrac_k(:,k).^2 + celldata(ic).Ytrac_k(:,k).^2);
[tx,ty] = meshgrid(linspace(0,celldata(ic).vM.N,size(celldata(ic).Xgrid,2)*ninterp),linspace(0,celldata(ic).vM.M,size(celldata(ic).Xgrid,1)*ninterp));
trac_f = scatteredInterpolant(celldata(ic).Xvector,celldata(ic).Yvector,trac_mag,'natural');
trac_2 = trac_f(tx,ty);

fig_ex = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
imagesc(tx(:),ty(:),trac_2)
colormap(jet)
axis image
axis manual
axis off
hold on
% if ~isempty(celldata(ic).CB)
%     nCells = length(celldata(ic).CB);
%     for ic = 1:nCells
%         p_bd = plot(celldata(ic).CB(:,1),celldata(ic).CB(:,2),'-w','linewidth',2);
%     end
% end
p_bd = plot(celldata(ic).CB(:,1),celldata(ic).CB(:,2),'-w','linewidth',2);

%% interactive plot -- calculate specific forces (not working yet)
% fig_1 = figure;
% ax_1 = axes;
% im_1 = imagesc(MergedImage); 
% title('Select Points That the Cell Touches')
% hold on
% p_dot = plot(px(:,vM.uFrame),py(:,vM.uFrame),'.r','markersize',10); % plot x y positions found above
% p2_dot = plot(px(1,vM.uFrame),py(1,vM.uFrame),'or'); 
% p_dot_calc = plot(0,0,'.k','markersize',10);
% hold off
% 
% set(im_1,'hittest','off')
% set(p_dot,'hittest','off')
% set(p2_dot,'hittest','off')
% 
% fcn2_1 = @(a,b) set(fig_1,'UserData',find(sqrt((px(:,vM.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,vM.uFrame) - ax_1.CurrentPoint(1,2)).^2) == min(sqrt((px(:,vM.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,vM.uFrame) - ax_1.CurrentPoint(1,2)).^2))));
% fcn2_5 = @(a,b) set(p2_dot,'xdata',px(fig_1.UserData,vM.uFrame),'ydata',py(fig_1.UserData,vM.uFrame));
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
% set(p_dot_calc,'xdata',px(p_calc,vM.uFrame),'ydata',py(p_calc,vM.uFrame));
% 
% but_calc = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
%     'Position',[0.8 0 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
%     'String','CALC','Callback','p_calc(fig_1.UserData) = true; set(p_dot_calc,''xdata'',px(p_calc,vM.uFrame),''ydata'',py(p_calc,vM.uFrame));');
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
cellboundary = [];
fig_ex = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
if exist('vidFrames','var')
    img_ex = imagesc(vidFrames(:,:,1),[min(min(vidFrames(:,:,1))) max(max(vidFrames(:,:,1)))]);
else
    img_ex = imagesc(img_raw);
end
axis image
axis manual
colormap(gray)
hold on
ic = 1;
k = 1
s = celldata(ic).vM.Crop([1 2]);

% plot(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),'.w','markersize',8)
% plot(celldata(ic).Xvector(celldata(ic).celldots) + s(1),celldata(ic).Yvector(celldata(ic).celldots) + s(2),'ow')
plot(celldata(ic).Xvector(celldata(ic).real_points) + s(1),celldata(ic).Yvector(celldata(ic).real_points) + s(2),'.w')
p_ex = plot(celldata(ic).Xloc_k_filt(:,1),celldata(ic).Yloc_k_filt(:,1),'ow');
p_bd = plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'--g');
sf = 0.005;
q_ex = quiver(celldata(ic).Xloc_k_filt(:,1) + s(1),celldata(ic).Yloc_k_filt(:,1) + s(2),celldata(ic).Xdisp_k_dn(:,1),celldata(ic).Ydisp_k_dn(:,1),1,'-w','linewidth',1)
% % q_fttc = quiver(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),sf*celldata(ic).Xforce_k(:,k),sf*celldata(ic).Yforce_k(:,k),0,'-c','linewidth',1)
% q_trpf = quiver(FA.Centroids(:,1),FA.Centroids(:,2),XforceTRPF_k(:,1),YforceTRPF_k(:,1),1,'-c');
tl = title('0');
hold off
set(gca,'YDir','reverse')

k = 0;
while true
    k = k + 1;
    if exist('vidFrames','var')
        set(img_ex,'CData',vidFrames(:,:,k))
    end
    set(p_ex,'XData',celldata(ic).Xloc_k_filt(:,k)+ s(1),'YData',celldata(ic).Yloc_k_filt(:,k) + s(2)) % should be s(1) - 1?
    set(q_ex,'XData',celldata(ic).Xloc_k_filt(:,k) + s(1),'YData',celldata(ic).Yloc_k_filt(:,k) + s(2),'UData',celldata(ic).Xdisp_k_dn(:,k),'VData',celldata(ic).Ydisp_k_dn(:,k))
% %     set(q_fttc,'XData',celldata(ic).Xloc_k_filt(:,k) + s(1),'YData',celldata(ic).Yloc_k_filt(:,k) + s(2),'UData',sf*celldata(ic).Xforce_k(:,k),'VData',sf*celldata(ic).Yforce_k(:,k))
%     set(q_trpf,'UData',XforceTRPF_k(:,k),'VData',YforceTRPF_k(:,k));
    set(tl,'String',sprintf('%d',k))
    drawnow
%     pause(0.1)
    if k >= vM.nFrames
        k = 0;
    end
    if vM.nFrames < 2
        break
    end
end

%% polar histogram
ic = 1;
bin_edges = linspace(0,2*pi,9);
force_rose_sum = zeros(length(bin_edges)-1,vM.nFrames);
[force_angles,force_mag] = cart2pol(celldata(ic).Xforce_k,celldata(ic).Yforce_k);
fangle_in_cell = force_angles(celldata(ic).celldots(:),:) + pi;
fmag_in_cell = force_mag(celldata(ic).celldots(:),:);
for k = 1:vM.nFrames
    for ibin = 1:length(bin_edges)-1
        fmag_in_bin = fmag_in_cell(fangle_in_cell(:,k) > bin_edges(ibin) & fangle_in_cell(:,k) < bin_edges(ibin+1),k);
        fangle_in_bin = fangle_in_cell(fangle_in_cell(:,k) > bin_edges(ibin) & fangle_in_cell(:,k) < bin_edges(ibin+1),k);
    %     [bin_edges(ibin)*180/pi bin_edges(ibin+1)*180/pi length(fmag_in_bin)]
        force_rose_sum(ibin,k) = sum(fmag_in_bin,1);
    end
end

sc = 1e-3; % scale factor for arrows
figure
pax = polaraxes;
for k = 1:vM.nFrames
% for k = celldata(ic).vM.uFrame
    polarhistogram(pax,'BinEdges',bin_edges,'BinCounts',1e-3*force_rose_sum(end:-1:1,k))
    pax.RLim = [0 1e-3*max(force_rose_sum(:))];
    pax.RLimMode = 'manual';
    hold on
%     [x,y] = pol2cart(fangle_in_bin,fmag_in_bin);
%     quiver(zeros(size(x)),zeros(size(y)),x,y,0,'-k')
    polarplot(-[fangle_in_cell(:,k),fangle_in_cell(:,k)]', sc*[zeros(size(fmag_in_cell(:,k))),fmag_in_cell(:,k)]','-k')
    hold off
    title(sprintf('Frame #: %i',k));
    drawnow
end

%% force transient over time
figure('Position',[700,100,600,200])
plot(celldata(ic).vM.Time,10^-3*celldata(ic).total_force,'-k','linewidth',2)
xlabel('Time [s]')
ylabel('Total Force [nN]')
box off
set(gca,'linewidth',1.5,'tickdir','out','XColor','k','YColor','k')

%% plot image with forces for publication
ic = 1;
% k = vM.uFrame;
% k = 100;
k = 1;

s = celldata(ic).vM.Crop([1 2]);

pos = [celldata(ic).Yvector + s(2), celldata(ic).Xvector + s(1)];
loc = [celldata(ic).Yloc_k(:,k) + s(2), celldata(ic).Xloc_k(:,k) + s(1)];
disp = [celldata(ic).Ydisp_k(:,k), celldata(ic).Xdisp_k(:,k)];
forc = [celldata(ic).Yforce_k(:,k), celldata(ic).Xforce_k(:,k)];
trac = [celldata(ic).Ytrac_k(:,k), celldata(ic).Xtrac_k(:,k)];

fig_disp = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
ax_disp = axes(fig_disp,'units','normalized','position',[0 0 1 1]);
% all_cell_dots = any(cat(3,celldata(ic).celldots),3);
imagesc(img_raw,[min(img_raw(:)), max(img_raw(:))])
% imagesc(img_filt)
colormap(gray*[1 0 0;0 130/255 0;0 0 0])
axis image
axis manual
axis off
% axis([435 795 207 575])
hold on
p_bd = plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'-w');

color_max = [1 1 0];
color_min = [1 0 1];

arrowscale = 0.002; % 0.002
arrow_start = [loc(:,2),loc(:,1)];
arrow_end = arrow_start + arrowscale*[forc(:,2), forc(:,1)];

arrowlength = sqrt(sum((arrow_end - arrow_start).^2,2))/2/2;
arrowwidth = arrowlength*5/15;
arrow_color = mat2gray(arrowlength).*(color_max - color_min) + color_min;
q_trac = arrow(arrow_start,arrow_end,...
    'BaseAngle',90,'TipAngle',25,'Width',arrowwidth,'Length',arrowlength,...
    'LineWidth',0.25,'EdgeColor','none','FaceVertexCData',arrow_color,'FaceColor','flat'); % 'FaceColor','w','EdgeColor','none'

ax_scalebars = axes(fig_disp,'position',get(ax_disp,'position'),'hittest','off','pickableparts','none');

scalebar_length = 20/celldata(ic).vM.Calibration;
scalebar_disp = floor(20/arrowscale)*arrowscale/celldata(ic).vM.Calibration;
% scalebar_force = floor(20/celldata(ic).vM.Calibration/(arrowscale/0.001)/10)*(arrowscale/0.001)*10;
% scalebar_force = scalebar_force*2/10;
scalebar_force = 50;
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
patch(ax_scalebars,'XData',[xlims(2) - scalebar_length - 0.2*diff(xlims),xlims(2),xlims(2),xlims(2) - scalebar_length - 0.2*diff(xlims)],...
    'YData',[ylims(2),ylims(2),ylims(2) - 0.15*diff(ylims),ylims(2) - 0.15*diff(ylims)],...
    'EdgeColor','none','FaceColor','k','FaceAlpha',1)

plot(ax_scalebars,[xlims(2) - scalebar_length - 0.15*diff(xlims), xlims(2) - 0.15*diff(xlims)],...
    [ylims(2) - 0.05*diff(ylims), ylims(2) - 0.05*diff(ylims)],'w','linewidth',5)
text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.05*diff(ylims),...
    sprintf('%.0f um',scalebar_length*celldata(ic).vM.Calibration),'color','w','VerticalAlignment','Middle')

arrow([xlims(2) - arrowscale*scalebar_force/0.001 - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
    [xlims(2) - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
    'BaseAngle',90,'TipAngle',25,'Width',max(arrowwidth),'Length',max(arrowlength),'LineWidth',0.25,'FaceColor','w','EdgeColor','none')
% text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
%     sprintf('%.0f nN',1/arrowscale*scalebar_force*0.001),'color','w','VerticalAlignment','Middle')
text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
    sprintf('%.0f nN',scalebar_force),'color','w','VerticalAlignment','Middle')

% quiver(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),celldata(ic).Xtrac_k(:,k),celldata(ic).Ytrac_k(:,k),0.5,'-c','linewidth',1)

% print('arrow_forc','-dsvg','-painters')

%% plot rotated image with forces for publication
ic = 1;
% k = vM.uFrame;
% k = 100;
k = 1;
aa = celldata(ic).vM.rot_angle;
s = celldata(ic).vM.Crop([1 2]);
img_rotated = imrotate(img_raw,aa*180/pi,'crop');

R = [cos(aa), sin(aa); -sin(aa), cos(aa)];
% R = [1 0; 0 1];
centerX = floor(celldata(ic).vM.N/2+1) + s(1);
centerY = floor(celldata(ic).vM.M/2+1) + s(2);
pos = [celldata(ic).Yvector + s(2), celldata(ic).Xvector + s(1)];
loc = [celldata(ic).Yloc_k(:,k) + s(2), celldata(ic).Xloc_k(:,k) + s(1)];
disp = [celldata(ic).Ydisp_k(:,k), celldata(ic).Xdisp_k(:,k)];
forc = [celldata(ic).Yforce_k(:,k), celldata(ic).Xforce_k(:,k)];
pos_rot = (pos - [centerY, centerX])*R + [centerY, centerX];
loc_rot = (loc - [centerY, centerX])*R + [centerY, centerX];
disp_rot = disp*R;
forc_rot = forc*R;

fig_disp = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
ax_disp = axes(fig_disp,'units','normalized','position',[0 0 1 1]);
% all_cell_dots = any(cat(3,celldata(ic).celldots),3);
imagesc(img_rotated,[min(img_raw(:)), max(img_raw(:))])
% imagesc(img_filt)
colormap(gray*[1 0 0;0 130/255 0;0 0 0])
axis image
axis manual
axis off
% axis([435 795 207 575])
hold on
p_bd = plot(celldata(ic).CB(:,1) + s(1),celldata(ic).CB(:,2) + s(2),'--g');
% plot(pos_rot(:,2),pos_rot(:,1),'+k','markersize',10)
% plot(loc_rot(:,2),loc_rot(:,1),'.k','markersize',25)
% quiver(pos_rot(:,2),pos_rot(:,1),disp_rot(:,2),disp_rot(:,1),0,'-w','linewidth',1)
% arrow([pos_rot(:,2),pos_rot(:,1)],[pos_rot(:,2),pos_rot(:,1)]+[disp_rot(:,2),disp_rot(:,1)],...
%     'Length',5)

arrowscale = 0.002;
arrow_start = [loc_rot(:,2),loc_rot(:,1)];
arrow_end = arrow_start + arrowscale*[forc_rot(:,2), forc_rot(:,1)];

arrowlength = sqrt(sum((arrow_end - arrow_start).^2,2))/2;
arrowwidth = arrowlength*5/15;
q_trac = arrow(arrow_start,arrow_end,...
    'BaseAngle',90,'TipAngle',25,'Width',arrowwidth,'Length',arrowlength,...
    'LineWidth',0.25,'FaceColor','w','EdgeColor','none');

ax_scalebars = axes(fig_disp,'position',get(ax_disp,'position'),'hittest','off','pickableparts','none');

scalebar_length = 20/celldata(ic).vM.Calibration;
scalebar_disp = floor(20/arrowscale)*arrowscale/celldata(ic).vM.Calibration;
% scalebar_force = floor(20/celldata(ic).vM.Calibration/(arrowscale/0.001)/10)*(arrowscale/0.001)*10;
% scalebar_force = scalebar_force*2/10;
scalebar_force = 20;
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
patch(ax_scalebars,'XData',[xlims(2) - scalebar_length - 0.2*diff(xlims),xlims(2),xlims(2),xlims(2) - scalebar_length - 0.2*diff(xlims)],...
    'YData',[ylims(2),ylims(2),ylims(2) - 0.15*diff(ylims),ylims(2) - 0.15*diff(ylims)],...
    'EdgeColor','none','FaceColor','k','FaceAlpha',1)

plot(ax_scalebars,[xlims(2) - scalebar_length - 0.15*diff(xlims), xlims(2) - 0.15*diff(xlims)],...
    [ylims(2) - 0.05*diff(ylims), ylims(2) - 0.05*diff(ylims)],'w','linewidth',5)
text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.05*diff(ylims),...
    sprintf('%.0f um',scalebar_length*celldata(ic).vM.Calibration),'color','w','VerticalAlignment','Middle')

arrow([xlims(2) - arrowscale*scalebar_force/0.001 - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
    [xlims(2) - 0.15*diff(xlims), ylims(2) - 0.1*diff(ylims)],...
    'BaseAngle',90,'TipAngle',25,'Width',max(arrowwidth),'Length',max(arrowlength),'LineWidth',0.25,'FaceColor','w','EdgeColor','none')
% text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
%     sprintf('%.0f nN',1/arrowscale*scalebar_force*0.001),'color','w','VerticalAlignment','Middle')
text(ax_scalebars,xlims(2) - 0.13*diff(xlims), ylims(2) - 0.1*diff(ylims),...
    sprintf('%.0f nN',scalebar_force),'color','w','VerticalAlignment','Middle')

% quiver(celldata(ic).Xvector + s(1),celldata(ic).Yvector + s(2),celldata(ic).Xtrac_k(:,k),celldata(ic).Ytrac_k(:,k),0.5,'-c','linewidth',1)

% print('arrow_forc','-dsvg','-painters')