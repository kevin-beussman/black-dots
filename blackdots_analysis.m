%% black dots analysis
clearvars; clc; close all;

uM.force_user_vals = false; % set this to 'true' to if you want to use the below values instead of the ND2 values
uM.analyzeVideo = false; % are you analyzing a video or an image? set true for video

% these following settings are used only if they cannot be pulled from the video data
uM.CameraPixelSize = 3.63; % physical pixel size [um], from camera manufacturer V3=6.5, Flash2.8=3.63, Andor=6.45
uM.Objective = 60; % e.g. 60x objective
uM.Binning = 1; % e.g. 1x1 binning
uM.FrameRate = 20; % [Hz]
uM.CouplerRatio = 1; % e.g. 0.7x coupler
uM.MagMultiplier = 1; % e.g. 1.5x magnification toggle

% video settings (ignored if analyzeVideo is false)
% uM.Frames = 1:50; % video frames to analyze. comment this out to use all frames
uM.uFrame = 1; % which video frame to use for relaxed frame -- this is relative to uM.Frames

% image settings (ignored if analyzeVideo is true)
uM.BDchannel = 2; % channel to use for black dots
uM.CBchannel = 1; % optional channel to use for cell boundary. comment this out if there is only one channel

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

% image boundary and cropping
uM.Crop = false; % true or give crop boundary[xmin ymin, width height]
% uM.crop_around_cells = true;
uM.manual_boundary = true;

% pattern tracking and force calculation settings
uM.bright_dots = false; % this code tracks dark dots. Set "true" here to invert image before analyzing
uM.track_method = 2; % method 1 = more manual method, method 2 = more automated method
uM.uPoints = 2; % number of outer undeformed points in each direction to use for grid estimation. comment out to use all points
uM.useLcurve = false; % calculate regularization parameter from L curve (takes a while)
uM.regParam = 5e-8; % regularization parameter guess

%% Choose file to analyze
file_script = mfilename('fullpath');
[path_script,~,~] = fileparts(file_script);
addpath(path_script)
addpath([path_script filesep 'bftools'])
addpath([path_script filesep 'bfmatlab'])
addpath([path_script filesep 'track'])
addpath([path_script filesep 'blackdots_functions'])

[file_BD, path_BD] = uigetfile({'*.nd2','NIS Elements';'*.cxd','HCImage Live';'*.avi','(AVI) Audio Video Interleave';'*.tif;*.tiff','TIFF Image Stacks'},'Select file to analyze. NOTE: .nd2 files require bioformats plugin.','MultiSelect','off');
cd(path_BD)

if uM.analyzeVideo
    [img_BD,meta_BD] = read_video([path_BD, file_BD],uM);
else
    if isfield(uM,'CBchannel')
        [img_BD,meta_BD,img_CB] = read_image([path_BD, file_BD],uM);
    else
        [img_BD,meta_BD] = read_image([path_BD, file_BD],uM);
    end
end
if uM.bright_dots
    img_BD = imcomplement(img_BD); % use this if dots are bright instead of dark
end
img_BD_raw = img_BD;
meta_BD_raw = meta_BD;

[~,filename,~] = fileparts(file_BD);

% create save folder(s)
cur_time = char(datetime('now','format','yyyyMMdd_HHmmss'));

path_save = [path_BD filename '_blackdots_analysis' filesep];
if ~exist(path_save,'dir')
    mkdir(path_save)
end

path_save = [path_save cur_time filesep];
if ~exist(path_save,'dir')
    mkdir(path_save)
end

file_save = ['blackdots_data_' cur_time '.mat'];

%% Find cell boundary
[file_CB_data, path_CB_data] = uigetfile({'*.mat'},'Select mat file containing cell boundary coordinates (a previous blackdots_analysis file).');
if file_CB_data ~= 0
    load([path_CB_data file_CB_data],'CB_uncrop','img_CB');
    
    if meta_BD.Crop ~= false
        CB = cellfun(@plus,CB_uncrop,{-meta_BD.Crop([1 2])},'UniformOutput',false);
    else
        CB = CB_uncrop;
    end
    
     save([path_save file_save],'CB_uncrop','img_CB')
else
    if meta_BD.analyzeVideo || ~isfield(meta_BD,'CBchannel')  % look for another image file to get the cell boundary
        [file_CB, path_CB] = uigetfile({'*.nd2','NIS Elements';'*.cxd','HCImage Live';'*.avi','(AVI) Audio Video Interleave';'*.tif;*.tiff','TIFF Image Stacks'},'Select image file to use for tracing cell boundary.','MultiSelect','off');
        if file_CB ~= 0
            [img_CB,meta_CB] = read_image([path_CB file_CB]);
            if meta_BD.manual_boundary
                CB = get_cell_boundary(img_CB); % manual clicking boundary
            else
                CB = get_cell_boundary2(img_CB,meta_BD); % semi-automatic segmentation
            end
        else % if no other image is selected, default crop is to just use the whole image
            CB{1} = [1,1; meta_BD.N,1; meta_BD.N,meta_BD.M; 1,meta_BD.M];
            img_CB = img_BD_raw(:,:,1);
            img_CB(:,:) = 0;
        end
        CBraw = CB;
        img_CB_raw = img_CB;
        if meta_CB.Binning ~= meta_BD.Binning % CB and BD images have different binning/sizes, do some correction
            for ic = 1:length(CB)
                CB{ic} = CB{ic}*(meta_CB.Binning/meta_BD.Binning);
            end
            img_CB = imresize(img_CB,meta_CB.Binning/meta_BD.Binning);
        end
    else % we already got the cell boundary image from current video/image in CBchannel
        if meta_BD.manual_boundary
            CB = get_cell_boundary(img_CB); % manual clicking boundary
        else
            CB = get_cell_boundary2(img_CB,meta_BD); % semi-automatic segmentation
    %         stats = get_cell_boundary_zizhen(cbFrame,meta_BD); % doesnt work on kb2 computer
        end
    end

    if meta_BD.Crop ~= false
        CB_uncrop = cellfun(@plus,CB,{meta_BD.Crop([1 2])},'UniformOutput',false);
    else
        CB_uncrop = CB;
    end
    
    save([path_save file_save],'CB_uncrop','img_CB')
end

%% Bandpass filter image or all video frames
fprintf('Filtering black dots image(s)...')

bpass_noise = 1; % bandpass characteristic noise length (leave at 1 usually)
dotsize_px = 2*round((meta_BD.DotSize/meta_BD.Calibration + 1)/2) - 1; % dot size in pixels
img_BD_filt = zeros(size(img_BD));
for k = 1:meta_BD.nFrames
    img_BD_filt(:,:,k) = mat2gray(imcomplement(bpass_kb2(imcomplement(img_BD(:,:,k)),bpass_noise,dotsize_px)));
end

fprintf('DONE\n')

%% start analysis for each cell
nCells = length(CB);
fprintf('\t%i objects (cells) detected\n',nCells)
skipped = false(nCells,1);
celldata(nCells) = struct();
for ic = 1:nCells
    %% crop to selected cell
    fprintf('Cropping black dots image(s) around cell...')

    % expands crop region to include some undeformed dots outside cell
    celldata(ic).crop_tight = [min(CB{ic}), max(CB{ic})-min(CB{ic})];
    celldata(ic).crop = celldata(ic).crop_tight;
    last_crop = celldata(ic).crop;

    expand = meta_BD.DotSize/meta_BD.Calibration*[-1 -1 2 2]; % crop boundary expands by this much every time "expand" is clicked
    
    fig_cellselect = figure('Units','Normalized','Position',[0.2 0.1 0.6 0.6]);
    ax1 = subplot(1,2,1);
    im1 = imagesc(img_BD(:,:,meta_BD.uFrame));
    axis image
    hold on
    plot(CB{ic}(:,1),CB{ic}(:,2),'-r')
    hold off
    colormap(ax1,gray*[1 0 0; 0 0.75 0; 0 0 0])
    
    ax2 = subplot(1,2,2);
    imagesc(img_CB);
    axis image
    hold on
    plot(CB{ic}(:,1),CB{ic}(:,2),'-r')
    hold off
    colormap(ax2,gray*[0 0 0; 0 0.75 0; 0 0 1])
    linkaxes([ax1 ax2])
    axis(celldata(ic).crop_tight([1 1 2 2]) + [0, celldata(ic).crop_tight(3), 0, celldata(ic).crop_tight(4)] - 0.5)
    
    but_expand = uicontrol(fig_cellselect,'Units','Normalized','Position',[0 0 0.2 0.1],'Style','PushButton','String','expand',...
        'Callback',['celldata(ic).crop = celldata(ic).crop + expand;', ...
                    'if celldata(ic).crop(1) < 1, celldata(ic).crop(1) = 1; end;', ...
                    'if celldata(ic).crop(2) < 1, celldata(ic).crop(2) = 1; end;', ...
                    'if sum(celldata(ic).crop([2 4])) > meta_BD.M, celldata(ic).crop(4) = meta_BD.M - celldata(ic).crop(2) + 1; end;', ...
                    'if sum(celldata(ic).crop([1 3])) > meta_BD.N, celldata(ic).crop(3) = meta_BD.N - celldata(ic).crop(1) + 1; end;', ...
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
    
    if skipped(ic) % if this cell is to be skipped, immediately move to next cell
        continue
    end
    
    % note: using imcrop is slow compared to this
    % xrange and yrange are just the range of rows/columns to keep
    % (identical to what you get out of imcrop)
    % celldata(ic).crop is float, so we need to round to nearest pixel
    xrange = celldata(ic).crop([2,4]); xrange(2) = round(xrange(2) + xrange(1) - 1); xrange(1) = round(xrange(1));
    yrange = celldata(ic).crop([1,3]); yrange(2) = round(yrange(2) + yrange(1) - 1); yrange(1) = round(yrange(1));
    
    img_BD_crop = img_BD(xrange(1):xrange(2),yrange(1):yrange(2),:);
    img_BD_filt_crop = img_BD_filt(xrange(1):xrange(2),yrange(1):yrange(2),:);
    
    img_REFBD = img_BD(:,:,meta_BD.uFrame);
    img_REFBD_crop = img_BD_crop(:,:,meta_BD.uFrame);
    img_REFBD_filt_crop = img_BD_filt_crop(:,:,meta_BD.uFrame);
    
    celldata(ic).M = size(img_BD_filt_crop,1);
    celldata(ic).N = size(img_BD_filt_crop,2);
    celldata(ic).CB = CB{ic} - celldata(ic).crop([1 2]);
    
    fprintf('DONE\n')
    
    %% Start analysis
    if meta_BD.track_method == 1 % more manual method of tracking dots
        %% Measure rotation angle from image
        fprintf('Measuring rotation angle...')

        DotSpacing_px = meta_BD.DotSpacing/meta_BD.Calibration;
        
        celldata(ic).rot_angle = get_rot_from_img(img_REFBD_filt_crop,DotSpacing_px);
    
        fprintf('DONE\n')
        
        %% Choose the sampling dots
        % allows user to estimate dot grid positions
        fprintf('Choosing dots...')
    
        [px_sample,py_sample] = choose_dots(img_REFBD_filt_crop,celldata(ic),meta_BD);
        
        fprintf('DONE\n')
        
        if isempty(px_sample)
            skipped(ic) =  true;
            continue
        end
                
        %% Update rotation angle
        fprintf('Updating rotation angle...')
    
        celldata(ic).rot_angle = get_rot_from_gridpts(cat(3,px_sample,py_sample));
    
        fprintf('DONE\n')
        
        %% Extend the sampling dots
        % fits a grid to the sampling dots, and extends the grid of dots
        % to fill the full image
        fprintf('Extending dots to full image...')
        
        celldata.real_points = []; % for debugging
        [px,py,real_points,img_REFBD_bw] = extend_dots(px_sample,py_sample,img_REFBD_filt_crop,celldata(ic),meta_BD);

        celldata(ic).real_points = real_points;

        fprintf('DONE\n')

        %% Characterize dots
        fprintf('Characterizing dots...')
        
        img = img_REFBD_filt_crop;
        img_filt = imgaussfilt(img,1);
        img_filt_bw = imfill(imcomplement(imbinarize(img_filt)),'holes');
        img_REFBD_bw = bwselect(img_filt_bw,px(real_points),py(real_points));
        
        BD = calc_dot_size_spacing(px,py,img_REFBD_bw,celldata(ic),meta_BD);
        
        celldata(ic).BD = BD;
        
        fprintf('DONE\n')
        
        %% Find undeformed dot positions
        fprintf('Finding undeformed dot centroids...')
    
        [px,py,px0,py0] = find_undeformed(px,py,celldata(ic),meta_BD);
        
        fprintf('DONE\n')
        
        %% Track dots across all frames
        % just basically repeats find_centroids for each frame
        fprintf('Finding dot centroids in each image frame...')
    
        % this should really call find_centroids for each frame
        if meta_BD.analyzeVideo
            [px_k,py_k,real_points] = track_dots_across_frames(img_BD_filt_crop,px,py,celldata(ic),meta_BD);

            celldata(ic).real_points = real_points; %just in case some dots are lost when analyzing
        else
            px_k = px(:);
            py_k = py(:);
        end

%         figure
%         for k = meta_BD.Frames
%             imagesc(img_BD_filt_crop(:,:,k))
%             hold on
%             plot(px_k(:,:,k),py_k(:,:,k),'.r')
%             hold off
%             xlim([100,200])
%             ylim([100,200])
%             title(num2str(k))
%             pause
%         end
    
        fprintf('DONE\n')
        
    elseif meta_BD.track_method == 2
        % this method uses the "track" package from http://site.physics.georgetown.edu/matlab/
        if ~exist('track','file')
            error('''track'' not found, visit http://site.physics.georgetown.edu/matlab/')
        end

        %% Find object centers
        fprintf('Finding dot centers in each image frame...')
    
        dotsize = 2*round((meta_BD.DotSize/meta_BD.Calibration + 1)/2) - 1; % pixels, nearest odd integer
        dotspacing = round(meta_BD.DotSpacing/meta_BD.Calibration);
    
        bpass_noise = 1; % bandpass noise: 1 pixel is typical      % 1*meta_BD.Calibration
        pkfnd_threshold = 0.15; % object peak threshold: find by trial and error, 0.15 is good
    
        
        centers = []; % centers has 3 columns: [x, y, frame#]. we don't know how many objects there will be
        pct = 0;
        fprintf('[%-20s] %3.0f%%\n',repmat('|',1,round(pct/100*20)),pct)
        for k = 1:meta_BD.nFrames
            pct = k/meta_BD.nFrames;
            fprintf([repmat('\b',1,28) '[%-20s] %3.0f%%\n'],repmat('|',1,round(pct*20)),pct*100)
            I1 = img_BD_crop(:,:,k);
            I2 = bpass_kb2(imcomplement(I1),bpass_noise,dotsize); % spatial bandpass filter
            pk = pkfnd_kb2(mat2gray(I2),pkfnd_threshold,dotsize); % find objects whose intensity is above threshold
            dsz = ceil(dotsize*1.2); if ~mod(dsz,2), dsz = dsz + 1; end % KB2 edit 9/3/2020
            cnt = cntrd_kb2(I2,pk,dsz); % KB2 edit 9/3/2020
    
            centers = [centers; cnt(:,1:2), ones(size(cnt,1),2)*k];
        end
        fprintf(repmat('\b',1,28))
    
%         figure
%         imagesc(img_REFBD_crop)
%         hold on
%         plot(centers(:,1),centers(:,2),'.r','markersize',10)
%         hold off
        
        fprintf('DONE\n')
        
        %% Track objects across video frames
        fprintf('Tracking centers...')
    
        if meta_BD.analyzeVideo
    
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
    
        %% Fill in gaps in time for tracked points
        fprintf('Filling in temporal gaps...')
        
        if meta_BD.analyzeVideo
            [px,py] = fill_in_temporal_gaps(res,celldata(ic),meta_BD);
        else
            px = res(:,1);
            py = res(:,2);
        end
    
        fprintf('DONE\n')
        
        %% Manual point removal
        fig_1 = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
        ax_1 = axes;
        im_1 = imagesc(imadjust(mat2gray(img_BD_crop(:,:,meta_BD.uFrame))));
        axis image
    %     colormap(gray)
        hold on
        p_dot = plot(px(:,meta_BD.uFrame),py(:,meta_BD.uFrame),'.r','markersize',10);
        p2_dot = plot(px(1,meta_BD.uFrame),py(1,meta_BD.uFrame),'or');
        p_dot_remove = plot(0,0,'.k','markersize',10);
        hold off
    
        set(im_1,'hittest','off')
        set(p_dot,'hittest','off')
        set(p2_dot,'hittest','off')
    
        fcn1_1 = @(a,b) set(fig_1,'UserData',find(sqrt((px(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,2)).^2) == min(sqrt((px(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,2)).^2))));
        fcn1_2 = @(a,b) set(p2_dot,'xdata',px(fig_1.UserData,meta_BD.uFrame),'ydata',py(fig_1.UserData,meta_BD.uFrame));
        fcn1 = @(a,b) cellfun(@feval,{fcn1_1 fcn1_2});
    
        fcn2_1 = @(a,b) set(fig_1,'UserData',[fig_1.UserData; find(sqrt((px(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,2)).^2) == min(sqrt((px(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,1)).^2 + (py(:,meta_BD.uFrame) - ax_1.CurrentPoint(1,2)).^2)))]);
        fcn2_2 = @(a,b) set(p2_dot,'xdata',px(fig_1.UserData,meta_BD.uFrame),'ydata',py(fig_1.UserData,meta_BD.uFrame));
        fcn2 = @(a,b) cellfun(@feval,{fcn2_1 fcn2_2});
    
        fcn3 = @(a,b) set(ax_1,'buttondownfcn',fcn2);
        fcn4 = @(a,b) set(ax_1,'buttondownfcn',fcn1);
    
        set(ax_1,'buttondownfcn',fcn1)
        set(fig_1,'KeyPressFcn',fcn3,'KeyReleaseFcn',fcn4)
    
        p_remove = false(size(px));
        but_remove = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
            'Position',[0.8 0 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
            'String','REMOVE','Callback','p_remove(fig_1.UserData) = true; set(p_dot_remove,''xdata'',px(p_remove,meta_BD.uFrame),''ydata'',py(p_remove,meta_BD.uFrame));');
    
        but_done = uicontrol('parent',fig_1,'units','normalized','Style','pushbutton',...
            'Position',[0.8 0.95 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
            'String','DONE','Callback','uiresume');
    
        uiwait
    
        px(p_remove,:) = [];
        py(p_remove,:) = [];
    
        close(fig_1)
        
        %% get rotation angle based on dot locations
        fprintf('Finding rotation angle...')
        
        DotSpacing_px = meta_BD.DotSpacing/meta_BD.Calibration;
        
        if meta_BD.analyzeVideo
            celldata(ic).rot_angle = get_rot_from_pts([px(:,meta_BD.uFrame),py(:,meta_BD.uFrame)],DotSpacing_px);
        else
            celldata(ic).rot_angle = get_rot_from_pts([px,py],DotSpacing_px);
        end
        fprintf('DONE\n')
        
        %% find undeformed state of grid
        fprintf('Finding undeformed state of grid...')
        
        [px_k,py_k,px0,py0,real_points] = find_undeformed_grid(px,py,celldata(ic),meta_BD); % IF YOU HAVE A LOT OF DISPLACEMENT, CHANGE INF TO 2 to find fewer points
        px = px_k(:,:,meta_BD.uFrame);
        py = py_k(:,:,meta_BD.uFrame);
        celldata(ic).real_points = real_points;
        
        fprintf('DONE\n')
        
        %% calculate dot size and spacing from image
        fprintf('Characterizing dots...')
        
        BD = calc_dot_size_spacing(px,py,img_REFBD_bw,celldata(ic),meta_BD);
        
        celldata(ic).BD = BD;
        
        fprintf('DONE\n')
        
    end
    
    %% Filter locations
    fprintf('Filtering dot positions...')
    
    Xloc_k = reshape(px_k,[],meta_BD.nFrames);
    Xvector = reshape(px0,[],1);
    
    Yloc_k = reshape(py_k,[],meta_BD.nFrames);
    Yvector = reshape(py0,[],1);
    
    % filter
    uHz = 4; % lowpass cutoff
    % lHz = 0.5;
    if meta_BD.nFrames >= 2
    %     [B,A] = butter(3,[lHz/(meta_BD.FrameRate/2), uHz/(meta_BD.FrameRate/2)],'bandpass');
        [B,A] = butter(3,uHz/(meta_BD.FrameRate/2),'low');
    
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
    %     totaldisp(f) = norm(XYdisp(:,f))*meta_BD.Calibration;
    % end
    
    fprintf('DONE\n')
    
    %% find cell points
    
    % CB = CBraw;
    % nCells = length(CBraw);
    % fprintf('\t%i cells detected\n',nCells)
    % edge_tooclose = 2*mean(BD.DotSpacings)/meta_BD.Calibration;
    dot_consider = 0.75*mean(celldata(ic).BD.DotSpacings)/meta_BD.Calibration;
    
    tb = celldata(ic).CB; % temporary cell boundary data
    %     fprintf('\tCell %i...',nc2)
    
    %     %check that cell is not too close to outside
    %     if any([(min(tb(:,1)) - 1), (min(tb(:,2)) - 1), (meta_BD.N - max(tb(:,1))), (meta_BD.M - max(tb(:,2)))] < edge_tooclose)
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
    fprintf('Removing points below noise...')
    
    Xdisp_all = Xdisp_k_filt;
    Ydisp_all = Ydisp_k_filt;
    % Xdisp_out = Xdisp_k_filt(real_points(:) & ~celldots(:));
    % Ydisp_out = Ydisp_k_filt(real_points(:) & ~celldots(:));
    
    disp_all = sqrt(Xdisp_all.^2 + Ydisp_all.^2)*meta_BD.Calibration;
    % disp_out = sqrt(Xdisp_out.^2 + Ydisp_out.^2)*meta_BD.Calibration;
    
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
    
    fprintf('DONE (%i points of %i removed)\n',nnz(noise_points),nnz(celldata(ic).real_points))
    
    %% calculate traction forces
    fprintf('Calculating tractions...')
    
    [n_row,n_col] = size(Xgrid);
    
    E = meta_BD.YoungsModulus;
    nu = meta_BD.Poisson;
    d = mean(celldata(ic).BD.DotSpacings);
    
    Xtrac_k = zeros(size(Xdisp_k_filt));
    Ytrac_k = zeros(size(Xdisp_k_filt));
    Xforce_k = zeros(size(Xdisp_k_filt));
    Yforce_k = zeros(size(Xdisp_k_filt));
    
    ind_k = [(meta_BD.uFrame:-1:1) (meta_BD.uFrame+1:meta_BD.nFrames)];
    
    % FTTC with Regularization
    % for k = 1:meta_BD.nFrames
    for kk = 1:length(ind_k)
        k = ind_k(kk);
        u_x = reshape(Xdisp_k_dn(:,k),n_row,n_col)*meta_BD.Calibration;
        u_y = reshape(Ydisp_k_dn(:,k),n_row,n_col)*meta_BD.Calibration;
        u = cat(3,u_x,u_y);
        if kk == 1
            if meta_BD.useLcurve
                L = logspace(-3,-8,100);
    %             [~,resid_norm(k,:),soln_norm(k,:)] = calcforce_regFTTC(u,E,nu,d,L,real_points);
    %             [reg_optimal(k),i_reg_optimal(k)] = find_L(resid_norm(k,:),soln_norm(k,:),L,meta_BD.regParam,'optimal');
    %             [reg_corner(k),i_reg_corner(k)] = find_L(resid_norm(k,:),soln_norm(k,:),L,meta_BD.regParam,'corner');
                [~,resid_norm,soln_norm] = calcforce_regFTTC(u,E,nu,d,L,celldata(ic));
                [reg_optimal,i_reg_optimal] = find_L(resid_norm,soln_norm,L,meta_BD.regParam,'optimal');
                [reg_corner,i_reg_corner] = find_L(resid_norm,soln_norm,L,meta_BD.regParam,'corner');
    %             L = reg_optimal(k);
                L = meta_BD.regParam;
                
    %             loglog(resid_norm,soln_norm,'-k',resid_norm(i_reg_optimal),soln_norm(i_reg_optimal),'or',resid_norm(i_reg_corner),soln_norm(i_reg_corner),'ob')
            else
                L = meta_BD.regParam;
                reg_optimal = NaN;
                reg_corner = NaN;
            end
        end
    
        trac = calcforce_regFTTC(u,E,nu,d,L,celldata(ic));
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
    
    fprintf('DONE\n')
    
    %% calculate data for each cell
    fprintf('Calculating data for each cell...')
    
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
    if meta_BD.useLcurve
        celldata(ic).resid_norm = resid_norm;
        celldata(ic).soln_norm = soln_norm;
    end
    
    celldata(ic).celldots = celldots;
    
    celldata(ic).total_trac = sum(sqrt(Xtrac_k(celldots(:),:).^2 + Ytrac_k(celldots(:),:).^2)); % pN*um^-2
    
    celldata(ic).total_force = sum(sqrt(Xforce_k(celldots(:),:).^2 + Yforce_k(celldots(:),:).^2)); % pN
    celldata(ic).net_force = sqrt(sum(Xforce_k(celldots(:),:)).^2 + sum(Yforce_k(celldots(:),:)).^2); % pN
    
    celldata(ic).total_disp = sum(sqrt(Xdisp_k_dn(celldots(:),:).^2 + Ydisp_k_dn(celldots(:),:).^2)); % um
    celldata(ic).net_disp = sqrt(sum(Xdisp_k_dn(celldots(:),:)).^2 + sum(Ydisp_k_dn(celldots(:),:)).^2); % um
    
    celldata(ic).nDots = nnz(celldots);
    celldata(ic).area = polyarea(tb(:,1),tb(:,2))*meta_BD.Calibration.^2; % [um^2]
    
    celldata(ic).img_ref = img_REFBD_crop;
    celldata(ic).img_ref_filt = img_REFBD_filt_crop;
    % figure
    % plot(px,py,'.k')
    % set(gca,'ydir','reverse')
    % hold on
    % plot(CB{1}(:,1),CB{1}(:,2),'-r')
    % plot(px(cell_dots{nc}),py(cell_dots{nc}),'ok')
    % hold off
    fprintf('DONE\n')

end % do next cell

%% save the data
fprintf('Saving...')

% cur_time = char(datetime('now','format','yyyyMMdd_HHmmss'));

% file_save = 'blackdots_data.mat';
% file_save = ['blackdots_data_' cur_time '.mat'];

save([path_save file_save],...
    'file_BD','path_BD','uM','meta_BD',...
    'img_REFBD*',...
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
%     'cellboundary','file_cb','path_cb','meta_BD_cb',...
%     'tform_data','file_stain','path_stain','meta_BD_stain')

% path_backup = [path_save 'backup\'];
% if ~exist(path_backup,'dir')
%     mkdir(path_backup)
% end
% 
% file_backup = ['blackdots_data_' cur_time '.mat'];
% 
% copyfile([path_save file_save],[path_backup file_backup])

fprintf('DONE\n')

%% Plotting results
arrowscale = 0.005; % 0.002 might need to play around with this one
scalebar_length = 5; % microns
scalebar_force = 10; % nN

plot_forces
