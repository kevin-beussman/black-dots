function [vidFrames,vidMetadata] = read_video(file_raw,userMetadata)
    [path_name, file_name, file_ext] = fileparts(file_raw);
    
    vidMetadata = userMetadata;
    vidMetadata.PathName = path_name;
    vidMetadata.FileName = [file_name, file_ext];

    if ~isfield(vidMetadata,'MagMultiplier')
        vidMetadata.MagMultiplier = 1;
    end

    if ~isfield(vidMetadata,'verbose')
        vidMetadata.verbose = 1;
    end
    
    fprintf('Reading file < %s >\n',[file_name, file_ext])
    
    load_image

    crop_image

    get_camera_pixel_size
    
    get_coupler_ratio
    
    get_objective
    
    get_binning
    
    get_timestamps

    get_calibration

    display_metadata

    %% functions listed above
    function load_image
        % load video frames
        if (strcmp(file_ext,'.nd2') || strcmp(file_ext,'.cxd'))
            if exist('bfopen_kb2','file')
                if isfield(vidMetadata,'Frames')
                    rawdata = bfopen_kb2(file_raw,vidMetadata.Frames);
                else
                    rawdata = bfopen_kb2(file_raw);
                end
            elseif exist('bfopen','file')
                rawdata = bfopen(file_raw);
            else
                error('CANNOT READ FILE: Please ensure you have the bioformats functions installed.')
            end
            rawFrames = rawdata{1};
            rawFrames = rawFrames(~cellfun(@isempty,rawFrames(:,1)),1);
            % get raw metadata
            vidMetadata.rawMetadata = rawdata{2};
            
            vidMetadata.nFrames = size(rawFrames,1);
            for k = 1:vidMetadata.nFrames
                vidFrames(:,:,k) = rawFrames{k,1};
            end
        
        elseif strcmp(file_ext,'.tif') || strcmp(file_ext,'.tiff')
            if isfield(vidMetadata,'Frames')
                for k = 1:length(vidMetadata.Frames)
                    vidFrames(:,:,k) = imread(file_raw,vidMetadata.Frames(k));
                end
                vidMetadata.nFrames = size(vidFrames,3);
            else
                vid_info = imfinfo(file_raw);
                vidMetadata.nFrames = size(vid_info,1);
                for k = 1:vidMetadata.nFrames
                    vidFrames(:,:,k) = imread(file_raw,k);
                end
            end
        
        elseif strcmp(file_ext,'.avi')
            v = VideoReader(file_raw);
            if isfield(vidMetadata,'Frames')
                i = 1;
                k = 0;
                while hasFrame(v)
                    k = k + 1;
                    temp = readFrame(v);
                    if vidMetadata.Frames(i) == k
                        vidFrames(:,:,i) = temp;
                    end
                end
            else
                k = 0;
                while hasFrame(v)
                    k = k + 1;
                    vidFrames(:,:,k) = readFrame(v);
                end
            end
            vidMetadata.nFrames = size(vidFrames,3);
            if ~vidMetadata.force_user_vals
                vidMetadata.FrameRate = v.FrameRate;
            end
        else
            error('UNSUPPORTED VIDEO FORMAT: Please use .nd2, .cxd, .avi, or .tif stacks.')
        end
        [vidMetadata.M,vidMetadata.N] = size(vidFrames(:,:,1));
    end
    
    function crop_image
        if isfield(vidMetadata,'Crop') && all(vidMetadata.Crop ~= false)
            if vidMetadata.Crop == true
                figure(1)
                title('CROP IMAGE')
                [~,vidMetadata.Crop] = imcrop(mat2gray(vidFrames(:,:,1)));
                close(1)
            end
            
            for k = 1:vidMetadata.nFrames
                vidFrames_crop(:,:,k) = imcrop(vidFrames(:,:,k),vidMetadata.Crop);
            end
            vidFrames = vidFrames_crop;
            [vidMetadata.M,vidMetadata.N] = size(vidFrames(:,:,1));
        end
    end
    
    function get_camera_pixel_size
        % get physical camera pixel value
        if (strcmp(file_ext,'.nd2') || strcmp(file_ext,'.cxd'))

            vidMetadata.CameraName = vidMetadata.rawMetadata.get('Global Camera Name');
            if isempty(vidMetadata.CameraName)
                vidMetadata.CameraName = vidMetadata.rawMetadata.get('Global CameraUserName');
            end
            if isempty(vidMetadata.CameraName)
                vidMetadata.CameraName = vidMetadata.rawMetadata.get('Global CameraName');
            end
            if isempty(vidMetadata.CameraName)
                vidMetadata.CameraName = vidMetadata.rawMetadata.get('Camera Name');
            end
            
            if ~vidMetadata.force_user_vals
                if strcmp(vidMetadata.CameraName,'Andor Clara DR-1357') % Andor
                    vidMetadata.CameraPixelSize = 6.45; % physical pixel size in userMetadata
                elseif strcmp(vidMetadata.CameraName,'C11440-10C') % Hamamatsu Flash2.8
                    vidMetadata.CameraPixelSize = 3.63; % physical pixel size in userMetadata
                elseif strcmp(vidMetadata.CameraName,'Flash4.0, SN:301638') % Hamamatsu Flash4.0 300017
                    vidMetadata.CameraPixelSize = 6.5; % physical pixel size in userMetadata
                elseif strcmp(vidMetadata.CameraName,'Nikon A1plus') % Garvey Confocal A1
                    vidMetadata.CameraPixelSize = 3.8529; % physical pixel size in userMetadata
                else
                    vidMetadata.CameraName = 'Unknown Camera';
                    % use user submitted CameraPixelSize
                end
            end
        else
            vidMetadata.CameraName = 'Unknown Camera';
            % use user submitted CameraPixelSize
        end
    end
    
    function get_coupler_ratio
        % get coupler ratio value
        if (strcmp(file_ext,'.nd2') || strcmp(file_ext,'.cxd'))
            
            if ~vidMetadata.force_user_vals
                coupler = vidMetadata.rawMetadata.get('Global dZoom');
                if isempty(coupler)
                    coupler = vidMetadata.rawMetadata.get('Global dRelayLensZoom');
                end
                if ~isempty(coupler)
                    vidMetadata.CouplerRatio = coupler;
                end
            end
        end

        if ~isfield(vidMetadata,'CouplerRatio')
            vidMetadata.CouplerRatio = 1;
        end
    end
    
    function get_objective
        % get objective value
        if (strcmp(file_ext,'.nd2') || strcmp(file_ext,'.cxd'))
            if ~vidMetadata.force_user_vals
                objective = vidMetadata.rawMetadata.get('Global wsObjectiveName');
                if ~isempty(objective)
                    objective = char(regexp(objective,'[0-9]*x','match'));
                    objective = str2double(char(regexp(objective,'[0-9]*','match')));
                    vidMetadata.Objective = objective;
                end
            end
        end
    end
    
    function get_binning
        % get binning value
        if (strcmp(file_ext,'.nd2') || strcmp(file_ext,'.cxd'))
            if ~vidMetadata.force_user_vals
                binning = vidMetadata.rawMetadata.get('Global Binning');
                if isempty(binning)
                    binning = vidMetadata.rawMetadata.get('Global Binning #1');
                end
                if isempty(binning)
                    binning = vidMetadata.rawMetadata.get('Global dBinningX');
                end
                if ~isempty(binning)
                    vidMetadata.Binning = binning;
                end
            end
        end
        if ~isempty(regexp(vidMetadata.Binning,'x','match'))
            vidMetadata.Binning = char(regexp(vidMetadata.Binning,'[0-9]*x','match'));
            vidMetadata.Binning = str2double(char(regexp(vidMetadata.Binning,'[0-9]*','match')));
        end
    end
    
    function get_timestamps
        % get timestamps of each frame
        if (strcmp(file_ext,'.nd2') || strcmp(file_ext,'.cxd'))
            rawMetadata_keys = char(vidMetadata.rawMetadata.keySet.toString);
            
            timestamps = regexp(rawMetadata_keys,'timestamp #\d*','match');
            if isempty(timestamps)
                timestamps = regexp(rawMetadata_keys,'Global Field \d* Time_From_Start','match');
            end
            vidMetadata.Time = zeros(vidMetadata.nFrames,1);
            for kt = 1:length(timestamps)
                index = str2double(char(regexp(timestamps{kt},'\d*','match')));
                vidMetadata.Time(index) = vidMetadata.rawMetadata.get(timestamps{kt});
            end
            if length(vidMetadata.Time) >= 2 % correction for some video files
                if vidMetadata.Time(end) < vidMetadata.Time(end-1)
                    vidMetadata.Time(end) = vidMetadata.Time(end-1) + mean(diff(vidMetadata.Time(1:end-1)));
                end
            end
            vidMetadata.Time = vidMetadata.Time - vidMetadata.Time(1);

            if nnz(vidMetadata.Time == 0) > 2
                % time is messed up, create it from the user input FrameRate
                vidMetadata.Time = linspace(0,vidMetadata.nFrames/vidMetadata.FrameRate,vidMetadata.nFrames)';
            else
                vidMetadata.FrameRate = 1/mean(diff(vidMetadata.Time));
            end
        end
        % create time stamp from framerate
        if ~isfield(vidMetadata,'Time')
            vidMetadata.Time = linspace(0,vidMetadata.nFrames/vidMetadata.FrameRate,vidMetadata.nFrames)';
        end
    end
    
    function get_calibration
        % get calibration value
        if strcmp(vidMetadata.CameraName,'Nikon A1plus')
            % confocal microscopes can zoom in while retaining number of pixels
            % the equation below does not account for this
            vidMetadata.Calibration = vidMetadata.rawMetadata.get('Global dCalibration');
        else
            vidMetadata.Calibration = vidMetadata.CameraPixelSize*vidMetadata.Binning/vidMetadata.Objective/vidMetadata.CouplerRatio/vidMetadata.MagMultiplier;
            % units are um/pixel
        end
    end

    function display_metadata
        if vidMetadata.verbose == 0
            return
        elseif vidMetadata.verbose == 1
            % display video metadata in command window
            fprintf('Video Information:\n')
            fields = fieldnames(vidMetadata);
            for i = 1:length(fields)
                if ~strcmp(fields{i},'rawMetadata')
                    if strcmp(fields{i},'Frames')
                        fprintf('%20s\t%g-%g\n',fields{i},vidMetadata.(fields{i})(1),vidMetadata.(fields{i})(end))
                    elseif strcmp(fields{i},'Time')
                        fprintf('%20s\t%g\n',fields{i},vidMetadata.(fields{i})(end))
                    else
                        if ischar(vidMetadata.(fields{i}))
                            fprintf('%20s\t%s\n',fields{i},vidMetadata.(fields{i}))
                        else
                            fprintf('%20s\t%g\n',fields{i},vidMetadata.(fields{i}))
                        end
                    end
                end
            end
        elseif vidMetadata.verbose == 2
            fig_analysis = figure('units','normalized','position',[0.3 0.1 0.4 0.7]);
            
            txt_title = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[0 0.9 1 0.1],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
                'String','Video Data');
            
            st = 0.02;
            wt = 0.4;
            sv = wt + 0.04;
            wv = 0.2;
            su = wt + wv + 0.04;
            wu = 0.2;

            txt_frames = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.85 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Frames:');
            val_frames = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[sv 0.85 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String',sprintf('%i',vidMetadata.nFrames));
            
            txt_framerate = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.8 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Framerate:');
            val_framerate = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.8 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.FrameRate,16),...
                'Callback',@callback_framerate);
            txt_framerate_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.8 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','fps');
            
            txt_time = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.75 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Time:');
            val_time = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.75 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.Time(end),16),...
                'Callback',@callback_time);
            txt_time_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.75 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','sec');
            


            txt_camerapixels = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.65 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Physical Pixel Size:');
            val_camerapixels = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.65 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.CameraPixelSize,16),'Callback',...
                @callback_calibration);
            txt_camerapixels_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.65 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','um');
            
            txt_objective = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.6 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Objective:');
            val_objective = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.6 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.Objective),'Callback',...
                @callback_calibration);
            txt_objective_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.6 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','x');

            txt_binning = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.55 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Binning:');
            val_binning = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.55 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.Binning),...
                'Callback', ...
                @callback_calibration);
            txt_binning_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.55 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String',['x' num2str(vidMetadata.Binning)]);

            txt_coupler = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.5 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Coupler:');
            val_coupler = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.5 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.CouplerRatio),'Callback',...
                @callback_calibration);
            txt_coupler_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.5 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','x');

            txt_magmultiplier = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.45 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Mag Multiplier:');
            val_magmultiplier = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.45 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.MagMultiplier),'Callback',...
                @callback_calibration);
            txt_magmultiplier = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.45 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','x');
            
            txt_calibration = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.4 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Calibration:');
            val_calibration = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[sv 0.4 wv 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.Calibration,16));
            txt_calibration_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.4 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','um/px');
            




            txt_dotsize = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.3 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Dot Size:');
            val_dotsize = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.3 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.DotSize,16));
            txt_dotsize_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.3 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','um');

            txt_dotspacing = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.25 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Dot Spacing:');
            val_dotspacing = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.25 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.DotSpacing,16));
            txt_dotspacing_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.25 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','um');

            txt_youngs = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.2 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Young''s Modulus:');
            val_youngs = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.2 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.YoungsModulus,16));
            txt_youngs_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.2 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','Pa');

            txt_pdmspct = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.15 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','PDMS %:');
            val_pdmspct = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.15 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','--', 'Callback', @callback_pdms);
            txt_pdmspct_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.15 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','%');

            txt_poisson = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.1 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Poisson''s Ratio:');
            val_poisson = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.1 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(vidMetadata.Poisson,16));

            
            



            but_done = uicontrol('parent',fig_analysis,'units','normalized','Style','pushbutton',...
                'Position',[0.25 0 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
                'String','DONE','Callback','uiresume');
            
%             but_reset = uicontrol('parent',fig_analysis,'units','normalized','Style','pushbutton',...
%                 'Position',[0.55 0.1 0.2 0.05],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
%                 'String','Reset',...
%                 'Callback',['val_camerapixels.String = num2str(vM.CameraPixelSize,16);',...
%                 'val_objective.String = num2str(vM.Objective);',...
%                 'val_binning.String = num2str(vM.Binning);',...
%                 'txt_binning_units.String = [''x'' num2str(vM.Binning)];',...
%                 'val_calibration.String = num2str(vM.Calibration,16);',...
%                 'val_framerate.String = num2str(vM.FrameRate,16);',...
%                 'val_time.String = num2str(vM.Time(end),16);']);
            
            check_time = false; % can change to true by callback above
            uiwait
            
            vidMetadata.CameraPixelSize = str2double(val_camerapixels.String);
            vidMetadata.Objective = str2double(val_objective.String);
            vidMetadata.CouplerRatio = str2double(val_coupler.String);
            vidMetadata.MagMultiplier = str2double(val_magmultiplier.String);
            vidMetadata.Binning = str2double(val_binning.String);
            vidMetadata.Calibration = str2double(val_calibration.String);
            vidMetadata.FrameRate = str2double(val_framerate.String);
            vidMetadata.DotSize = str2double(val_dotsize.String);
            vidMetadata.DotSpacing = str2double(val_dotspacing.String);
            vidMetadata.YoungsModulus = str2double(val_youngs.String);
            vidMetadata.Poisson = str2double(val_poisson.String);
            if check_time
                vidMetadata.Time = linspace(0,str2double(val_time.String),vM.nFrames)';
            end
            
            close(fig_analysis)

        end
        
        
        function callback_pdms(src, event)
            % palchesko:
            %   <20% y = 0.3236*x^2 + 2.0606*x + 5
            %   >20% y = 18.5x - 156.87
            %   x = [0, 1/11, 1/6, 0.5, 5/6, 1]*100
            %   y = [5, 50, 130, 830, 1340, 1720]

            % our data:
            %   x = [0, 2.5, 5, 10]
            %   y = [5, 7.7, 13.5, 46.7]
            % cfit: y = 0.4125*x^2 + 5
            
            x = [0, 2.5, 5, 10];
            y = [5, 7.7, 13.5, 46.7]*1000;

            val_youngs.String = num2str(interp1(x, y, str2double(val_pdmspct.String)),16);

        end
        function callback_framerate(src, event)
            val_time.String = num2str(vidMetadata.nFrames/str2double(val_framerate.String),16);
            check_time = true;
        end
        function callback_time(src, event)
            val_framerate.String = num2str(vidMetadata.nFrames/str2double(val_time.String),16);
            check_time = true;
        end
        function callback_calibration(src, event)
            txt_binning_units.String = ['x' val_binning.String];
            val_calibration.String = num2str(str2double(val_camerapixels.String)*str2double(val_binning.String)/str2double(val_objective.String)/str2double(val_coupler.String)/str2double(val_magmultiplier.String),16);
        end
    end
end