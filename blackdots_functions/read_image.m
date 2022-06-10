function [BDFrame,BDMetadata,CBFrame] = read_image(file_raw,varargin)
    [path_name, file_name, file_ext] = fileparts(file_raw);
    
    if ~isempty(varargin)
        userMetadata = varargin{1};
    else
        userMetadata = struct();
        userMetadata.force_user_vals = false;
    end
    BDMetadata = userMetadata;
    BDMetadata.PathName = path_name;
    BDMetadata.FileName = [file_name, file_ext];

    if ~isfield(BDMetadata,'MagMultiplier')
        BDMetadata.MagMultiplier = 1;
    end
    
    fprintf('Reading file < %s >\n',[file_name, file_ext])
    
    load_image

    crop_image

    get_camera_pixel_size
    
    get_coupler_ratio
    
    get_objective
    
    get_binning
    
    get_calibration

    display_metadata

    %% functions listed above
    function load_image
        % load video frames
        if (strcmp(file_ext,'.nd2'))
            if exist('bfopen_kb2','file')
                if isfield(BDMetadata,'BDchannel')
                    rawdata = bfopen_kb2(file_raw,BDMetadata.BDchannel);
                else
                    rawdata = bfopen_kb2(file_raw,1);
                end
                
                if isfield(BDMetadata,'CBchannel')
                    rawdataCB = bfopen_kb2(file_raw,BDMetadata.CBchannel);
                end

            elseif exist('bfopen','file')
                rawdata = bfopen(file_raw);
            else
                error('CANNOT READ FILE: Please ensure you have the bioformats functions installed.')
            end
            rawFrames = rawdata{1};
            rawFrames = rawFrames(~cellfun(@isempty,rawFrames(:,1)),1);
            % get raw metadata
            BDMetadata.rawMetadata = rawdata{2};
            
            BDFrame = rawFrames{1,1};

            if exist('rawdataCB','var')
                CBFrame = rawdataCB{1}{1,1};
            else
                CBFrame = BDFrame;
            end
        
        elseif strcmp(file_ext,'.tif') || strcmp(file_ext,'.tiff')
            if isfield(BDMetadata,'BDchannel')
                BDFrame = imread(file_raw,BDMetadata.BDchannel);
            else
                BDFrame = imread(file_raw,1);
            end
        elseif strcmp(file_ext,'.jpg') || strcmp(file_ext,'.png') || strcmp(file_ext,'.bmp')
            BDFrame = imread(file_raw);
            if size(BDFrame,3) > 1
                BDFrame = rgb2gray(BDFrame);
            end
        else
            error('UNSUPPORTED IMAGE FORMAT: Please use .nd2, .tif, .jpg, .png, or .bmp.')
        end
        [BDMetadata.M,BDMetadata.N] = size(BDFrame);
        BDMetadata.uFrame = 1; % these are needed for the black dots analysis
        BDMetadata.nFrames = 1;
    end
    
    function crop_image
        if isfield(BDMetadata,'Crop') && all(BDMetadata.Crop ~= false)
            if BDMetadata.Crop == true
                figure(1)
                title('CROP IMAGE')
                [~,BDMetadata.Crop] = imcrop(mat2gray(BDFrame));
                close(1)
            end

            BDFrame = imcrop(BDFrame,BDMetadata.Crop);
            [BDMetadata.M,BDMetadata.N] = size(BDFrame);
        end
    end
    
    function get_camera_pixel_size
        % get physical camera pixel value
        if strcmp(file_ext,'.nd2')

            BDMetadata.CameraName = BDMetadata.rawMetadata.get('Global Camera Name');
            if isempty(BDMetadata.CameraName)
                BDMetadata.CameraName = BDMetadata.rawMetadata.get('Global CameraUserName');
            end
            if isempty(BDMetadata.CameraName)
                BDMetadata.CameraName = BDMetadata.rawMetadata.get('Global CameraName');
            end
            if isempty(BDMetadata.CameraName)
                BDMetadata.CameraName = BDMetadata.rawMetadata.get('Camera Name');
            end
            
            if ~BDMetadata.force_user_vals
                if strcmp(BDMetadata.CameraName,'Andor Clara DR-1357') % Andor
                    BDMetadata.CameraPixelSize = 6.45; % physical pixel size in userMetadata
                elseif strcmp(BDMetadata.CameraName,'C11440-10C') % Hamamatsu Flash2.8
                    BDMetadata.CameraPixelSize = 3.63; % physical pixel size in userMetadata
                elseif strcmp(BDMetadata.CameraName,'Flash4.0, SN:301638') % Hamamatsu Flash4.0 300017
                    BDMetadata.CameraPixelSize = 6.5; % physical pixel size in userMetadata
                elseif strcmp(BDMetadata.CameraName,'Nikon A1plus') % Garvey Confocal A1
                    BDMetadata.CameraPixelSize = 3.8529; % physical pixel size in userMetadata
                else
                    BDMetadata.CameraName = 'Unknown Camera';
                    % use user submitted CameraPixelSize
                end
            end
        else
            BDMetadata.CameraName = 'Unknown Camera';
            % use user submitted CameraPixelSize
        end
    end
    
    function get_coupler_ratio
        % get coupler ratio value
        if strcmp(file_ext,'.nd2')
            
            if ~BDMetadata.force_user_vals
                coupler = BDMetadata.rawMetadata.get('Global dZoom');
                if isempty(coupler)
                    coupler = BDMetadata.rawMetadata.get('Global dRelayLensZoom');
                end
                if ~isempty(coupler)
                    BDMetadata.CouplerRatio = coupler;
                end
            end
        end

        if ~isfield(BDMetadata,'CouplerRatio')
            BDMetadata.CouplerRatio = 1;
        end
    end
    
    function get_objective
        % get objective value
        if strcmp(file_ext,'.nd2')
            if ~BDMetadata.force_user_vals
                objective = BDMetadata.rawMetadata.get('Global wsObjectiveName');
                if ~isempty(objective)
                    objective = char(regexp(objective,'[0-9]*x','match'));
                    objective = str2double(char(regexp(objective,'[0-9]*','match')));
                    BDMetadata.Objective = objective;
                end
            end
        end
    end
    
    function get_binning
        % get binning value
        if strcmp(file_ext,'.nd2')
            if ~BDMetadata.force_user_vals
                binning = BDMetadata.rawMetadata.get('Global Binning');
                if isempty(binning)
                    binning = BDMetadata.rawMetadata.get('Global Binning #1');
                end
                if isempty(binning)
                    binning = BDMetadata.rawMetadata.get('Global dBinningX');
                end
                if ~isempty(binning)
                    BDMetadata.Binning = binning;
                end
            end
        end
        if ~isempty(regexp(BDMetadata.Binning,'x','match'))
            BDMetadata.Binning = char(regexp(BDMetadata.Binning,'[0-9]*x','match'));
            BDMetadata.Binning = str2double(char(regexp(BDMetadata.Binning,'[0-9]*','match')));
        end
    end
    
    function get_calibration
        % get calibration value
        if strcmp(BDMetadata.CameraName,'Nikon A1plus')
            % confocal microscopes can zoom in while retaining number of pixels
            % the equation below does not account for this
            BDMetadata.Calibration = BDMetadata.rawMetadata.get('Global dCalibration');
        else
            BDMetadata.Calibration = BDMetadata.CameraPixelSize*BDMetadata.Binning/BDMetadata.Objective/BDMetadata.CouplerRatio/BDMetadata.MagMultiplier;
            % units are um/pixel
        end
    end

    function display_metadata
        if BDMetadata.verbose == 0
            return
        elseif BDMetadata.verbose == 1
            % display video metadata
            fprintf('Video Information:\n')
            fields = fieldnames(BDMetadata);
            for i = 1:length(fields)
                if ~strcmp(fields{i},'rawMetadata')
                    if ischar(BDMetadata.(fields{i}))
                        fprintf('%20s\t%s\n',fields{i},BDMetadata.(fields{i}))
                    else
                        fprintf('%20s\t%g\n',fields{i},BDMetadata.(fields{i}))
                    end
                end
            end
        elseif BDMetadata.verbose == 2
            fig_analysis = figure('units','normalized','position',[0.3 0.1 0.4 0.7]);
            
            txt_title = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[0 0.9 1 0.1],'HorizontalAlignment','Center','FontUnits','normalized','FontSize',0.8,...
                'String','Image Data');
            
            st = 0.02;
            wt = 0.4;
            sv = wt + 0.04;
            wv = 0.2;
            su = wt + wv + 0.04;
            wu = 0.2;

%             ax1 = axes(fig_analysis,'units','normalized','position',[su 0.8 wu 0.05]);
%             im1 = imagesc(ax1, BDFrame);
% 
%             ax2 = axes(fig_analysis,'units','normalized','position',[su 0.75 wu 0.05]);
%             im2 = imagesc(ax1, CBFrame);

            txt_channelBD = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.8 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Black Dots Channel:');
            val_channelBD = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.8 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',3,...
                'Callback',@callback_channelBD);
            
            txt_channelCB = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.75 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Cell Boundary Channel:');
            val_channelCB = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.75 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',2,...
                'Callback',@callback_channelCB);

            txt_camerapixels = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.65 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Physical Pixel Size:');
            val_camerapixels = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.65 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(BDMetadata.CameraPixelSize,16),'Callback',...
                @callback_calibration);
            txt_camerapixels_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.65 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','um');
            
            txt_objective = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.6 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Objective:');
            val_objective = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.6 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(BDMetadata.Objective),'Callback',...
                @callback_calibration);
            txt_objective_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.6 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','x');

            txt_binning = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.55 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Binning:');
            val_binning = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.55 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(BDMetadata.Binning),...
                'Callback', ...
                @callback_calibration);
            txt_binning_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.55 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String',['x' num2str(BDMetadata.Binning)]);

            txt_coupler = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.5 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Coupler:');
            val_coupler = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.5 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(BDMetadata.CouplerRatio),'Callback',...
                @callback_calibration);
            txt_coupler_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.5 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','x');

            txt_magmultiplier = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.45 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Mag Multiplier:');
            val_magmultiplier = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.45 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(BDMetadata.MagMultiplier),'Callback',...
                @callback_calibration);
            txt_magmultiplier = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.45 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','x');
            
            txt_calibration = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.4 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Calibration:');
            val_calibration = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[sv 0.4 wv 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(BDMetadata.Calibration,16));
            txt_calibration_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.4 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','um/px');
            




            txt_dotsize = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.3 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Dot Size:');
            val_dotsize = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.3 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(BDMetadata.DotSize,16));
            txt_dotsize_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.3 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','um');

            txt_dotspacing = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.25 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Dot Spacing:');
            val_dotspacing = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.25 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(BDMetadata.DotSpacing,16));
            txt_dotspacing_units = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[su 0.25 wu 0.05],'HorizontalAlignment','Left','FontUnits','normalized','FontSize',0.8,...
                'String','um');

            txt_youngs = uicontrol('parent',fig_analysis,'units','normalized','Style','text',...
                'Position',[st 0.2 wt 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String','Young''s Modulus:');
            val_youngs = uicontrol('parent',fig_analysis,'units','normalized','Style','edit',...
                'Position',[sv 0.2 wv 0.05],'HorizontalAlignment','Right','FontUnits','normalized','FontSize',0.8,...
                'String',num2str(BDMetadata.YoungsModulus,16));
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
                'String',num2str(BDMetadata.Poisson,16));

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
            
            BDMetadata.BDchannel = str2double(val_channelBD.String);
            BDFrame = imread(file_raw,BDMetadata.BDchannel);
            BDMetadata.CBchannel = str2double(val_channelCB.String);
            CBFrame = imread(file_raw,BDMetadata.CBchannel);
            BDMetadata.CameraPixelSize = str2double(val_camerapixels.String);
            BDMetadata.Objective = str2double(val_objective.String);
            BDMetadata.CouplerRatio = str2double(val_coupler.String);
            BDMetadata.MagMultiplier = str2double(val_magmultiplier.String);
            BDMetadata.Binning = str2double(val_binning.String);
            BDMetadata.Calibration = str2double(val_calibration.String);
            BDMetadata.DotSize = str2double(val_dotsize.String);
            BDMetadata.DotSpacing = str2double(val_dotspacing.String);
            BDMetadata.YoungsModulus = str2double(val_youngs.String);
            BDMetadata.Poisson = str2double(val_poisson.String);
            
            close(fig_analysis)
        end
        
        function callback_channelBD(src, event)
%             BDFrame = imread(file_raw,str2double(val_channelBD.String));
%             set(im1,'CData',BDFrame)
            return
        end
        function callback_channelCB(src, event)
%             CBFrame = imread(file_raw,str2double(val_channelCB.String));
%             set(im2,'CData',CBFrame)
            return
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
        function callback_calibration(src, event)
            txt_binning_units.String = ['x' val_binning.String];
            val_calibration.String = num2str(str2double(val_camerapixels.String)*str2double(val_binning.String)/str2double(val_objective.String)/str2double(val_coupler.String)/str2double(val_magmultiplier.String),16);
        end
    end

end