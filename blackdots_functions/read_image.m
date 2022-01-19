function [BDFrame,BDMetadata,CBFrame] = read_image(file_raw,varargin)
    [path_name, file_name, file_ext] = fileparts(file_raw);
    
    if ~isempty(varargin)
        userMetadata = varargin(0);
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
    end

end