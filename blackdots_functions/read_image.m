function [imFrame,imMetadata] = read_image(file_raw,userMetadata)
    [path_name, file_name, file_ext] = fileparts(file_raw);
    
    imMetadata = userMetadata;
    imMetadata.PathName = path_name;
    imMetadata.FileName = [file_name, file_ext];

    if ~isfield(imMetadata,'MagMultiplier')
        imMetadata.MagMultiplier = 1;
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
        if (strcmp(file_ext,'.nd2')
            if exist('bfopen_kb2','file')
                if isfield(imMetadata,'BDchannel')
                    rawdata = bfopen_kb2(file_raw,imMetadata.BDchannel);
                else
                    rawdata = bfopen_kb2(file_raw,1);
                end
            elseif exist('bfopen','file')
                rawdata = bfopen(file_raw);
            else
                error('CANNOT READ FILE: Please ensure you have the bioformats functions installed.')
            end
            rawFrames = rawdata{1};
            rawFrames = rawFrames(~cellfun(@isempty,rawFrames(:,1)),1);
            % get raw metadata
            imMetadata.rawMetadata = rawdata{2};
            
            imFrame = rawFrames{1,1};
        
        elseif strcmp(file_ext,'.tif') || strcmp(file_ext,'.tiff')
            if isfield(imMetadata,'BDchannel')
                imFrame = imread(file_raw,imMetadata.BDchannel);
            else
                imFrame = imread(file_raw,1);
            end
        elseif strcmp(file_ext,'.jpg') || strcmp(file_ext,'.png') || strcmp(file_ext,'.bmp')
            imFrame = imread(file_raw);
            if size(imFrame,3) > 1
                imFrame = rgb2gray(imFrame);
            end
        else
            error('UNSUPPORTED IMAGE FORMAT: Please use .nd2, .tif, .jpg, .png, or .bmp.')
        end
        [imMetadata.M,imMetadata.N] = size(imFrame);
    end
    
    function crop_image
        if isfield(imMetadata,'Crop') && all(imMetadata.Crop ~= false)
            if imMetadata.Crop == true
                figure(1)
                title('CROP IMAGE')
                [~,imMetadata.Crop] = imcrop(mat2gray(imFrame));
                close(1)
            end

            imFrame = imcrop(imFrame,imMetadata.Crop);
            [imMetadata.M,imMetadata.N] = size(imFrame);
        end
    end
    
    function get_camera_pixel_size
        % get physical camera pixel value
        if strcmp(file_ext,'.nd2')

            imMetadata.CameraName = imMetadata.rawMetadata.get('Global Camera Name');
            if isempty(imMetadata.CameraName)
                imMetadata.CameraName = imMetadata.rawMetadata.get('Global CameraUserName');
            end
            if isempty(imMetadata.CameraName)
                imMetadata.CameraName = imMetadata.rawMetadata.get('Global CameraName');
            end
            if isempty(imMetadata.CameraName)
                imMetadata.CameraName = imMetadata.rawMetadata.get('Camera Name');
            end
            
            if ~imMetadata.force_user_vals
                if strcmp(imMetadata.CameraName,'Andor Clara DR-1357') % Andor
                    imMetadata.CameraPixelSize = 6.45; % physical pixel size in userMetadata
                elseif strcmp(imMetadata.CameraName,'C11440-10C') % Hamamatsu Flash2.8
                    imMetadata.CameraPixelSize = 3.63; % physical pixel size in userMetadata
                elseif strcmp(imMetadata.CameraName,'Flash4.0, SN:301638') % Hamamatsu Flash4.0 300017
                    imMetadata.CameraPixelSize = 6.5; % physical pixel size in userMetadata
                elseif strcmp(imMetadata.CameraName,'Nikon A1plus') % Garvey Confocal A1
                    imMetadata.CameraPixelSize = 3.8529; % physical pixel size in userMetadata
                else
                    imMetadata.CameraName = 'Unknown Camera';
                    % use user submitted CameraPixelSize
                end
            end
        else
            imMetadata.CameraName = 'Unknown Camera';
            % use user submitted CameraPixelSize
        end
    end
    
    function get_coupler_ratio
        % get coupler ratio value
        if strcmp(file_ext,'.nd2')
            
            if ~imMetadata.force_user_vals
                coupler = imMetadata.rawMetadata.get('Global dZoom');
                if isempty(coupler)
                    coupler = imMetadata.rawMetadata.get('Global dRelayLensZoom');
                end
                if ~isempty(coupler)
                    imMetadata.CouplerRatio = coupler;
                end
            end
        end

        if ~isfield(imMetadata,'CouplerRatio')
            imMetadata.CouplerRatio = 1;
        end
    end
    
    function get_objective
        % get objective value
        if strcmp(file_ext,'.nd2')
            if ~imMetadata.force_user_vals
                objective = imMetadata.rawMetadata.get('Global wsObjectiveName');
                if ~isempty(objective)
                    objective = char(regexp(objective,'[0-9]*x','match'));
                    objective = str2double(char(regexp(objective,'[0-9]*','match')));
                    imMetadata.Objective = objective;
                end
            end
        end
    end
    
    function get_binning
        % get binning value
        if strcmp(file_ext,'.nd2')
            if ~imMetadata.force_user_vals
                binning = imMetadata.rawMetadata.get('Global Binning');
                if isempty(binning)
                    binning = imMetadata.rawMetadata.get('Global Binning #1');
                end
                if isempty(binning)
                    binning = imMetadata.rawMetadata.get('Global dBinningX');
                end
                if ~isempty(binning)
                    imMetadata.Binning = binning;
                end
            end
        end
        if ~isempty(regexp(imMetadata.Binning,'x','match'))
            imMetadata.Binning = char(regexp(imMetadata.Binning,'[0-9]*x','match'));
            imMetadata.Binning = str2double(char(regexp(imMetadata.Binning,'[0-9]*','match')));
        end
    end
    
    function get_calibration
        % get calibration value
        if strcmp(imMetadata.CameraName,'Nikon A1plus')
            % confocal microscopes can zoom in while retaining number of pixels
            % the equation below does not account for this
            imMetadata.Calibration = imMetadata.rawMetadata.get('Global dCalibration');
        else
            imMetadata.Calibration = imMetadata.CameraPixelSize*imMetadata.Binning/imMetadata.Objective/imMetadata.CouplerRatio/imMetadata.MagMultiplier;
            % units are um/pixel
        end
    end

    function display_metadata
        % display video metadata
        fprintf('Video Information:\n')
        fields = fieldnames(imMetadata);
        for i = 1:length(fields)
%             if length(getfield(vidMetadata,fields{i})) <= 1
            if ~strcmp(fields{i},'rawMetadata') && length(imMetadata.(fields{i})) <= 1
%                 fprintf('%20s\t%g\n',fields{i},getfield(vidMetadata,fields{i}))
                fprintf('%20s\t%g\n',fields{i},imMetadata.(fields{i}))
            end
        end
    end

end