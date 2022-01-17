function [vidFrames,vidMetadata] = read_image(file_raw,userMetadata)
    if isfield(userMetadata,'Frames')
        frames_to_read = userMetadata.Frames;
    elseif isfield(userMetadata,'BDchannel')
        frames_to_read = userMetadata.BDchannel;
    else
        frames_to_read = [];
    end
    
    [path_name, file_name, file_ext] = fileparts(file_raw);
    
    vidMetadata = userMetadata;
    vidMetadata.PathName = path_name;
    vidMetadata.FileName = [file_name, file_ext];

    if ~isfield(vidMetadata,'MagMultiplier')
        vidMetadata.MagMultiplier = 1;
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
                if ~isempty(frames_to_read)
                    rawdata = bfopen_kb2(file_raw,frames_to_read);
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
            if ~isempty(frames_to_read)
                for k = 1:length(frames_to_read)
                    vidFrames(:,:,k) = imread(file_raw,frames_to_read(k));
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
            if ~isempty(frames_to_read)
                i = 1;
                k = 0;
                while hasFrame(v)
                    k = k + 1;
                    temp = readFrame(v);
                    if frames_to_read(i) == k
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
        % display video metadata
        fprintf('Video Information:\n')
        fields = fieldnames(vidMetadata);
        for i = 1:length(fields)
%             if length(getfield(vidMetadata,fields{i})) <= 1
            if ~strcmp(fields{i},'rawMetadata') && length(vidMetadata.(fields{i})) <= 1
%                 fprintf('%20s\t%g\n',fields{i},getfield(vidMetadata,fields{i}))
                fprintf('%20s\t%g\n',fields{i},vidMetadata.(fields{i}))
            end
        end
    end

end