%% analyze transient waveform and output results to command window

% controls for how a beat "start" and base is defined
start_pct = 0.03;
base_pct = 0.1;

ic = 1; % which cell to analyze

% code is borrowed from calcium analysis, so that is why everything is
% named 'ca_'
ca_rel = celldata(ic).total_force'*10^-6;

caI_peaks_init = peakfinder(ca_rel); % see end of script for peakfinder
if caI_peaks_init(end) == meta_BD.nFrames
    caI_peaks_init = caI_peaks_init(1:end-1);
end

if length(caI_peaks_init) == 1
    [~,caI_troughs(1)] = min(ca_rel(1:caI_peaks_init(1)));
    [~,caI_troughs(2)] = min(ca_rel(caI_peaks_init(1):end));
    
    caI_troughs(2) = caI_troughs(2) + caI_peaks_init(1) - 1;
    
    caI_peaks = caI_peaks_init;
elseif length(caI_peaks_init) == 2
    [~,caI_troughs(1)] = min(ca_rel(1:caI_peaks_init(1)));
    [~,caI_troughs(2)] = min(ca_rel(caI_peaks_init(1):caI_peaks_init(2)));
    [~,caI_troughs(3)] = min(ca_rel(caI_peaks_init(2):end));
    
    caI_troughs(2) = caI_troughs(2) + caI_peaks_init(1) - 1;
    caI_troughs(3) = caI_troughs(3) + caI_peaks_init(2) - 1;
    
    caI_peaks = caI_peaks_init;
else
    caI_troughs = zeros(length(caI_peaks_init),1);
    for np = 1:length(caI_peaks_init)-1
        [~,caI_troughs(np,1)] = min(ca_rel(caI_peaks_init(np):caI_peaks_init(np + 1)));
        caI_troughs(np,1) = caI_troughs(np) + caI_peaks_init(np) - 1;
    end
    np = np + 1;
    [~,caI_troughs(np,1)] = min(ca_rel(caI_peaks_init(np):end));
    caI_troughs(np,1) = caI_troughs(np) + caI_peaks_init(np) - 1;

    if (caI_peaks_init(1) >= 0.9*mean(caI_peaks_init(2:end-1) - caI_troughs(1:end-2)))
        [~,first_trough] = min(ca_rel(1:caI_peaks_init(1)));

        if first_trough == 1 %  get rid of first beat if trough is first point
            caI_peaks = caI_peaks_init(2:end);
        else
            caI_peaks = caI_peaks_init(1:end);
            caI_troughs = [first_trough; caI_troughs];
        end
    else % remove first beat if it is too close to starting time
        caI_peaks = caI_peaks_init(2:end);
    end

    if ((meta_BD.nFrames - caI_peaks(end)) < 0.9*mean(caI_troughs(2:end-1) - caI_peaks(1:end-1)))
        % remove last beat if it is too close to ending time
        caI_peaks = caI_peaks(1:end-1);
        caI_troughs = caI_troughs(1:end-1);
    end
end

for np = 1:length(caI_peaks)
    i_L = caI_troughs(np);
    i_R = caI_troughs(np + 1);
    temp_ca_rel = ca_rel(i_L:caI_peaks(np));
    caV_baselines_pre(np,1) = median(temp_ca_rel(temp_ca_rel < (base_pct*(max(temp_ca_rel) - min(temp_ca_rel)) + min(temp_ca_rel))));
    temp_ca_rel = ca_rel(caI_peaks(np):i_R);
    caV_baselines_post(np,1) = median(temp_ca_rel(temp_ca_rel < (base_pct*(max(temp_ca_rel) - min(temp_ca_rel)) + min(temp_ca_rel))));
    
    caI_start(np,1) = find(ca_rel(i_L:caI_peaks(np)) < (start_pct*(ca_rel(caI_peaks(np)) - caV_baselines_pre(np)) + caV_baselines_pre(np)),1,'Last');
    caI_start(np,1) = caI_start(np) + i_L - 1;
    
    if (caI_start(np) - i_L) < 0.2*mean(diff(caI_peaks))
        caV_baselines_pre(np,1) = ca_rel(i_L);
        caI_start(np,1) = find(ca_rel(i_L:caI_peaks(np)) < (start_pct*(ca_rel(caI_peaks(np)) - caV_baselines_pre(np)) + caV_baselines_pre(np)),1,'Last');
        caI_start(np,1) = caI_start(np) + i_L - 1;
        if (np > 1)
            caV_baselines_post(np-1,1) = ca_rel(i_L);
        end
    end
    
    if np == length(caI_peaks)
        if all(caV_baselines_post(1:length(caI_peaks)-1) == ca_rel(caI_troughs(2:length(caI_peaks-1)))) % if all previous beat starts were on a trough
            caV_baselines_post(np,1) = ca_rel(i_R);
        end
    end
% %     while ca_rel(caI_start(np) - 1) < ca_rel(caI_start(np))
% %         caI_start(np) = caI_start(np) - 1;
% %     end
    
% %     caI_start(np,1) = find((ca_rel(i_L:caI_peaks(np)) - ca_rel(i_L)) < start_pct*(ca_rel(caI_peaks(np)) - ca_rel(i_L)),1,'last');
%     caI_start(np,1) = find((ca_rel(i_L:caI_peaks(np)) - linspace(ca_rel(i_L),(ca_rel(i_R) - ca_rel(i_L))*(caI_peaks(np) - i_L)/(i_R - i_L) + ca_rel(i_L),(caI_peaks(np) - i_L + 1))') < start_pct*(ca_rel(caI_peaks(np)) - ca_rel(i_L)),1,'last');
    
    caI_50R(np,1) = find((ca_rel(caI_peaks(np):i_R) - ca_rel(i_R)) < 0.5*(ca_rel(caI_peaks(np)) - ca_rel(i_R)),1,'first');
    caI_50R(np,1) = caI_50R(np) + caI_peaks(np) - 1;
    caT_50R_i(np,1) = interp1(ca_rel([caI_50R(np)-1, caI_50R(np)]) - ca_rel(i_R),meta_BD.Time([caI_50R(np)-1, caI_50R(np)]),0.5*(ca_rel(caI_peaks(np)) - ca_rel(i_R)));
    
    caI_90R(np,1) = find((ca_rel(caI_peaks(np):i_R) - ca_rel(i_R)) < 0.1*(ca_rel(caI_peaks(np)) - ca_rel(i_R)),1,'first');
    caI_90R(np,1) = caI_90R(np) + caI_peaks(np) - 1;
    caT_90R_i(np,1) = interp1(ca_rel([caI_90R(np)-1, caI_90R(np)]) - ca_rel(i_R),meta_BD.Time([caI_90R(np)-1, caI_90R(np)]),0.1*(ca_rel(caI_peaks(np)) - ca_rel(i_R)));
    
    caI_10R(np,1) = find((ca_rel(caI_peaks(np):i_R) - ca_rel(i_R)) < 0.9*(ca_rel(caI_peaks(np)) - ca_rel(i_R)),1,'first');
    caI_10R(np,1) = caI_10R(np) + caI_peaks(np) - 1;
    caT_10R_i(np,1) = interp1(ca_rel([caI_10R(np)-1, caI_10R(np)]) - ca_rel(i_R),meta_BD.Time([caI_10R(np)-1, caI_10R(np)]),0.9*(ca_rel(caI_peaks(np)) - ca_rel(i_R)));
    
%     decay_c = lsqcurvefit(@(a,x) a(1)+a(2)*exp(-a(3)*x),[100 100 10],meta_BD.Time(caI_50R(np):caI_90R(np))-meta_BD.Time(caI_50R(np)),ca_rel(caI_50R(np):caI_90R(np)),[],[],...
%         optimoptions('lsqcurvefit','Display','off','MaxFunctionEvaluations',1000,'OptimalityTolerance',1e-6));
%     decay_fits(np,:) = decay_c;
%     tau(np,1) = 1/decay_c(3);
    
    % accurate peak fitting
    caI_80A(np,1) = find((ca_rel(i_L:caI_peaks(np)) - ca_rel(i_L)) < 0.8*(ca_rel(caI_peaks(np)) - ca_rel(i_L)),1,'last');
    caI_80A(np,1) = caI_80A(np) + i_L - 1;
    caI_20R(np,1) = find((ca_rel(caI_peaks(np):i_R) - ca_rel(i_R)) < 0.8*(ca_rel(caI_peaks(np)) - ca_rel(i_R)),1,'first');
    caI_20R(np,1) = caI_20R(np) + caI_peaks(np) - 1;

    peak_c = polyfit(meta_BD.Time(caI_80A(np):caI_20R(np)),ca_rel(caI_80A(np):caI_20R(np)),3);
    peak_fits(np,:) = peak_c;
    roots_peaks = roots(polyder(peak_c));
    caT_peaks_i(np,1) = roots_peaks(polyval(polyder(polyder(peak_c)),roots_peaks) < 0);
    caV_peaks_i(np,1) = polyval(peak_c,caT_peaks_i(np));
    
    % max up and downstroke
    [~,caI_RMaxUp(np,1)] = max(diff(ca_rel(caI_start(np):caI_peaks(np))));
    caI_RMaxUp(np,1) = caI_RMaxUp(np) + caI_start(np) - 1;
    RMaxUp(np,1) = (ca_rel(caI_RMaxUp(np) + 1) - ca_rel(caI_RMaxUp(np)))/(meta_BD.Time(caI_RMaxUp(np) + 1) - meta_BD.Time(caI_RMaxUp(np)));
    
    [~,caI_RMaxDn(np,1)] = min(diff(ca_rel(caI_peaks(np):caI_90R(np))));
    caI_RMaxDn(np,1) = caI_RMaxDn(np) + caI_peaks(np) - 1;
    RMaxDn(np,1) = (ca_rel(caI_RMaxDn(np) + 1) - ca_rel(caI_RMaxDn(np)))/(meta_BD.Time(caI_RMaxDn(np) + 1) - meta_BD.Time(caI_RMaxDn(np)));
end

caT_start = meta_BD.Time(caI_start);
caT_peaks = meta_BD.Time(caI_peaks);

caV_start = ca_rel(caI_start);
caV_peaks = ca_rel(caI_peaks);
caV_50R_i = (0.5*(ca_rel(caI_peaks) - ca_rel(caI_troughs(2:end))) + ca_rel(caI_troughs(2:end)));
caV_90R_i = (0.1*(ca_rel(caI_peaks) - ca_rel(caI_troughs(2:end))) + ca_rel(caI_troughs(2:end)));
caV_10R_i = (0.9*(ca_rel(caI_peaks) - ca_rel(caI_troughs(2:end))) + ca_rel(caI_troughs(2:end)));

caR_peaks = (caV_peaks - caV_start)./(caT_peaks - caT_start);
caR_50R_i = (caV_peaks - caV_50R_i)./(caT_50R_i - caT_peaks);
caR_90R_i = (caV_peaks - caV_90R_i)./(caT_90R_i - caT_peaks);
caR_10R_i = (caV_peaks - caV_10R_i)./(caT_10R_i - caT_peaks);

Tpeak_mean = mean(caT_peaks - caT_start);
T50R_mean = mean(caT_50R_i - caT_peaks);
T90R_mean = mean(caT_90R_i - caT_peaks);
T10R_mean = mean(caT_10R_i - caT_peaks);

Rpeak_mean = mean(caR_peaks);
R50R_mean = mean(caR_50R_i);
R90R_mean = mean(caR_90R_i);
R10R_mean = mean(caR_10R_i);

freq_mean = (length(caT_peaks) - 1)./(caT_peaks(end)-caT_peaks(1));

% tau_mean = mean(tau);

RMaxUp_mean = mean(RMaxUp);
RMaxDn_mean = mean(RMaxDn);

fprintf('Success: Waveform analysis completed.\n')
fprintf('\tAverage beat frequency = %0.4f Hz.\n',freq_mean)

% figure
% hold on
% plot(meta_BD.Time,ca_rel,'-k','linewidth',1);
% plot(caT_start,caV_start,'>','markeredgecolor','none','markerfacecolor',[0 0.75 0])
% plot(caT_peaks,caV_peaks,'^','markeredgecolor','none','markerfacecolor','r')
% plot(caT_50R_i,caV_50R_i,'v','markeredgecolor','none','markerfacecolor',[0 0 1])
% plot(caT_90R_i,caV_90R_i,'v','markeredgecolor','none','markerfacecolor',[0.6 0 1])
% hold off

figure('Position',[700,100,600,200])
hold on
plot(meta_BD.Time,ca_rel,'-k','linewidth',2);
plot(caT_start,caV_start,'>','markeredgecolor','none','markerfacecolor',[0 0.75 0])
plot(caT_peaks,caV_peaks,'^','markeredgecolor','none','markerfacecolor','r')
plot(caT_50R_i,caV_50R_i,'v','markeredgecolor','none','markerfacecolor',[0 0 1])
plot(caT_90R_i,caV_90R_i,'v','markeredgecolor','none','markerfacecolor',[0.6 0 1])
% plot(meta_BD.Time,10^-6*celldata(ic).total_force,'-k','linewidth',2)
hold off
xlabel('Time [s]')
ylabel('Total Force [uN]')
box off
set(gca,'linewidth',1.5,'tickdir','out','XColor','k','YColor','k')

print(sprintf('plot_force_transient_labeled_cell%i', ic),'-dpng')
print(sprintf('plot_force_transient_labeled_cell%i', ic),'-dsvg','-vector')

fprintf(['Results\n',...
    '%s\n',...
    '%s\n',...
    '%12s\t%12s\t',...
    '%12s\t%12s\t%12s\t',...
    '%12s\t%12s\t%12s\t',...
    '%12s\t%12s\t%12s\n',...
    '%12.2f\t%12.2f\t',...
    '%12.2f\t%12.2f\t%12.2f\t',...
    '%12.3f\t%12.3f\t%12.3f\t',...
    '%12.3f\t%12.3f\t%12.3f\n'],...
    'File',...
    [meta_BD.PathName filesep meta_BD.FileName], ...
    'Area(um2)','Freq(Hz)',...
    'F0(uN)','Fmax(uN)','Fmax/F0',...
    'Tpeak(s)','T50R(s)','T90R(s)',...
    'Rpeak(uN/s)','R50R(uN/s)','R90R(uN/s)',...
    celldata(ic).area, freq_mean,...
    mean(caV_baselines_pre), mean(caV_peaks), mean(caV_peaks./caV_baselines_pre),...
    Tpeak_mean, T50R_mean, T90R_mean,...
    Rpeak_mean, R50R_mean, R90R_mean)

%% code from the internet
function varargout = peakfinder(x0, sel, thresh, extrema)
    %PEAKFINDER Noise tolerant fast peak finding algorithm
    %   INPUTS:
    %       x0 - A real vector from the maxima will be found (required)
    %       sel - The amount above surrounding data for a peak to be
    %           identified (default = (max(x0)-min(x0))/4). Larger values mean
    %           the algorithm is more selective in finding peaks.
    %       thresh - A threshold value which peaks must be larger than to be
    %           maxima or smaller than to be minima.
    %       extrema - 1 if maxima are desired, -1 if minima are desired
    %           (default = maxima, 1)
    %   OUTPUTS:
    %       peakLoc - The indicies of the identified peaks in x0
    %       peakMag - The magnitude of the identified peaks
    %
    %   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
    %       are at least 1/4 the range of the data above surrounding data.
    %
    %   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
    %       that are at least sel above surrounding data.
    %
    %   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local 
    %       maxima that are at least sel above surrounding data and larger
    %       (smaller) than thresh if you are finding maxima (minima).
    %
    %   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
    %       data if extrema > 0 and the minima of the data if extrema < 0
    %
    %   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
    %       local maxima as well as the magnitudes of those maxima
    %
    %   If called with no output the identified maxima will be plotted along
    %       with the input data.
    %
    %   Note: If repeated values are found the first is identified as the peak
    %
    % Ex:
    % t = 0:.0001:10;
    % x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
    % x(1250:1255) = max(x);
    % peakfinder(x)
    %
    % Copyright Nathanael C. Yoder 2011 (nyoder@gmail.com)

    % Perform error checking and set defaults if not passed in
    narginchk(1,4);
    nargoutchk(0,2);

    s = size(x0);
    flipData =  s(1) < s(2);
    len0 = numel(x0);
    if len0 ~= s(1) && len0 ~= s(2)
        error('PEAKFINDER:Input','The input data must be a vector')
    elseif isempty(x0)
        varargout = {[],[]};
        return;
    end
    if ~isreal(x0)
        warning('PEAKFINDER:NotReal','Absolute value of data will be used')
        x0 = abs(x0);
    end

    if nargin < 2 || isempty(sel)
        sel = (max(x0)-min(x0))/4;
    elseif ~isnumeric(sel) || ~isreal(sel)
        sel = (max(x0)-min(x0))/4;
        warning('PEAKFINDER:InvalidSel',...
            'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
    elseif numel(sel) > 1
        warning('PEAKFINDER:InvalidSel',...
            'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
        sel = sel(1);
    end

    if nargin < 3 || isempty(thresh)
        thresh = [];
    elseif ~isnumeric(thresh) || ~isreal(thresh)
        thresh = [];
        warning('PEAKFINDER:InvalidThreshold',...
            'The threshold must be a real scalar. No threshold will be used.')
    elseif numel(thresh) > 1
        thresh = thresh(1);
        warning('PEAKFINDER:InvalidThreshold',...
            'The threshold must be a scalar.  The first threshold value in the vector will be used.')
    end

    if nargin < 4 || isempty(extrema)
        extrema = 1;
    else
        extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
        if extrema == 0
            error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
        end
    end

    x0 = extrema*x0(:); % Make it so we are finding maxima regardless
    thresh = thresh*extrema; % Adjust threshold according to extrema.
    dx0 = diff(x0); % Find derivative
    dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
    ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign

    % Include endpoints in potential peaks and valleys
    x = [x0(1);x0(ind);x0(end)];
    ind = [1;ind;len0];

    % x only has the peaks, valleys, and endpoints
    len = numel(x);
    minMag = min(x);


    if len > 2 % Function with peaks and valleys

        % Set initial parameters for loop
        tempMag = minMag;
        foundPeak = false;
        leftMin = minMag;

        % Deal with first point a little differently since tacked it on
        % Calculate the sign of the derivative since we taked the first point
        %  on it does not neccessarily alternate like the rest.
        signDx = sign(diff(x(1:3)));
        if signDx(1) <= 0 % The first point is larger or equal to the second
            ii = 0;
            if signDx(1) == signDx(2) % Want alternating signs
                x(2) = [];
                ind(2) = [];
                len = len-1;
            end
        else % First point is smaller than the second
            ii = 1;
            if signDx(1) == signDx(2) % Want alternating signs
                x(1) = [];
                ind(1) = [];
                len = len-1;
            end
        end

        % Preallocate max number of maxima
        maxPeaks = ceil(len/2);
        peakLoc = zeros(maxPeaks,1);
        peakMag = zeros(maxPeaks,1);
        cInd = 1;
        % Loop through extrema which should be peaks and then valleys
        while ii < len
            ii = ii+1; % This is a peak
            % Reset peak finding if we had a peak and the next peak is bigger
            %   than the last or the left min was small enough to reset.
            if foundPeak
                tempMag = minMag;
                foundPeak = false;
            end

            % Make sure we don't iterate past the length of our vector
            if ii == len
                break; % We assign the last point differently out of the loop
            end

            % Found new peak that was lager than temp mag and selectivity larger
            %   than the minimum to its left.
            if x(ii) > tempMag && x(ii) > leftMin + sel
                tempLoc = ii;
                tempMag = x(ii);
            end

            ii = ii+1; % Move onto the valley
            % Come down at least sel from peak
            if ~foundPeak && tempMag > sel + x(ii)
                foundPeak = true; % We have found a peak
                leftMin = x(ii);
                peakLoc(cInd) = tempLoc; % Add peak to index
                peakMag(cInd) = tempMag;
                cInd = cInd+1;
            elseif x(ii) < leftMin % New left minima
                leftMin = x(ii);
            end
        end

        % Check end point
        if x(end) > tempMag && x(end) > leftMin + sel
            peakLoc(cInd) = len;
            peakMag(cInd) = x(end);
            cInd = cInd + 1;
        elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
            peakLoc(cInd) = tempLoc;
            peakMag(cInd) = tempMag;
            cInd = cInd + 1;
        end

        % Create output
        peakInds = ind(peakLoc(1:cInd-1));
        peakMags = peakMag(1:cInd-1);
    else % This is a monotone function where an endpoint is the only peak
        [peakMags,xInd] = max(x);
        if peakMags > minMag + sel
            peakInds = ind(xInd);
        else
            peakMags = [];
            peakInds = [];
        end
    end

    % Apply threshold value.  Since always finding maxima it will always be
    %   larger than the thresh.
    if ~isempty(thresh)
        m = peakMags>thresh;
        peakInds = peakInds(m);
        peakMags = peakMags(m);
    end



    % Rotate data if needed
    if flipData
        peakMags = peakMags.';
        peakInds = peakInds.';
    end



    % Change sign of data if was finding minima
    if extrema < 0
        peakMags = -peakMags;
        x0 = -x0;
    end
    % Plot if no output desired
    if nargout == 0
        if isempty(peakInds)
            disp('No significant peaks found')
        else
            figure;
            plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
        end
    else
        varargout = {peakInds,peakMags};
    end
end