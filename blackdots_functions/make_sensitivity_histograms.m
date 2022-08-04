%% load the data
clearvars; clc;

bd_files = dir('./0pct_1/*_blackdots_analysis'); % replace this with whichever cell you want to analyze

all_disp_in = {};
all_disp_out = {};
% for f = 1
nc = 0;
for f = 1:length(bd_files)
    latest = dir([bd_files(f).folder filesep bd_files(f).name]);
	if length(latest) > 2
        latest = latest(end);
    else
        continue
    end
    latest_file = dir([latest.folder filesep latest.name]);
    if length(latest_file) > 2
        latest_file = latest_file(end);
    else
        continue
    end
    
    disp(bd_files(f).name)
    
    load([latest_file.folder filesep latest_file.name])
    for ic = 1:length(celldata)
        if ~isempty(celldata(ic).celldots)
            nc = nc + 1;
            Xdisp_all = celldata(ic).Xdisp_k;
            Ydisp_all = celldata(ic).Ydisp_k;
            Xdisp_out = celldata(ic).Xdisp_k(celldata(ic).real_points(:) & ~celldata(ic).celldots(:));
            Ydisp_out = celldata(ic).Ydisp_k(celldata(ic).real_points(:) & ~celldata(ic).celldots(:));
            Xdisp_in = celldata(ic).Xdisp_k(celldata(ic).real_points(:) & celldata(ic).celldots(:));
            Ydisp_in = celldata(ic).Ydisp_k(celldata(ic).real_points(:) & celldata(ic).celldots(:));
            
            disp_all = sqrt(Xdisp_all.^2 + Ydisp_all.^2)*celldata(ic).vM.Calibration;
            disp_out = sqrt(Xdisp_out.^2 + Ydisp_out.^2)*celldata(ic).vM.Calibration;
            disp_in = sqrt(Xdisp_in.^2 + Ydisp_in.^2)*celldata(ic).vM.Calibration;
            
%             dist_out_boot = bootstrp(2000,@median,disp_out);
%             noise_limit = quantile(disp_out,0.95); % this is the level that includes the lower 95% of displacements outside the cell;
            
            all_disp_in{nc} = disp_in;
            all_disp_out{nc} = disp_out;
            
        end
    end
end

all_disp_in2 = cell2mat(all_disp_in');
all_disp_out2 = cell2mat(all_disp_out');

noise_limit = quantile(all_disp_out2,0.95)

%% make the histogram

nbins = 100;

% bin_edges = linspace(min([all_disp_in2; all_disp_out2]),max([all_disp_in2; all_disp_out2]),nbins);
bin_edges = bin_edges_0pct;

in_counts = histcounts(all_disp_in2,bin_edges);
out_counts = histcounts(all_disp_out2,bin_edges);

% in_counts_norm = in_counts/max(in_counts);
% out_counts_norm = out_counts/max(out_counts);

in_counts_norm = in_counts/sum(in_counts);
out_counts_norm = out_counts/sum(out_counts);

fig_1 = figure;
ax_1 = axes;
hold on
hout = histogram('BinCounts',out_counts_norm,'BinEdges',bin_edges,'EdgeColor','none');
hin = histogram('BinCounts',in_counts_norm,'BinEdges',bin_edges,'EdgeColor','none');
hold off
ylabel('Count (Normalized)')
xlabel('Displacement (\mum)')
box on
legend('No Cell','Outside Cell','Inside Cell')
