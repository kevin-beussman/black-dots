function [rot_angle] = get_rot_from_pts(input,spacing)
px = input(:,1);
py = input(:,2);
sub = 1; % sub-divide pixels to get more accurate FFT after gaussian blur
dots_img = zeros(max(py),max(px));
for np = 1:size(px,1)
    dots_img(round(py(np)),round(px(np))) = 1;
end

% max_bound = max([vM.M vM.N]);
max_bound = size(dots_img);
dots_img_fft = fftshift(abs(fft2(dots_img,2^nextpow2(2*max_bound),2^nextpow2(2*max_bound))).^2);
mid_pt = [floor(size(dots_img_fft,1)/2)+1,floor(size(dots_img_fft,2)/2)+1];
% mid_val = dots_img_fft(mid_pt(1),mid_pt(2));
dots_img_fft(mid_pt(1),mid_pt(2)) = mean(mean(dots_img_fft(mid_pt(1)-1:mid_pt(1)+1,mid_pt(2)-1:mid_pt(2)+1).*[1 1 1;1 0 1;1 1 1]));
dots_img_fft_sub = repelem(dots_img_fft,sub,sub);
dots_img_fft_filt = imgaussfilt(dots_img_fft_sub,1);

[Mfft,Nfft] = size(dots_img_fft);

mid_pt_sub = [Mfft*sub/2 + 1 + 0.5*(sub - 1),Nfft*sub/2 + 1 + 0.5*(sub - 1)];
max_sp = 1.25;
min_sp = 0.75;
% dots_fft_mask = false(Mfft*sub,Nfft*sub); OLD METHOD
% for iy = 1:Mfft*sub
%     for jx = 1:Nfft*sub
%         if (1/sqrt(((mid_pt_sub(1)-iy)/(Mfft*sub*vM.Calibration))^2 + ((mid_pt_sub(2)-jx)/(Nfft*sub*vM.Calibration))^2) <= max_sp*vM.DotSpacing) && ...
%                 (1/sqrt(((mid_pt_sub(1)-iy)/(Mfft*sub*vM.Calibration))^2 + ((mid_pt_sub(2)-jx)/(Nfft*sub*vM.Calibration))^2) >= min_sp*vM.DotSpacing)
%             dots_fft_mask(iy,jx) = true;
%         end
%     end
% end
outer_x = mid_pt_sub(2) + Nfft*sub/(min_sp*spacing)*cos(0:2/Mfft:2*pi);
outer_y = mid_pt_sub(1) + Mfft*sub/(min_sp*spacing)*sin(0:2/Nfft:2*pi);
inner_x = mid_pt_sub(2) + Nfft*sub/(max_sp*spacing)*cos(0:2/Mfft:2*pi);
inner_y = mid_pt_sub(1) + Mfft*sub/(max_sp*spacing)*sin(0:2/Nfft:2*pi);
dots_fft_mask = poly2mask(outer_x,outer_y,Mfft*sub,Nfft*sub) - poly2mask(inner_x,inner_y,Mfft*sub,Nfft*sub);

dots_img_fft_filt_half = dots_img_fft_filt(floor(mid_pt_sub(1)+0.5*sub):end,:);
dots_img_fft_filt_half_masked = dots_img_fft_filt_half.*dots_fft_mask(floor(mid_pt_sub(1)+0.5*sub):end,:);

[ypeak,xpeak] = find(dots_img_fft_filt_half_masked == max(dots_img_fft_filt_half_masked(:)));
ypeak = ypeak + (mid_pt_sub(1) - 1 + 0.5*(sub-1));

xpeak = mean(xpeak);
ypeak = mean(ypeak);

% len = 1/sqrt(((mid_pt_sub(1)-ypeak)/(Mfft*sub*vM.Calibration))^2 + ((mid_pt_sub(2)-xpeak)/(Nfft*sub*vM.Calibration))^2);
rot_angle = atan((ypeak - mid_pt(1))/(xpeak - mid_pt(2)));