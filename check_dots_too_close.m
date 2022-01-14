fig_1 = figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
ax_1 = axes;
im_1 = imagesc(imadjust(vidFrames_crop_raw(:,:,1)));
axis image
hold on
p_dot = plot(px(:,vMcrop.uFrame),py(:,vMcrop.uFrame),'.r','markersize',10);

k = vM.uFrame;
for idot = 1:length(px)
    dist = sqrt(((px(:,k) - px(idot,k)).^2 + (py(:,k) - py(idot,k)).^2));
    if nnz(dist < 0.5*dotspacing) > 1
        [idot,px(idot,k),py(idot,k)]
        plot(px(idot,k),py(idot,k),'ok','markersize',10)
        1;
    end
end