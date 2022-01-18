ic = 1;

k = vM.uFrame;

L = 0;

[n_row,n_col] = size(celldata(ic).Xgrid);

Xvector = celldata(ic).Xvector;
Yvector = celldata(ic).Yvector;

Xdisp_k_dn = celldata(ic).Xdisp_k;
Ydisp_k_dn = celldata(ic).Ydisp_k;
% Xdisp_k_dn = celldata(ic).Xdisp_k_dn;
% Ydisp_k_dn = celldata(ic).Ydisp_k_dn;

vMcrop = celldata(ic).vM;

real_points = celldata(ic).real_points;
celldots = celldata(ic).celldots;

E = vMcrop.YoungsModulus;
nu = vMcrop.Poisson;
% d = mean(BD.DotSpacings);
d = 5.89;

u_x = reshape(Xdisp_k_dn(:,k),n_row,n_col)*vMcrop.Calibration;
u_y = reshape(Ydisp_k_dn(:,k),n_row,n_col)*vMcrop.Calibration;
u = cat(3,u_x,u_y);

pos_x = reshape(Xvector,n_row,n_col)*vMcrop.Calibration;
pos_y = reshape(Yvector,n_row,n_col)*vMcrop.Calibration;
pos = cat(3,pos_x,pos_y);

test_unconstrained = calcforce_regFTTC(u,E,nu,d,L,real_points);
test = calcforce_regFTTC_constrained(u,E,nu,d,L,real_points,celldots,pos);
test_mag = vecnorm(test,2,3);
test_unconstrained_mag = vecnorm(test_unconstrained,2,3);

%%
figure
subplot(1,2,1)
imagesc(test_unconstrained_mag)
hold on
plot(y,x,'sr')
hold off
axis image
title('Unconstrained FTTC')
subplot(1,2,2)
imagesc(test_mag)
hold on
[x,y] = find(celldots);
plot(y,x,'sr')
hold off
axis image
title('Constrained FTTC')