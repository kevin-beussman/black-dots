function [L_out,i_out,i_all] = find_L(resid_norm,soln_norm,L,init_L,method)
    if nargin < 5
        method = 'corner';
    end

    x = log(resid_norm);
    y = log(soln_norm);
    
    L_slope = L(1:end-4);
    x_slope = x(1:end-4);
%     slope = diff(y)./diff(x);
    for k = 1:length(x_slope)
        pca_coeff = pca([x(k:k+4)', y(k:k+4)']);
        slope(k) = pca_coeff(2,1)/pca_coeff(1,1);
    end

    L_curv = L_slope(1:end-2);
    x_curv = x_slope(1:end-2);
%     curv = diff(slope)./diff(x_slope);
    for k = 1:length(x_curv)
        pca_coeff = pca([x_slope(k:k+2)', slope(k:k+2)']);
        curv(k) = pca_coeff(2,1)/pca_coeff(1,1);
    end
    
    if strcmp(method,'corner')
        nSections = 3;
        nPoints = length(curv);
        p=0;
        i_maxcurv = [];
        for ii=1:nSections
            [~, i_mc] = max(curv(floor((ii-1)*nPoints/nSections)+1:floor(ii*nPoints/nSections)));
            % this is right at the L-corner which is usually over-smoothing
            i_mc = i_mc + floor((ii-1)*nPoints/nSections);
            % check if this is truly local maximum
            if i_mc>1 && i_mc<nPoints && (curv(i_mc)>curv(i_mc-1) && curv(i_mc)>curv(i_mc+1))
                p=p+1;
                i_maxcurv(p) = i_mc;
            end
        end

        if length(i_maxcurv)==1
            i_maxcurv = i_maxcurv(1);
        elseif length(i_maxcurv)>1
            % pick the one which is closer to initial lambda
            [~,Idx_close] = min(abs(log(L_curv(i_maxcurv))-log(init_L)));
            i_maxcurv = i_maxcurv(Idx_close);
        elseif isempty(i_maxcurv)
    %         disp('There is no local maximum in curvature in the input lambda range.Using global maximum instead ...');
            [~, i_maxcurv] = max(curv);
        end
        
        i_out = i_maxcurv + 1;
        if i_out > length(L_curv)
            i_out = length(L_curv);
        elseif i_out < 1
            i_out = 1;
        end
        L_out = L_curv(i_out);
        
    elseif strcmp(method,'optimal')
        nPoints = length(curv);
        
        % find all zero-crossing points
        i_zerocurv = find(diff(sign(curv)));
        i_all = i_zerocurv + 1;

        if length(i_zerocurv)==1
            i_zerocurv = i_zerocurv(1);
        elseif length(i_zerocurv)>1
            % pick the one which is closer to initial lambda
            [~,Idx_close] = min(abs(log(L_curv(i_zerocurv))-log(init_L)));
            i_zerocurv = i_zerocurv(Idx_close);
        elseif isempty(i_zerocurv)
            [~, i_zerocurv] = min(abs(curv));
        end
        
        i_out = i_zerocurv + 1;
        if i_out > length(L_curv)
            i_out = length(L_curv);
        elseif i_out < 1
            i_out = 1;
        end
        L_out = L_curv(i_out);
        
    end
    
end