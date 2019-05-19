% This function thresholds the spatial maps
%
% Input:
%    Ain              Spatial filters matrix (num pixels x num ROIs)
%    dims             Vector of pixel dims of patch: [d1, d2]
%    thr              Fraction for thresholding (optional)
% 
% Output:
%    
%    A_thresh         Thresholded spatial filters matrix

function A_thresh = threshold_filters(Ain,dims,thr)

if nargin < 3 || isempty(thr)
    thr = 0.9;
end

d1 = dims(1);
d2 = dims(2);

%% Threshold each spatial filter

A_thresh = zeros(size(Ain));
for i = 1:size(Ain,2)
    A_temp = full(reshape(Ain(:,i),d1,d2));
    ix_nonzero = find(A_temp>0);
    A_temp(ix_nonzero) = A_temp(ix_nonzero) -min(A_temp(ix_nonzero));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    A_temp = max(A_temp,0);
    [temp,ind] = sort(A_temp(:).^2,'ascend');
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    ix = A_temp>A_temp(ind(ff));
    A_thresh(ix,i) = 1;
end

%% Plot thresholded 

figure, imagesc(reshape(sum(A_thresh,2),d1,d2))
set(gca, 'xtick', []); set(gca, 'ytick', []); 
axis tight; axis equal;
caxis([0,3])
