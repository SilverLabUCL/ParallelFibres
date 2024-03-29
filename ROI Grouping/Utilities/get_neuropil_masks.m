% This function makes neuropil masks to compare signal in rois vs neuropil
% fluorescence, i.e., for Supp figure 4
% (NOT the filters used for grouping criterion - for that see
% get_neuropil_filters.m)
%
% Input:
%    Ain              Spatial filters
%    dims             Vector of pixel dims of patch: [d1, d2]
%    Pixel_size       Size of each pixel in um
%    npl_msk_size     1/2 length of square mask (optional)
% 
% Output:
%    Ain_npl       Spatial mask for fake spatial filters

function Ain_npl = get_neuropil_masks(Ain,Cn,dims,Pixel_size,npl_msk_size)


if nargin < 5 || isempty(npl_msk_size)
    npl_msk_size = round(10 / Pixel_size);
end


d1 = dims(1);
d2 = dims(2);

N_ROIs = size(Ain,2);

% Average radius of detected varicosity
r = sqrt(mean(sum(Ain,1))/pi);

% Threshold on correlation image
rho_thresh = prctile(Cn(:),95);

% First generate one mask that defines where each neuropil mask will avoid
% i.e., any pixel within 2 * r of a varicosity centre
Ain_mask = zeros(d1,d2);
centers = zeros(N_ROIs,2);

for roi = 1:N_ROIs
    % Get center of that varicosity
    Ain_temp = reshape(Ain(:,roi),d1,d2);
    c = regionprops(Ain_temp,'centroid'); 
    c = round(c.Centroid);
    centers(roi,:) = c;

    % Any pixel within 2*r of the varicosity center
    for i = max(1,c(2)-floor(2*r)) : min(d1,c(2)+ceil(2*r))
        for j = max(1,c(1)-floor(2*r)) : min(d2,c(1)+ceil(2*r))
            dist_temp = sqrt((c(2)-i)^2 + (c(1)-j)^2);
            if dist_temp <= 2*r
                Ain_mask(i,j) = 1; 
            end
        end
    end
end

for i = 1:d1
    for j = 1:d2
        if Cn(i,j) >= rho_thresh
            Ain_mask(i,j) = 1; 
        end
    end
end

% Next generate the actual neuropixel masks for each ROI
Ain_npl = zeros(size(Ain));

for roi = 1:N_ROIs
    % Reload varicosity center
    c = centers(roi,:);

    % Generate mask - square of size npl_msk_size, removing any pixel
    % inside the mask 
    Ain_temp = zeros(d1,d2);
    for i = max(1,c(2)-npl_msk_size) : min(d1,c(2)+npl_msk_size)
        for j = max(1,c(1)-npl_msk_size) : min(d2,c(1)+npl_msk_size)
            if Ain_mask(i,j)==0
                Ain_temp(i,j) = 1;
            end
        end
    end

    Ain_npl(:,roi) = Ain_temp(:);
end

if any(sum(Ain_npl,1) < 50)
    error('Neuropil mask is very small')
end
