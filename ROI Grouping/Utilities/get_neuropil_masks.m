% This function makes filters for for neuropil subtraction
%
% Input:
%    Ain              Spatial filters
%    dims             Vector of pixel dims of patch: [d1, d2]
%    Pixel_size       Size of each pixel in um
%    npl_msk_size     Size of square mask (optional)
% 
% Output:
%    Ain_npl       Spatial mask for fake spatial filters

function Ain_npl = get_neuropil_masks(Ain,dims,Pixel_size,npl_msk_size)


if nargin < 4 || isempty(npl_msk_size)
    npl_msk_size = round(10 / Pixel_size);
end


d1 = dims(1);
d2 = dims(2);

N_ROIs = size(Ain,2);

% Average radius of detected varicosity
r = sqrt(mean(sum(Ain,1))/pi);

% First generate one mask that defines where each neuropil mask will avoid
% i.e., any pixel within 3 * r of a varicosity centre
Ain_mask = zeros(d1,d2);
centers = zeros(N_ROIs,2);

for roi = 1:N_ROIs
    % Get center of that varicosity
    Ain_temp = reshape(Ain(:,roi),d1,d2);
    c = regionprops(Ain_temp,'centroid'); 
    c = round(c.Centroid);
    centers(roi,:) = c;

    % Any pixel within 3*r of the varicosity center
    for i = max(1,c(2)-floor(3*r)) : min(d1,c(2)+ceil(3*r))
        for j = max(1,c(1)-floor(3*r)) : min(d2,c(1)+ceil(3*r))
            dist_temp = sqrt((c(2)-i)^2 + (c(1)-j)^2);
            if dist_temp <= 3*r
                Ain_mask(i,j) = 1; 
            end
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
    
