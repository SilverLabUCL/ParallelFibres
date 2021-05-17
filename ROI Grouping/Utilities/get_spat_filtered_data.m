% This function returns data and correlation image after light spatial filtering
% 
% Input:
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    dims             Vector of pixel dims of patch: [d1, d2]
%    filter_pix       Spatial scale (in pixels)
%    win_pix          Size of square window (in pixels) to draw box around
%                     each ROI (approx size of bouton or soma)
% 
% Output:
%    HY               Spatially filtered data
%    Cn               Correlation matrix


function [HY,Cn] = get_spat_filtered_data(Y,dims,filter_pix,win_pix)

    d1 = dims(1);
    d2 = dims(2);

    % Gaussian filter within local window
    psf = fspecial('gaussian',win_pix, filter_pix);

    % Tensor (d1 x d2 x T) of each frame spacially smoothed by Gaussian filter
    HY = imfilter(reshape(Y, d1,d2,[]), psf, 'replicate');
    HY = reshape(HY, d1*d2, []);

    % Remove median from each pixel, removes background
    HY = bsxfun(@minus, HY, median(HY, 2));

    % Correlation image
    Cn = correlation_image(HY, [1,2], d1,d2);
