% This function detects ROIs using code adapted from CNMFE initialisation
%
% Input:
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    dims             Vector of pixel dims of patch: [d1, d2]
%    filter_pix       Spatial scale (in pixels)
%    win_pix          Size of square window (in pixels) to draw box around
%                     each ROI (approx size of bouton or soma)
% 
% Output:
%    Ain              *Thresholded* spatial filters matrix (num pixels x num ROIs)
%    Cn               Correlation matrix


function [Ain,Cn] = detect_ROIs(Y, dims, filter_pix, win_pix, thresh, manual)


if nargin < 6 || isempty(manual)
    manual = 0;
end

if nargin < 5 || isempty(thresh)
    thresh = 1;
end

min_Cn = .83;

%% Load size parameters

T = size(Y,2);

d1 = dims(1);
d2 = dims(2);
if d1*d2 ~= size(Y,1)
    error('Dimensions are incorrect');
end

% Make sure win_pix is an integer
win_pix = round(win_pix);

%% Spatially filtered data
[HY,Cn] = get_spat_filtered_data(Y,[d1,d2],filter_pix,win_pix);

%% Look for local maxima

% screen seeding pixels as center of the neuron
v_search = Cn;

% Matrix showing whether each pixel has been searched before
ind_search = false(d1*d2,1);  

% Spatial filter and an extra value to avoid repeated seed pixels within one ROI
[ii, jj] = meshgrid(1:d2, 1:d1);
pixel_v = (ii*10+jj)*(1e-15);
v_search = v_search+pixel_v; 

% To display when manually looking for seed points
v_search_display = v_search;

% ignore pixels with small correlations or low peak-noise-ratio
min_pixel = 8;
min_v_search = median(v_search(:));
ind_search(v_search<min_v_search) = true; 
v_search(ind_search) = 0;

%%
tmp_d = max(3, round(win_pix/4));
v_max = ordfilt2(v_search, tmp_d^2, true(tmp_d));
v_max(v_max < min_Cn) = 0;
ind_localmax = find(and(v_search(:)==v_max(:), v_max(:)>0));
[r_peak, c_peak] = ind2sub([d1,d2],ind_localmax);

% Remove points on the boundary
ix = find(r_peak ~= 1 & r_peak ~= d1 & c_peak ~= 1 & c_peak ~= d2);
ind_localmax = ind_localmax(ix);
r_peak = r_peak(ix);
c_peak = c_peak(ix);

% Plot and hand select more
if manual
    figure('papersize',[d2,d1]/40); init_fig;
    imagesc(v_search_display), colormap(gray)
    hold on, plot(c_peak,r_peak,'.r')
    caxis(median(v_search_display(:))+[-1,2]*std(v_search_display(:)))

    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    title('Click additional ROI centers.');
    xlabel('Press enter to stop.', 'color', 'r');
    for kk = 1:100
        [tmp_x, tmp_y] = ginput(1);
        if numel(tmp_x) == 0
            break
        end
        plot(tmp_x, tmp_y, '*r', 'linewidth', 2);
        ind_localmax = [ind_localmax; sub2ind([d1,d2], round(tmp_y), round(tmp_x))];
    end
end
% Sort according to v_search
[~, ind_sort] = sort(v_search_display(ind_localmax), 'descend');
ind_localmax = ind_localmax(ind_sort);

%% Initialize components

searching_flag = true; 
k = 0;      %number of found components

K = floor(sum(v_search(:)>min_v_search)/10); % estimate max number of components
Ain = zeros(d1*d2, K);  % spatial components
Cin = zeros(K, T);      % temporal components
Sin = zeros(K, T);    % spike counts
Cin_raw = zeros(K, T);
kernel_pars = cell(K,1);    % parameters for the convolution kernels of all neurons
center = zeros(K, 2);   % center of the initialized components

%% Try initialization over all local maximums

for mcell = 1:length(ind_localmax)
    % find the starting point
    ind_p = ind_localmax(mcell);
    max_v = v_search(ind_p);
    if mcell==1
        img_clim = [0, max_v];
    end
    ind_search(ind_p) = true; % indicating that this pixel has been searched.

    [r, c]  = ind2sub([d1, d2], ind_p);

    % roughly check whether this is a good starting point
    y0 = HY(ind_p, :);
    y0_std = std(diff(y0));
    if max(diff(y0))< 3*y0_std % signal is weak
        continue;
    end

    % select its neighbours for estimation of ai and ci, the box size is
    rsub = max(1, -win_pix+r):min(d1, win_pix+r);
    csub = max(1, -win_pix+c):min(d2, win_pix+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    HY_box = HY(ind_nhood, :);      % extract temporal component from HY_box
    Y_box = Y(ind_nhood, :);    
    %ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % subscripts of the center

    % Linear regression to guess spatial filter
    ai = (y0*y0')\(y0*HY_box');
    ai = ai';
    
    if any(isnan(ai)) || sum(ai)<=1 || sum(ai(:)>0)<min_pixel
        continue;
    else
        k = k+1;

        Ain(ind_nhood, k) = ai;
        center(k, :) = [r, c];
       
        % avoid searching nearby pixels
        ind_search(ind_nhood(ai>max(ai)*0.5)) = true;

        % Remove component from filtered data
        HY_box = HY(ind_nhood, :) - ai*y0;
        HY(ind_nhood, :) = HY_box;

    end

end

Ain = Ain(:, 1:k); 
center = center(1:k,:);
%% Threshold filters
if thresh
    Ain = threshold_filters(Ain,[d1,d2]);
end
%%

figure 
plot_contours(Ain,Cn,.8,true,[],[],2);
colormap('gray')
set(gca, 'xtick', []); set(gca, 'ytick', []); 
axis tight; axis equal;