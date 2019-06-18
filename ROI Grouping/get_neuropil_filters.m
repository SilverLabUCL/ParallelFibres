% This function makes fake ROIs to determine distribution of SNR in the
% neuropil
%
% Input:
%    Ain              Spatial filters
%    Cn               Correlation image (after smoothing)
%    dims             Vector of pixel dims of patch: [d1, d2]
%    manual           Set to 1 to choose centres manually, 0 otherwise
% 
% Output:
%    A_neuropil       Spatial mask for fake spatial filters

function A_neuropil = get_neuropil_filters(Ain,Cn,dims,manual)

if nargin < 4 || isempty(manual)
    manual = 0;
end

d1 = dims(1);
d2 = dims(2);


% Approximate radius of each fake ROI
r = sqrt(mean(sum(Ain,1))/pi);
[xx, yy] = meshgrid(1:d2, 1:d1);

% Plot correlation image of raw data
figure, plot_contours(Ain,Cn,.9);
colormap(gray);% caxis(median(Cn(:))+[-1,1]*std(Cn(:)))
set(gca, 'xtick', []); set(gca, 'ytick', []);
set(gca,'FontSize',18)
axis tight; axis equal;

if manual == 1
    % Choose neuropil centers and make fake ROIs
    title('Click centers to calculate neuropil.');
    xlabel('Press enter to stop.', 'color', 'r');
    ind_np_center = [];
    A_neuropil = zeros(d1*d2,1);
    for kk = 1:100
        [tmp_x, tmp_y] = ginput(1);
        if numel(tmp_x) == 0
            break
        end
        plot(tmp_x, tmp_y, '*r', 'linewidth', 2);
        %ind_np_center = [ind_np_center; sub2ind([d1,d2], tmp_y, tmp_x)];
        % Add 
        A_ = sqrt((xx-tmp_x).^2 + (yy-tmp_y).^2) < r;
        A_neuropil(:,kk) = A_(:);
    end
else
    A_neuropil = zeros(d1*d2,100);
    Ain_all = sum(Ain,2)>0;
    A_neuropil_all = sum(A_neuropil,2)>0;
    
    A_ok = ones(d1*d2,1) - Ain_all; 
    count = 1; % number of attempts
    ix = 1; % number successful neuropil fake ROIs 
    while ix <= 100 && count <= 5000
        % Randomly sample some centre
        ind_ctr = randsample(find(A_ok),1);
        [tmp_y, tmp_x] = ind2sub([d1,d2],ind_ctr);
        A_ = sqrt((xx-tmp_x).^2 + (yy-tmp_y).^2) < r;
        A_ = A_(:);
        
        % Overlap with other ROIs as fraction of its size
        overlap_rois = sum(A_+ Ain_all > 1)/sum(A_);
        overlap_npl_rois =  sum(A_+ A_neuropil_all > 1)/sum(A_);
        
        % If < 10% overlap with other ROI, keep it
        if overlap_rois <= .1 &&  overlap_npl_rois <= .1 
            A_neuropil(:,ix) = A_;
            A_neuropil_all = sum(A_neuropil,2)>0;
            ix = ix+1;
        end
        count = count+1;
    end
    A_neuropil = A_neuropil(:,1:ix-1);    
end

