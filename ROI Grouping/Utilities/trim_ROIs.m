% This function trims ROIs in the following order:
%   1. Trim overlap between ROIs
%   2. Trim pixels that are directly adjacent to another ROI
%   3. Trim disconnected regions of ROIs
%   4. Remove ROIs with too few pixels 
%
% Input:
%    Ain              *Thresholded* spatial filters matrix (num pixels x num ROIs)
%    dims             Vector of pixel dims of patch: [d1, d2]
%
% Output:
%    Ain              New spatial filters matrix 
%    

function Ain = trim_ROIs(Ain,dims)

d1 = dims(1);
d2 = dims(2);

N_ROIs = size(Ain,2);

% Smallest size allowable for ROI
min_pix = 8;

%% Remove overlaps

% Set all direct overlaps to zero
ix = sum(Ain,2)>1;
Ain(ix,:) = 0;

% Remove pixels that are adjacent to another ROI
A_new = Ain;
for k = 1:size(Ain,2)
    A_otherROIs = sum(Ain,2) - Ain(:,k);
    
    pix = find(Ain(:,k));
    for kk = 1:length(pix)
        
        % Get each pixel and corresponding neighborhood
        [r, c]  = ind2sub([d1, d2], pix(kk));
        box = [r,c+1; r,c-1; r+1,c; r-1,c];
        box(box(:,1) == 0,:) = [];
        box(box(:,1) == d1+1,:) = [];
        box(box(:,2) == 0,:) = [];
        box(box(:,2) == d2+1,:) = [];
        box_pix = sub2ind([d1,d2],box(:,1),box(:,2));
        
        % Check if pixel abuts neighbouring ROI
        ix = find(A_otherROIs(box_pix));
        if numel(ix) > 0
            A_new(pix(kk),k) = 0;
        end
    end
end
Ain = A_new;

%% Trim disconnected ROIs
for k = 1:size(Ain,2)
    
    % Get number of disconnected parts
    CC = bwconncomp(reshape(Ain(:,k),d1,d2),4);
    
    % If more than 1 connected regions, only keep the largest
    if CC.NumObjects > 1
        [~,large_roi] = max(cellfun('length',CC.PixelIdxList));
        A_new = zeros(d1*d2,1);
        A_new(CC.PixelIdxList{large_roi}) = 1;
        Ain(:,k) = A_new;
    end
    
end

%% Remove ROIs of size < min_pixels
ix = (sum(Ain,1) < min_pix);
Ain(:,ix) = [];

%% Show warning if some ROIs were moved
if size(Ain,2) ~= N_ROIs
    disp('Warning: some ROIs have been removed.')
end

%% Plot results 

figure, imagesc(reshape(sum(Ain,2),d1,d2))
set(gca, 'xtick', []); set(gca, 'ytick', []); 
axis tight; axis equal;
caxis([0,3])

