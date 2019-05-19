% This function calculates the inter-bouton distances, useful to compare to
% literature to see if algorithm is in right ballpark
%
% Input:
%    Ain                Spatial filters matrix (num pixels x num ROIs)
%    dims               Vector of pixel dims of patch: [d1, d2]
%    ix_axons_to_rois   Cell mapping from axon id to rois
%    Pixel_size         Microns per pixel
% 
% Output:
%    dists              Vector of inter bouton distances (um)

function dists = get_interbouton_dist(Ain,dims,ix_axons_to_rois,Pixel_size)
    
    d1 = dims(1);
    d2 = dims(2);

    N_ROIs = size(Ain,2);
    N_axons = size(ix_axons_to_rois);
        
    dists = [];
    
    a = 1;
    N_ROIs_this_axon = numel(ix_axons_to_rois{a});
    
    while N_ROIs_this_axon > 1
        
        % Calculate centres for all rois in this axon
        centers = zeros(N_ROIs_this_axon,2);
        for roi = 1:N_ROIs_this_axon
            c = regionprops(reshape(Ain(:,ix_axons_to_rois{a}(roi)),d1,d2),'centroid'); 
            centers(roi,:) = c.Centroid;
        end
        
        % Sort rois by d1 (short side of patch) - fine as long as fibres
        % are always oriented in similar way
        [~,ix] = sort(centers(:,2));
        centers = centers(ix,:);
        
        % Now calculate distance between adjacent rois
        for roi = 1:(N_ROIs_this_axon-1)
            dist_temp = sqrt(sum((centers(roi+1,:) - centers(roi,:)).^2));
            dists = [dists; dist_temp];
        end
 
        % Go to next axon
        a = a+1;
        N_ROIs_this_axon = numel(ix_axons_to_rois{a});
    end
    
    dists = dists * Pixel_size;
    
    