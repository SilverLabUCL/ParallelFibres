
function [Ain_new,ix_axons_to_rois_new,axon_ids_new] = hand_merge(Ain,ax1,ax2,ix_axons_to_rois)

    % Merge axons and update indices
    ix_axons_to_rois_new = ix_axons_to_rois;
    ix_axons_to_rois_new{ax1} = [ix_axons_to_rois_new{ax1},ix_axons_to_rois_new{ax2}];
    ix_axons_to_rois_new(ax2) = [];

    % New number of axons
    num_axons = length(ix_axons_to_rois_new);

    % Get axon IDs
    axon_ids_new = zeros(1,size(Ain,2));
    for axon = 1:num_axons
        rois = ix_axons_to_rois_new{axon};
        axon_ids_new(rois) = axon;
    end

    Ain_new = zeros(size(Ain,1),num_axons);
    for axon = 1:num_axons
        for roi = ix_axons_to_rois_new{axon}    
            Ain_new(:,axon) = Ain_new(:,axon) + Ain(:,roi);
        end
    end