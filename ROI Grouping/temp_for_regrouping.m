ix_axons_to_rois_new = ix_axons_to_rois;

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
    
axon_ids=axon_ids_new;
Ain_axons = Ain_new;
clear Ain_new ix_axons_to_rois_new axon_ids_new

%%

i1= 13; i2= 15;
corrcoef(dFF([ix_axons_to_rois{i1},ix_axons_to_rois{i2}],:)')
figure, plot(dFF(ix_axons_to_rois{i1},:)), hold on, plot(1+dFF(ix_axons_to_rois{i2},:))
get_var_ratio(dFF(ix_axons_to_rois{i1},:),dFF(ix_axons_to_rois{i2},:),1)

%%

[Ain_axons,ix_axons_to_rois,axon_ids] = hand_merge(Ain,i1,i2,ix_axons_to_rois);