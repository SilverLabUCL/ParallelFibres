
% Cn - each cell contains the 'correlation image' for that patch, you can
% use this to visualize the ROIs, e.g.

figure, imagesc(Cn{1}), colormap(gray)

[d1,d2] = size(Cn{1}); % d1 and d2 are the size of the patch
T = size(dFF_axons{1},2); % Number of timepoints

%% The spatial and functional info for the grouped ROIs are contained in:
%     dFF_axons - each cell contains dFF after grouping within that patch
%     Ain_axons - each cell contains a matrix whose kth element is the mask
%              of that axon, but it has to be reshaped from a vector to a
%              matrix

%% For example, this plots the dFF of the 2nd axon in patch 1

figure, plot((1:T)/acquisition_rate,dFF_axons{1}(2,:))

% and this plots the corresponding spatial mask - note that we have to
% reshape the vector of length d1*d2, into a matrix of length d1 & width d2

figure, imagesc(reshape(Ain_axons{1}(:,2),d1,d2))
axis tight equal
%% The following variables have the same spatial & functional information for 
% ungrouped ROIs
%     dFF_rois
%     Ain_rois

% You can switch between the axon IDs and ROI IDs using these variables
%    ix_axons_to_rois - each cell contains a cell array where each element
%                    is the ID of all ROIs corresponding to that axon
%    axon_ids  - vector where each element is the axon ID for that
%                   particular ROI

% So for example to find which ROIs are associated with the 2nd axon in
% patch 1, 

roi_ids_for_axon2 = ix_axons_to_rois{1}{2};

% And now to plot all of the dFF for those ROIs in black, and the dFF for
% the entire axon in blue

figure, hold on
plot((1:T)/acquisition_rate,dFF_axons{1}(2,:),'b')

for k = 1:length(roi_ids_for_axon2)
    roi = roi_ids_for_axon2(k);
    plot((1:T)/acquisition_rate,dFF_rois{1}(roi,:) + k,'k')
end

%% That was just to get you acquainted with the variables. I've written the 
% function plot_grouped_rois to make this a lot easier to do automatically

% So again for patch 1
plot((1:T)/acquisition_rate,dFF_axons{1}(2,:),'b','LineWidth',1.5)

% beware that calling this function will close all figures that you have
% open!!!



