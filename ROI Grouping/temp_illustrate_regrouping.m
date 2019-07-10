
load('/Users/caycogajic/Documents/ParallelFibres/Data/FL87_180501_11_03_09/processed/FL87_180501_11_03_09_GroupedData.mat')


Ain_axons = Ain_axons{1};
Cn = Cn{1};
dFF_rois = dFF_rois{1};
ix_axons_to_rois = ix_axons_to_rois{1};

plot_grouped_rois(Ain_axons,Cn,dFF_rois,ix_axons_to_rois,acquisition_rate)

figure(1), 
hold on, plot([1,1+5/Pixel_size],[-1,-1],'k','LineWidth',3)
set(gca,'XColor','w','YColor','w')
%%

axons_to_illustrate = [2,4,5,6,10];

ix_reorder = [];
for i = axons_to_illustrate
    ix_reorder = [ix_reorder, ix_axons_to_rois{i}];
end

figure, imagesc(corrcoef(dFF_rois(ix_reorder,:)'))
caxis([-1,1])
colormap(bluewhitered)
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('ROI')
ylabel('ROI')
title('Grouped')
set(gca,'FontSize',18)

figure, imagesc(corrcoef(dFF_rois(sort(ix_reorder),:)'))
caxis([-1,1])
colormap(bluewhitered)
set(gca,'XTick',[])
set(gca,'YTick',[])
xlabel('ROI')
ylabel('ROI')
title('Ungrouped')
set(gca,'FontSize',18)
%%
for i = axons_to_illustrate
    figure
    c = 0;
    for j = (ix_axons_to_rois{i})
        hold on, plot((1:size(dFF_rois,2))/acquisition_rate,dFF_rois(j,:)+c); c= c+1.5;
    end
end