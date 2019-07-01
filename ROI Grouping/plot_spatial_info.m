% This script plots spatial information that wasn't used for regrouping, as
% a sanity check to compare with existing literature on parallel fibre
% bouton information (Pichitpornchai et al., 1994; Shepherd et al., 2002)

clear all; clc
addpath('From CNMF_E/')
addpath('Utilities/')

%basedir = '~/Documents/ParallelFibres/Data/';
basedir = 'C:\Users\SilverLab\Documents\Alex\ParallelFibres\Data\';

datasets = {'FL87_180501_11_03_09',...  1
            'FL87_180501_10_47_25',...  2
            'FL87_180501_10_36_14',...  3 
            'FL87_180220_10_38_55',...  4
            'FL77_180213_10_46_41',...  5
            'FL_S_170906_11_26_25',...  6
            'FL_S_170905_10_40_52',...  7            
            'FL45_170125_14_47_04'}; %  8

d_all = cell(8,1);


figure, hold on, plot(1:9,'k')
for dataset_ix = 1:8
    
    fname = datasets{dataset_ix};
    
    load([basedir,fname,'/processed/',fname,'_GroupedData.mat'],'Ain_rois','Cn','ix_axons_to_rois');
    
    load([basedir,fname,'/raw/Patch001.mat'],'d1','d2','Pixel_size')
    
    Numb_patches = size(Cn,1);

    d = [];
    for patch_no = 1:Numb_patches
        temp = get_interbouton_dist(Ain_rois{patch_no},[d1,d2],ix_axons_to_rois{patch_no},Pixel_size);
        d = [d; temp];
    end
    
    d_all{dataset_ix} = d;
    
    [dataset_ix, mean(d), std(d)]
    plot(mean(d),std(d),'o','MarkerEdgeColor','k','MarkerFaceColor','w')
end
set(gca,'FontSize',15)
xlabel('Mean distance (\mu m)')
ylabel('Standard deviation (\mu m)')
saveas(gcf,[basedir,'/figs/intervaricosity_distances_mean_std.png']);

figure, histogram(vertcat(d_all{:}))
set(gca,'FontSize',15)
xlabel('Intervaricosity distance (\mu m)')
ylabel('Number')
saveas(gcf,[basedir,'/figs/intervaricosity_distances_histogram.png']);

[mean(vertcat(d_all{:})), std(vertcat(d_all{:}))]
