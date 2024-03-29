% This script plots spatial information that wasn't used for regrouping, as
% a sanity check to compare with existing literature on parallel fibre
% bouton information (Pichitpornchai et al., 1994; Shepherd et al., 2002)

clear all; clc
addpath('From CNMF_E/')
addpath('Utilities/')

basedir = '~/Documents/ParallelFibres/Data/';

datasets = {'FL87_180501_11_03_09',...  1
            'FL87_180501_10_47_25',...  2
            'FL87_180501_10_36_14',...  3 
            'FL87_180220_10_38_55',...  4
            'FL77_180213_10_46_41',...  5
            'FL_S_170906_11_26_25',...  6
            'FL_S_170905_10_40_52',...  7            
             'FL95_180425_10_53_40',... 8
             'FL87_180413_11_00_55',... 9
             'FL87_180117_11_23_20',... 10
             'FL_S_171109_14_54_34',... 11
             'FL_S_171109_15_19_52',... 12
             'FL76_170913_10_57_06'}; % 13

%% Calculate intervaricosity distances for all experiments, every patch

d_all = cell(13,1);

figure, hold on, plot(1:9,'k')
for dataset_ix = 1:13
    
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

[mean(vertcat(d_all{:})), std(vertcat(d_all{:}))]

%% Info from Pichitpornchai et al Table 2

mean_per_zone = [3.71 4.35 4.74 5.46 5.72 5.46 6.45 7.44];
estimates_per_zone =  [2315 1974 1812 1573 1501 1570 1331 1154];
ste_per_zone = [0.06 0.07 0.07 0.09 0.12 0.11 0.14 0.23];

num_estimates_all = sum(estimates_per_zone);
mean_all = sum(mean_per_zone.*estimates_per_zone)/num_estimates_all;

std_per_zone = sqrt(estimates_per_zone).*ste_per_zone;
var_per_zone = std_per_zone.^2;

var_all = sum((var_per_zone + mean_per_zone.^2).*estimates_per_zone) / num_estimates_all - mean_all^2
sqrt(var_all)/sqrt(num_estimates_all)

%% Histogram of intervaricosity distances - compare to Pichitpornchai et al
% Figure 1d

figure, histogram(vertcat(d_all{:}),'FaceColor','r')
hold on, plot(mean(vertcat(d_all{:})),260,'vr','MarkerFaceColor','r')
hold on, plot(5.2,300,'vk','MarkerFaceColor','k')
hold on, plot([1,14],[310,310],'k','LineWidth',1.5)
hold on, plot([1,1],[305,315],'k','LineWidth',1.5)
hold on, plot([14,14],[305,315],'k','LineWidth',1.5)
set(gca,'box','off')
set(gca,'FontSize',15)
xlabel('Distance (\mu m)')
ylabel('Number')

%% Cumulative distribution - compare to Shepherd et al 

x = 0:.01:12;
d_temp = vertcat(d_all{:});
for k = 1:length(x)
    y(k) = sum(d_temp>x(k))/numel(d_temp);
end

figure, plot(x,y,'k','LineWidth',2);
set(gca,'FontSize',18)
xlabel('Distance (\mu m)')
ylabel('Probability')

%% Calculate number of ROIs per axon

N_ROIs = 0;
N_axons = 0;
for dataset_ix = 1:15
    fname = datasets{dataset_ix};
    
    load([basedir,fname,'/processed/',fname,'_GroupedData.mat'],'Ain_rois','Ain_axons')
    N_ROIs = N_ROIs + size(horzcat(Ain_rois{:}),2);
    N_axons = N_axons + size(horzcat(Ain_axons{:}),2);
end

N_ROIs / N_axons

