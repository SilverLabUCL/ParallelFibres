% This script groups ROIs onto putative axons, separately for each patch in
% each experiment

clear all; clc; close all
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
            'FL45_170125_14_47_04'}; %  8

% Choose dataset
dataset_ix = 6;
fname = datasets{dataset_ix};
disp(fname)

fname_fibreangles = [basedir,fname,'/processed/fibre_direction.mat'];
load(fname_fibreangles,'vector_mean');

fname_SNR = [basedir,fname,'/processed/SNR.mat'];
load(fname_SNR,'SNR_thresh');

fname_corr = [basedir,fname,'/processed/corr_histograms.mat'];
load(fname_corr,'C_inter');
rho_min = prctile(C_inter,95);
    

load([basedir,fname,'/',fname,'.mat'],'Numb_patches')


Cn_all = cell(Numb_patches,1);
dFF_axons_all = cell(Numb_patches,1);
Ain_axons_all = cell(Numb_patches,1);
dFF_rois_all = cell(Numb_patches,1);
Ain_rois_all = cell(Numb_patches,1);
ix_axons_to_rois_all = cell(Numb_patches,1);
axon_ids_all = cell(Numb_patches,1);

% Choose patch number
for patch_no = 1:Numb_patches
    %
    disp([num2str(patch_no),' / ',num2str(Numb_patches)])
    load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
    Y = double(Y);
    Cn = correlation_image(Y, [1,2], d1,d2);

    % Use CNMF initialization to estimate initial spatial filters
    [Ain,~] = detect_ROIs(Y, [d1,d2], [],[],[],1);
    close all

    % Remove low SNR ROIs
    [Ain,~,~,~] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,SNR_thresh); 

    % Trim ROIs to remove 
    Ain = trim_ROIs(Ain,[d1,d2]);

    % Remove low SNR ROIs again
    [Ain,SNR,A_bad,SNR_bad] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,SNR_thresh); 

    % Calculate dFF
    dFF = get_dFF(Ain,Y,acquisition_rate);

    % Plot all dFF
    figure, hold on
    for k = 1:size(dFF,1)
        plot((1:size(dFF,2))/acquisition_rate,dFF(k,:)+k)
    end
    
    % Manually delete ROIs due to slow drift
    [Ain,dFF] = delete_ROIs_manually(Ain,dFF,acquisition_rate);

    % Group axons
    [Ain_axons,ix_axons_to_rois,axon_ids] = ...
        get_axon_grouping(Ain,Y,[d1,d2],acquisition_rate,Pixel_size,vector_mean,[],rho_min,[],[],1);

    % Plot results
    plot_grouped_rois(Ain_axons,Cn,dFF,ix_axons_to_rois,acquisition_rate,[],[],[basedir,fname],patch_no);%,thr,display_numbers,basedir)

    disp('Press enter to accept')
    pause
    
    % Save figures in png
        figure(1), saveas(gcf,[basedir,fname,'/figs/spatial_grouped_Patch',sprintf('%03d',patch_no),'.png']);
        figure(2), saveas(gcf,[basedir,fname,'/figs/dFF_grouped_Patch',sprintf('%03d',patch_no),'.png']);
        figure(3), saveas(gcf,[basedir,fname,'/figs/spatial_loners_Patch',sprintf('%03d',patch_no),'.png']);
        figure(4), saveas(gcf,[basedir,fname,'/figs/dFF_loners_Patch',sprintf('%03d',patch_no),'.png']);
    close all
        
    % Save results
    Ain_rois = Ain; dFF_rois = dFF;
    dFF_axons = get_dFF(Ain_axons,Y,acquisition_rate);
    save([basedir,fname,'/processed/GroupedData_Patch',sprintf('%03d',patch_no),'.mat'],'Cn','dFF_axons','Ain_axons','dFF_rois','Ain_rois','ix_axons_to_rois','axon_ids','acquisition_rate');
    
    % To save all data later
    Cn_all{patch_no} = Cn;
    dFF_axons_all{patch_no} = dFF_axons;
    Ain_axons_all{patch_no} = Ain_axons;
    dFF_rois_all{patch_no} = dFF_rois;
    Ain_rois_all{patch_no} = Ain_rois;
    ix_axons_to_rois_all{patch_no} = ix_axons_to_rois;
    axon_ids_all{patch_no} = axon_ids;

end


%% Save all data

Cn = Cn_all;
dFF_axons = dFF_axons_all;
Ain_axons = Ain_axons_all;
dFF_rois = dFF_rois_all;
Ain_rois = Ain_rois_all;
ix_axons_to_rois = ix_axons_to_rois_all;
axon_ids = axon_ids_all;

save([basedir,fname,'/processed/GroupedData_AllPatches.mat'],...
    'Cn','dFF_axons','Ain_axons','dFF_rois','Ain_rois','ix_axons_to_rois','axon_ids','acquisition_rate');

disp('All data successfully saved.')