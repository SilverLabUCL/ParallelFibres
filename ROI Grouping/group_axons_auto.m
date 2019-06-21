% This script groups ROIs onto putative axons, separately for each patch in
% each experiment

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
            'FL45_170125_14_47_04'}; %  8

% Choose dataset
dataset_ix = 1;
fname = datasets{dataset_ix};
disp(fname)

fname_fibreangles = [basedir,fname,'/processed/fibre_direction.mat'];
load(fname_fibreangles,'vector_mean');

fname_SNR = [basedir,fname,'/processed/SNR.mat'];
load(fname_SNR,'SNR_thresh');

load([basedir,fname,'/',fname,'.mat'],'Numb_patches')

%% Choose patch number
for patch_no = 1:Numb_patches
    %%
    load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
    Y = double(Y);
    Cn = correlation_image(Y, [1,2], d1,d2);

    % Use CNMF initialization to estimate initial spatial filters
    [Ain,~] = detect_ROIs(Y, [d1,d2]);

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

%%    % Group axons
    [Ain_axons,ix_axons_to_rois,axon_ids] = ...
        get_axon_grouping(Ain,Y,[d1,d2],acquisition_rate,Pixel_size,vector_mean);
%%
    % Plot results
    plot_grouped_rois(Ain_axons,Cn,dFF,ix_axons_to_rois,acquisition_rate,[],[],[basedir,fname],patch_no);%,thr,display_numbers,basedir)

    pause
    
    % Save figures in png
        figure(1), saveas(gcf,[basedir,fname,'/figs/spatial_grouped_Patch',sprintf('%03d',patch_no),'.png']);
        figure(2), saveas(gcf,[basedir,fname,'/figs/dFF_grouped_Patch',sprintf('%03d',patch_no),'.png']);
        figure(3), saveas(gcf,[basedir,fname,'/figs/spatial_loners_Patch',sprintf('%03d',patch_no),'.png']);
        figure(4), saveas(gcf,[basedir,fname,'/figs/dFF_loners_Patch',sprintf('%03d',patch_no),'.png']);
    
    % Save results
    Ain_rois = Ain; dFF_rois = dFF;
    dFF_axons = get_dFF(Ain_axons,Y,acquisition_rate);
    save([basedir,fname,'/processed/GroupedData_Patch',sprintf('%03d',patch_no),'.mat'],'Cn','dFF_axons','Ain_axons','dFF_rois','Ain_rois','ix_axons_to_rois','axon_ids');
    
end



