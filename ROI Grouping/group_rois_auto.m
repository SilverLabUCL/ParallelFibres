% This script groups ROIs onto putative axons, separately for each patch in
% each experiment

clear all; clc; close all
addpath('From CNMF_E/')
addpath('Utilities/')

% basedir = '~/Documents/ParallelFibres/Data/';
% datasets = {'FL87_180501_11_03_09',...  1
%             'FL87_180501_10_47_25',...  2
%             'FL87_180501_10_36_14',...  3 
%             'FL87_180220_10_38_55',...  4
%             'FL77_180213_10_46_41',...  5
%             'FL_S_170906_11_26_25',...  6
%             'FL_S_170905_10_40_52',...  7            
%             'FL45_170w125_14_47_04',...  8
%              ...%%
%              'FL95_180425_10_53_40',...  9
%              'FL87_180413_11_00_55',...  10
%              'FL104_180725_10_42_37',... 11
%              'FL87_180117_11_23_20',...  12
%              'FL106_180807_10_52_25',... 13
%              'FL_S_171109_14_54_34',...  14
%              'FL_S_171109_15_19_52',...  15
%              ...%
%              'FL75_170912_10_33_29',... 16
%              'FL76_170913_10_57_06',... 17
%              'FL77_180113_10_58_50'};  %18

% basedir = '~/Documents/ParallelFibres/Data_Puff/';         
% datasets = {'FL41_161123_15_27_16',... 1
%             'FL46_161109_11_43_30',... 2
%             'FL46_170131_11_31_33',... 3
%             'FL65_170426_11_56_56',... 4
%             'FL65_170607_12_50_17',... 5
%             'FL74_170822_10_38_13',... 6
%             'FL75_170822_11_39_40',... 7
%             'FL_S_170905_11_16_10'};  %8

% Choose dataset
dataset_ix = 8;
fname = datasets{dataset_ix};
disp(fname)

fname_fibreangles = [basedir,fname,'/processed/fibre_direction.mat'];
load(fname_fibreangles,'vector_mean');

fname_SNR = [basedir,fname,'/processed/SNR.mat'];
load(fname_SNR,'SNR_thresh');

fname_corr = [basedir,fname,'/processed/corr_histograms.mat'];
load(fname_corr,'C_inter');
rho_min = prctile(C_inter,95)

load([basedir,fname,'/',fname,'.mat'],'Numb_patches','Numb_trials')

% Get time information
MatrixTime = [];
load([basedir,fname,'/',fname,'_time.mat'])

% Default smoothing window
smooth_win_s = [];

Cn_all = cell(Numb_patches,1);
dFF_axons_all = cell(Numb_patches,1);
Ain_axons_all = cell(Numb_patches,1);
dFF_rois_all = cell(Numb_patches,1);
Ain_rois_all = cell(Numb_patches,1);
ix_axons_to_rois_all = cell(Numb_patches,1);
axon_ids_all = cell(Numb_patches,1);
time_rois_all = cell(Numb_patches,1);
time_axons_all = cell(Numb_patches,1);

% Check if file already exists
if isfile([basedir,fname,'/processed/',fname,'_GroupedData.mat'])
    error('Data has already been grouped. Proceeding will cause this file to be overwritten.')
end

%% Choose patch number
for patch_no = 1:Numb_patches
    
    %
    disp([num2str(patch_no),' / ',num2str(Numb_patches)])
     
    load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
    %load([basedir,fname,'\raw\Patch',sprintf('%03d',patch_no),'.mat'])
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
    dFF = get_dFF(Ain,Y,acquisition_rate,smooth_win_s);

    % Plot all dFF
    figure, hold on
    for k = 1:size(dFF,1)
        plot((1:size(dFF,2))/acquisition_rate,dFF(k,:)+k)
    end
    
    % Manually delete ROIs due to slow drift
    [Ain,dFF] = delete_ROIs_manually(Ain,dFF,acquisition_rate);
    
    % Calculate time for rois
    time_rois = zeros(size(dFF));
    for roi = 1:size(dFF,1)
        Ain_temp = reshape(Ain(:,roi),d1,d2);
        c = regionprops(Ain_temp,'centroid'); c = round(c.Centroid);
        time_rois(roi,:) = get_times(MatrixTime_patch{patch_no}(c(2),c(1)),MatrixTime(end,end),Numb_cycle,Time_between_trial,Numb_trials);
    end

    % Group axons
    [Ain_axons,ix_axons_to_rois,axon_ids] = ...
        get_axon_grouping(Ain,Y,[d1,d2],acquisition_rate,smooth_win_s,Pixel_size,vector_mean,[],rho_min,[],[],1);
        
    plot_grouped_rois(Ain_axons,Cn,dFF,ix_axons_to_rois,acquisition_rate);
    disp('Press enter to continue')
    pause
    
    % Save figures in png
        figure(1), saveas(gcf,[basedir,fname,'/figs/spatial_grouped_Patch',sprintf('%03d',patch_no),'.png']);
        figure(2), saveas(gcf,[basedir,fname,'/figs/dFF_grouped_Patch',sprintf('%03d',patch_no),'.png']);
        figure(3), saveas(gcf,[basedir,fname,'/figs/spatial_loners_Patch',sprintf('%03d',patch_no),'.png']);
        figure(4), saveas(gcf,[basedir,fname,'/figs/dFF_loners_Patch',sprintf('%03d',patch_no),'.png']);
    close all
    
    % Save results
    Ain_rois = Ain; dFF_rois = dFF;
    dFF_axons = get_dFF(Ain_axons,Y,acquisition_rate,smooth_win_s);
    
    % Calculate time for axons
    time_axons = zeros(size(dFF_axons));
    for roi = 1:size(dFF_axons,1)
        Ain_temp = reshape(Ain_axons(:,roi),d1,d2);
        c = regionprops(Ain_temp,'centroid'); c = round(c.Centroid);
        time_axons(roi,:) = get_times(MatrixTime_patch{patch_no}(c(2),c(1)),MatrixTime(end,end),Numb_cycle,Time_between_trial,Numb_trials);
    end
    
    % To save all data later
    Cn_all{patch_no} = Cn;
    dFF_axons_all{patch_no} = dFF_axons;
    Ain_axons_all{patch_no} = Ain_axons;
    dFF_rois_all{patch_no} = dFF_rois;
    Ain_rois_all{patch_no} = Ain_rois;
    ix_axons_to_rois_all{patch_no} = ix_axons_to_rois;
    axon_ids_all{patch_no} = axon_ids;
	time_rois_all{patch_no} = time_rois;
    time_axons_all{patch_no} = time_axons;

end


%% Save all data

Cn = Cn_all;
dFF_axons = dFF_axons_all;
Ain_axons = Ain_axons_all;
dFF_rois = dFF_rois_all;
Ain_rois = Ain_rois_all;
ix_axons_to_rois = ix_axons_to_rois_all;
axon_ids = axon_ids_all;
time_rois = time_rois_all;
time_axons = time_axons_all;

save([basedir,fname,'/processed/',fname,'_GroupedData.mat'],...
    'Cn','dFF_axons','Ain_axons','dFF_rois','Ain_rois',...
    'time_axons','time_rois','ix_axons_to_rois','axon_ids','acquisition_rate');

disp('All data successfully saved.')

