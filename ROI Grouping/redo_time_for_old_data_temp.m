% This script groups ROIs onto putative axons, separately for each patch in
% each experiment

clear all; clc; close all
addpath('From CNMF_E/')
addpath('Utilities/')

basedir = '~/Documents/ParallelFibres/Data/';
%basedir = 'C:\Users\SilverLab\Documents\Alex\ParallelFibres\Data\';
datasets = {'FL87_180501_11_03_09',...  1
            'FL87_180501_10_47_25',...  2
            'FL87_180501_10_36_14',...  3 
            'FL87_180220_10_38_55',...  4
            'FL77_180213_10_46_41',...  5
            'FL_S_170906_11_26_25',...  6
            'FL_S_170905_10_40_52',...  7            
            'FL45_170125_14_47_04'};

% Choose dataset
dataset_ix = 7;
fname = datasets{dataset_ix};
disp(fname)

load([basedir,fname,'/',fname,'.mat'],'Numb_patches','Numb_trials')

% time
load([basedir,fname,'/',fname,'_time.mat'])
load([basedir,fname,'/processed/',fname,'_GroupedData.mat'])

smooth_win_s = [];


%% Choose patch number

time_rois = cell(Numb_patches,1);
time_axons = cell(Numb_patches,1);

for patch_no = 1:Numb_patches
    
    disp([num2str(patch_no),' / ',num2str(Numb_patches)])
     
    load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
    Y = double(Y);
    
    % Calculate dFF
    dFF_rois{patch_no} = get_dFF(Ain_rois{patch_no},Y,acquisition_rate,smooth_win_s);
    dFF_axons{patch_no} = get_dFF(Ain_axons{patch_no},Y,acquisition_rate,smooth_win_s);
    
    % Calculate time
    time_rois{patch_no} = zeros(size(dFF_rois{patch_no}));
    for roi = 1:size(dFF_rois{patch_no},1)
        Ain_temp = reshape(Ain_rois{patch_no}(:,roi),d1,d2);
        c = regionprops(Ain_temp,'centroid'); c = round(c.Centroid);
        time_rois{patch_no}(roi,:) = get_times(MatrixTime_patch{patch_no}(c(2),c(1)),MatrixTime(end,end),Numb_cycle,Time_between_trial,Numb_trials);
    end
    
    time_axons{patch_no} = zeros(size(dFF_axons{patch_no}));
    for roi = 1:size(dFF_axons{patch_no},1)
        Ain_temp = reshape(Ain_axons{patch_no}(:,roi),d1,d2);
        c = regionprops(Ain_temp,'centroid'); c = round(c.Centroid);
        time_axons{patch_no}(roi,:) = get_times(MatrixTime_patch{patch_no}(c(2),c(1)),MatrixTime(end,end),Numb_cycle,Time_between_trial,Numb_trials);
    end

end


%% Save all data

save([basedir,fname,'/processed/',fname,'_GroupedData.mat'],...
    'Cn','dFF_axons','Ain_axons','dFF_rois','Ain_rois','ix_axons_to_rois',...
    'time_rois','time_axons','axon_ids','acquisition_rate');

disp('All data successfully saved.')

