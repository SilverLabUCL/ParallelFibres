% This script groups ROIs onto putative axons, separately for each patch in
% each experiment

clear all; clc
addpath('From CNMF_E/')
addpath('Utilities/')

basedir = '~/Documents/FredPF/raw/offMC/';
datasets = {'FL87_180501_11_03_09',...  1
            'FL87_180501_10_47_25',...  2
            'FL87_180501_10_36_14',...  3 
            'FL87_180220_10_38_55',...  4
            'FL77_180213_10_46_41',...  5
            'FL_S_170906_11_26_25',...  6
            'FL_S_170905_10_40_52',...  7            
            'FL45_170125_14_47_04'}; %  8

% Choose dataset
dataset_ix = 4;

% Choose patch number
patch_no = 16;

fname = datasets{dataset_ix};
disp(fname)

load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])

%fname_fibreangles = ['~/Documents/FredPF/raw/',fname,'/Processed/fibre_angles.mat'];
fname_fibreangles = ['~/Documents/FredPF/raw/FL87_180501_11_03_09/Processed/fibre_angles.mat'];

Y = double(Y);

%% Use CNMF initialization to estimate initial spatial filters
[Ain,Cn] = detect_ROIs(Y, [d1,d2], 1, 5);

%% Check neuropil
[~,dFF_neuropil2] = get_neuropil_byhand(Ain,Y,[d1,d2],acquisition_rate);

%% Remove low SNR ROIs
[Ain,~,~,~] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,6); 


%% Trim ROIs to remove 
Ain = trim_ROIs(Ain,[d1,d2]);

%% Remove low SNR ROIs again
[Ain,SNR,A_bad,~] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,6); 

%% Calculate dFF
dFF = get_dFF(Ain,Y,acquisition_rate);

%% Plot all dFF
figure, hold on
for k = 1:size(dFF,1)
    plot((1:size(dFF,2))/acquisition_rate,dFF(k,:)+k)
end

%% Group axons
load(fname_fibreangles,'vector_mean');

%[Ain_axons,ix_axons_to_rois_new,axon_ids_new] = get_axon_grouping(Ain,Y,[d1,d2],acquisition_rate,Pixel_size,.5,vector_mean,1);
[Ain_axons,ix_axons_to_rois,axon_ids_new] = ...
    get_axon_grouping(Ain,Y,[d1,d2],acquisition_rate,Pixel_size,vector_mean,...
    [],[],[],[],1);

dFF_axons = get_dFF(Ain_axons,Y,acquisition_rate);

%% Plot results
figure, plot_grouped_rois(Ain_axons,Cn,ix_axons_to_rois,.95,true)


figure, plot_loners(Ain_axons,Cn,ix_axons_to_rois,.95,true)

%% Save video if want to check motion correction

make_video_dFF(Ain,Y,[d1,d2],acquisition_rate)%,[basedir,fname,'/Processed/Patch',sprintf('%03d',patch_no)])



