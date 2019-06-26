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
dataset_ix = 8;

% Choose patch number
patch_no = 3;

fname = datasets{dataset_ix};
disp(fname)

load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])

fname_fibreangles = [basedir,fname,'/processed/fibre_direction.mat'];
fname_SNR = [basedir,fname,'/processed/SNR.mat'];
fname_corr = [basedir,fname,'/processed/corr_histograms.mat'];

Y = double(Y);

Cn = correlation_image(Y, [1,2], d1,d2);
figure, imagesc(Cn)

%% Use CNMF initialization to estimate initial spatial filters
[Ain,~] = detect_ROIs(Y, [d1,d2]);

%% Remove low SNR ROIs
load(fname_SNR,'SNR_thresh');
[Ain,~,~,~] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,SNR_thresh); 

%% Trim ROIs to remove 
Ain = trim_ROIs(Ain,[d1,d2]);

%% Remove low SNR ROIs again
[Ain,SNR,A_bad,SNR_bad] = remove_bad_cells(Ain,Y,[d1,d2],acquisition_rate,SNR_thresh); 

%% Calculate dFF
dFF = get_dFF(Ain,Y,acquisition_rate);

% Plot all dFF
figure, hold on
for k = 1:size(dFF,1)
    plot((1:size(dFF,2))/acquisition_rate,dFF(k,:)+k)
end

%% Group axons
load(fname_fibreangles,'vector_mean');
load(fname_corr,'C_inter');
rho_min = prctile(C_inter,95);

%get_axon_grouping(Ain,Y,dims,acq_rate,um_per_pix,fibre_vec,...
%   angle_std,rho_min,var_ratio_max,pix_defl_max,manual_check)

[Ain_axons,ix_axons_to_rois,axon_ids_new] = ...
    get_axon_grouping(Ain,Y,[d1,d2],acquisition_rate,Pixel_size,vector_mean,...
    [],rho_min,[],[],1);

dFF_axons = get_dFF(Ain_axons,Y,acquisition_rate);

%% Plot results
plot_grouped_rois(Ain_axons,Cn,dFF,ix_axons_to_rois,acquisition_rate);

%% Save video if want to check motion correction

make_video_dFF(Ain,Y,[d1,d2],acquisition_rate)%,[basedir,fname,'/Processed/Patch',sprintf('%03d',patch_no)])



