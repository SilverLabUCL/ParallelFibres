
clear all; clc
addpath('CNMFE_code/')

patch_no = 2;

basedir = '~/Documents/FredPF/raw/offMC/';
fname = 'FL87_180501_11_03_09';
load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])

fname_fibreangles = ['~/Documents/FredPF/raw/',fname,'/Processed/fibre_angles.mat'];

Y = double(Y);

%% Use CNMF initialization to estimate initial spatial filters
[Ain,Cn] = detect_ROIs(Y, [d1,d2], 1, 5);

%% Before removing ROIs, check neuropil
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

make_video_dFF(Ain,Y,[d1,d2],acquisition_rate,[basedir,fname,'/Processed/Patch',sprintf('%03d',patch_no)])

%% Get inter-bouton distance

d = get_interbouton_dist(Ain,[d1,d2],ix_axons_to_rois,Pixel_size);
[mean(d),std(d)]

%% Save inter-bouton distance, 

