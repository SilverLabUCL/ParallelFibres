% Old code - see group_rois
clear all; clc

patch_no = 1;
acq_rate = 15;

basedir = '~/Documents/FredPF/raw/FL87_180501_11_03_09/';
fname_rawdata = strcat(basedir,'Exported TIFFs/GreenChannel_rectangle_',sprintf('%03d',patch_no),'.tif');
fname_save = strcat(basedir,'Processed/Patch',sprintf('%03d',patch_no),'.mat');
fname_fibreangles = [basedir,'Processed/fibre_angles.mat'];

% Load raw data
info = imfinfo(fname_rawdata);

T = numel(info);
d2 = info(1).Width;
d1 = info(1).Height;
Y = zeros(d1,d2,T);
for t = 1:T
    A = imread(fname_rawdata, t, 'Info', info);
    Y(:,:,t) = A;
end
Y = reshape(Y, d1*d2, []);

% Remove "dark noise"
Y = Y - min(Y(:));

um_per_pix = 0.341;

%% Use CNMF initialization to estimate initial spatial filters
[Ain,Cn] = detect_ROIs(Y, [d1,d2], 1, 5, 1);

%% Before removing ROIs, check neuropil
[~,dFF_neuropil2] = get_neuropil_byhand(Ain,Y,[d1,d2],acq_rate);

%% Trim ROIs to remove 
Ain = trim_ROIs(Ain,[d1,d2]);

%% Remove low SNR ROIs
[Ain,SNR,A_bad,~] = remove_bad_cells(Ain,Y,[d1,d2],acq_rate,6); 

%% Calculate dFF
dFF = get_dFF(Ain,Y,acq_rate);

%%

figure, hold on
for k = 1:size(dFF,1)
    plot((1:size(dFF,2))/acq_rate,dFF(k,:)+k)
end

%% Group axons

load(fname_fibreangles,'vector_mean','angle_std');

[Ain_axons,ix_axons_to_rois_new,axon_ids_new] = get_axon_grouping(Ain,Y,[d1,d2],acq_rate,um_per_pix,.5,vector_mean,0);
dFF_axons = get_dFF(Ain_axons,Y,acq_rate);

%% Plot results
figure, plot_grouped_rois(Ain_axons,Cn,ix_axons_to_rois_new,.95,true)

figure, plot_loners(Ain_axons,Cn,ix_axons_to_rois_new,.95,true)

