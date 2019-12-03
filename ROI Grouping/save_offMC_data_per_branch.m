% This script converts large (~4GB) mat files that Fred gave me into: 
%     1. a small mat file (~15MB) with metadata and behavioural data
%     2. a larger mat file (~140MB) for raw post-offMC data for each patch

clear all; clc

basedir = '~/Documents/ParallelFibres/Data/';
%'~/Documents/FredPF/raw/offMC/';
fname = 'FL76_170913_10_57_06';
%'FL77_180113_10_58_50';
%'FL75_170912_10_33_29';


% Only arousal data
%'FL92_180228_11_10_48'
%'FL92_180228_11_18_24'

% Data later added by Fred
%'FL95_180425_10_53_40';
%'FL87_180413_11_00_55';
%'FL104_180725_10_42_37';
%'FL87_180117_11_23_20';
%'FL106_180807_10_52_25';
%'FL_S_171109_14_54_34'; % check w fred if this has pupil !!!
%'FL_S_171109_15_19_52'; % check w fred if this has pupil !!!

% Original datasets
%'FL87_180501_11_03_09';
%'FL87_180501_10_47_25';
%'FL87_180501_10_36_14';
%'FL87_180220_10_38_55';
%'FL77_180213_10_46_41';
%'FL45_170125_14_47_04';
%'FL_S_170906_11_26_25';
%'FL_S_170905_10_40_52';

disp([basedir,fname])

%% Load data
load([basedir,'Crus1_patches_',fname,'_driftcorrected.mat'])
%load([basedir,'Crus1_patches_',fname,'_arousaldata.mat'])

disp('Data has successfully been loaded.')

%% Save metadata, behavioural data, etc.

%save([basedir,fname,'/',fname,'.mat'],'acquisition_rate','Pixel_size','Axon_dFF','TimeAxon','Axons_coordinates','dlc_whisk_angle','dlc_whisk_time','MI_facepad','Numb_frames','Numb_patches','Numb_trials','Patch_coordinates','pia','SpeedDataMatrix','SpeedTimeMatrix');
save([basedir,fname,'/',fname,'.mat'],'acquisition_rate','Pixel_size','dlc_whisk_angle','dlc_whisk_time','Numb_frames','Numb_patches','Numb_trials','Patch_coordinates','pia','SpeedDataMatrix','SpeedTimeMatrix');
%save([basedir,fname,'/',fname,'.mat'],'dlc_whisk_angle','dlc_whisk_time','SpeedDataMatrix','SpeedTimeMatrix');

if exist('Whiskers_time_0')
    save([basedir,fname,'/',fname,'.mat'],'Whiskers_angle_0','Whiskers_time_0','-append');
end
if exist('Pupil')
    save([basedir,fname,'/',fname,'_Pupil.mat'],'Pupil');
end
if exist('wheel_MI')
    save([basedir,fname,'/',fname,'.mat'],'wheel_MI','-append');
    save([basedir,fname,'/',fname,'_MIwheel.mat'],'wheel_MI');
end

if exist('MatrixTime')
    save([basedir,fname,'/',fname,'_time.mat'],'MatrixTime','MatrixTime_patch','Time_between_trial');
end
disp('Metadata has successfully been saved.')

%% Covert patch data into uint16

clearvars -except patches acquisition_rate Pixel_size fname basedir

[d1,d2,T] = size(patches{1});

darknoise = zeros(1,size(patches,2));
Y_all = cell(size(patches));
for p = 1:size(patches,2)
    p
    Y = reshape(patches{p},d1*d2,T);
    % First check that Y doesn't have any significant imaginary parts
    if imag(sum(sum((Y)))) > 10^-5
        error
    end
    Y = real(Y);
    % Then check that Y doesn't go out of range of uint16
    darknoise(p) = min(min(Y));
    if max(max(Y)) > 65535 || darknoise(p) < 0
        error
    end
    Y_all{p} = uint16(Y);
end

disp('Conversion of patch data is complete.')

%% Remove dark noise and save

darknoise = uint16(min(darknoise));

for p = 1:size(patches,2)
    p
    Y = Y_all{p} - darknoise;
    save([basedir,fname,'/raw/Patch',sprintf('%03d',p),'.mat'],'Y','d1','d2','acquisition_rate','Pixel_size')
    
end
disp('Patch data has been successfully saved.')