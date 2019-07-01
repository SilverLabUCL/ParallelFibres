

%% Test
clear all; clc; close all

basedir = '~/Documents/ParallelFibres/Data/';

datasets = {'FL87_180501_11_03_09',...  1
                'FL87_180501_10_47_25',...  2
                'FL87_180501_10_36_14',...  3
                'FL87_180220_10_38_55',...  4
                'FL77_180213_10_46_41',...  5
                'FL_S_170906_11_26_25',...  6
                'FL_S_170905_10_40_52',...  7
                'FL45_170125_14_47_04'}; %  8


            
dataset_ix = 2;

fname = datasets{dataset_ix};

load([basedir,fname,'/',fname,'.mat'])

load([basedir,fname,'/',fname,'_MIwheel.mat'])

[whisk_angle,whisk_sp,whisk_amp,whisk_phase] = get_whisking_vars(dlc_whisk_angle,dlc_whisk_time);

figure, plot(whisk_sp, whisk_amp,'.')

figure, plot(dlc_whisk_time,zscore(whisk_angle))
hold on, plot(dlc_whisk_time,zscore(whisk_sp),'r','LineWidth',2)

fs_MI = length(wheel_MI(:,2)) / ((wheel_MI(end,2)-wheel_MI(1,2)) * 1/1000);
fs_enc = length(SpeedTimeMatrix) / ((SpeedTimeMatrix(end)-SpeedTimeMatrix(1)) * 1/1000);

temp_MI = smoothdata(wheel_MI(:,1),'Gaussian',[0.200 * fs_MI 0]*2);
temp_enc = smoothdata(SpeedDataMatrix,'Gaussian',[0.200 * fs_enc 0]*2);

plot(wheel_MI(:,2),zscore(temp_MI)-5,'k','LineWidth',2)
plot(SpeedTimeMatrix-SpeedTimeMatrix(1),zscore(temp_enc)-10,':k','LineWidth',2)
%%
figure, plot(wheel_MI(:,2),zscore(temp_MI))
hold on, plot(SpeedTimeMatrix,zscore(temp_enc))

temp_enc = interp1(SpeedTimeMatrix,temp_enc,wheel_MI(:,2));
figure, plot(temp_enc,temp_MI,'.')