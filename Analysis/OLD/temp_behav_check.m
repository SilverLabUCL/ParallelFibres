

%% Test
clear all; clc; 

            
dataset_ix = 2;

% Load data
[dFF,time,acquisition_rate] = load_data(dataset_ix,0);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time,0.2);

%ix = find(~isnan(whisk_amp));
C = corr(whisk_amp,dFF');
figure, histogram(C)
%% Load meta information

% Load whisking info
[whisk_time,whisk_angle,whisk_set_point,whisk_amp,~] = get_whisking_vars(dataset_ix);




load([basedir,fname,'/',fname,'.mat'])
load([basedir,fname,'/',fname,'_MIwheel.mat'])

[whisk_angle,whisk_sp,whisk_amp,whisk_phase] = get_whisking_vars(dlc_whisk_angle,dlc_whisk_time);
[whisk_time,whisk_angle,whisk_set_point,whisk_amp,~] = get_whisking_vars(dataset_ix);

figure, plot(whisk_sp, whisk_amp,'.')

figure, plot(dlc_whisk_time,zscore(whisk_angle))
hold on, plot(dlc_whisk_time,zscore(whisk_sp),'r','LineWidth',2)
`
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