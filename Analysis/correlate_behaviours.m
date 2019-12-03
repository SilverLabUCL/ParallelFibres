% This script groups ROIs onto putative axons, separately for each patch in
% each experiment

clear all; clc; close all
addpath('Utilities/')

basedir = '~/Documents/ParallelFibres/Data/';
datasets = {'FL87_180501_11_03_09',...  1
            'FL87_180501_10_47_25',...  2
            'FL87_180501_10_36_14',...  3 
            'FL87_180220_10_38_55',...  4
            'FL77_180213_10_46_41',...  5
            'FL_S_170906_11_26_25',...  6
            'FL_S_170905_10_40_52',...  7            
            'FL45_170125_14_47_04',...  8
             ...%%
             'FL95_180425_10_53_40',... 9
             'FL87_180413_11_00_55',... 10
             'FL104_180725_10_42_37',...11
             'FL87_180117_11_23_20',... 12
             'FL106_180807_10_52_25',...13
             'FL_S_171109_14_54_34',... 14
             'FL_S_171109_15_19_52',... 15
             ...%%
             'FL92_180228_11_10_48',... 16
             'FL92_180228_11_18_24',... 17
             };

% Choose dataset
dataset_ix = 16;
fname = datasets{dataset_ix};
disp(fname)

[whisk_angle,whisk_set_point,whisk_amp,Speed_smooth,whisk_time,Speed_time] ...
    = load_behav_data(dataset_ix,[],[]);

disp('whisk set point and whisk amp')
corr(whisk_set_point,whisk_angle)

disp('speed and whisk set point')
ix = find(Speed_time >= whisk_time(1) & Speed_time <= whisk_time(end));
corr(Speed_smooth(ix),interp1(whisk_time,whisk_set_point,Speed_time(ix)))

disp('speed and whisk amp')
corr(Speed_smooth(ix),interp1(whisk_time,whisk_amp,Speed_time(ix)))

if isfile([basedir,fname,'/',fname,'_Pupil.mat'])
    load([basedir,fname,'/',fname,'_Pupil.mat'])
    
    Pupil_time = Pupil(:,1); Pupil_smooth = Pupil(:,2);
    figure, plot(Pupil_time,Pupil_smooth,'Color',[.7,.7,.7]), hold on
    
    % Filter at 1s (see Reimer 2014)
    fs = length(Pupil_time)/((Pupil_time(end)-Pupil_time(1))/1000);
    fc = 1; cutoff = fc/(fs/2);
    [b,a] = butter(4,cutoff);
    Pupil_smooth = filter(b,a,Pupil_smooth);
    Pupil_smooth = filter(b,a,Pupil_smooth(end:-1:1));
    Pupil_smooth = Pupil_smooth(end:-1:1);

    % Remove boundary effects
    Pupil_time(1:round(1/cutoff)) = [];
    Pupil_smooth(1:round(1/cutoff)) = [];
    Pupil_time(end-round(1/cutoff):end) = [];
    Pupil_smooth(end-round(1/cutoff):end) = [];
    plot(Pupil_time,Pupil_smooth,'r','LineWidth',1)
    
    % Convert to seconds
    Pupil_time = Pupil_time/1000;
    
    ix = find(Pupil_time >= whisk_time(1) & Pupil_time <= whisk_time(end));
    disp('Pupil and whisk amp')
    corr(Pupil_smooth(ix),interp1(whisk_time,whisk_amp,Pupil_time(ix)))
    disp('Pupil and whisk set point')
    corr(Pupil_smooth(ix),interp1(whisk_time,whisk_set_point,Pupil_time(ix)))
    disp('Pupil and speed')
    ix = find(Pupil_time >= Speed_time(1) & Pupil_time <= Speed_time(end));
    corr(Pupil_smooth(ix),interp1(Speed_time,Speed_smooth,Pupil_time(ix)))
    
end


%% Filter at 1s (see Reimer 2014)
    fs = length(Speed_time)/((Speed_time(end)-Speed_time(1)));
    fc = 1; cutoff = fc/(fs/2);
    [b,a] = butter(4,cutoff,'high');
    temp = filter(b,a,Speed_smooth);
    temp = filter(b,a,temp(end:-1:1));
    temp = temp(end:-1:1);
    temp = envelope(temp);
   
    figure, plot(Speed_time,zscore(Speed_smooth))
    hold on, plot(Speed_time,zscore(temp))

