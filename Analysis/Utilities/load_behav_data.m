% Load behavioural variables then smooths and interpolates to same rate
%
% Required input:
%    dataset_ix     Dataset number (1-8)
%    time           Interpolation time - if empty, just loads values w/out
%                     interpolation
%    smooth_win_s   Smoothing window for speed
% 
% Output:
%    whisk_angle         whisking angle
%    whisk_set_point     whisker set point
%    whisk_amp           Whisking amplitude
%    loco_smooth         Smoothed wheel motion index
%    Speed_smooth        Smoothed wheel MI
%    whisk_time          Native (non-interp) time of whisking
%    loco_time           Native (non_interp) time of WMI  
%    Speed_time          Native (non_interp) time of speed 

function [whisk_angle,whisk_set_point,whisk_amp,loco_smooth,speed_smooth,whisk_time,loco_time,speed_time] = load_behav_data(dataset_ix,time,smooth_win_s)

    if nargin < 2 
        time = [];
    end

    if nargin < 3 || isempty(smooth_win_s)
        smooth_win_s = 0.2;
    end
    
    define_dirs;
    
    % Whisking - get all vars except phase
    [whisk_time,whisk_angle,whisk_set_point,whisk_amp,~] = get_whisking_vars(dataset_ix);
    
    % Load wheel MI for locomotion
    fname = datasets{dataset_ix};
    load([basedir,fname,'/',fname,'_MIwheel.mat']);
    
    % Load encoder speed and wheelMI
    load([basedir,fname,'/',fname,'.mat'],'SpeedTimeMatrix','SpeedDataMatrix');

    [whisk_angle,whisk_set_point,whisk_amp,loco_smooth,speed_smooth,whisk_time,loco_time,speed_time] = ...
        convert_behav_vars(time,smooth_win_s,whisk_time,whisk_angle,whisk_set_point,whisk_amp,wheel_MI,SpeedTimeMatrix,SpeedDataMatrix);
    
    