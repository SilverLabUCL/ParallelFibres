% Load DLC whisking angle and extract whisking variables
%
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
%    Speed_smooth        Smoothed wheel MI
%    whisk_time          Native (non-interp) time of whisking
%    Speed_time          Native (non_interp) time of speed 

function [whisk_angle,whisk_set_point,whisk_amp,Speed_smooth,whisk_time,Speed_time] = load_behav_data(dataset_ix,time,smooth_win_s)

    if nargin < 2 
        time = [];
    end

    if nargin < 3 || isempty(smooth_win_s)
        smooth_win_s = 0.2;
    end
    
    define_dirs;
    
    % Whisking    
    [whisk_time,whisk_angle,whisk_set_point,whisk_amp,~] = get_whisking_vars(dataset_ix);
    
    fname = datasets{dataset_ix};
    load([basedir,fname,'/',fname,'_MIwheel.mat']);
    
    % Correct double frames by taking average - wheel
    ix_to_correct = find(diff(wheel_MI(:,2))==0);
    for ix = ix_to_correct
        wheel_MI(ix,1) = (wheel_MI(ix,1)+wheel_MI(ix+1,1))/2;
        wheel_MI(ix+1,:) = [];
    end
    
    % Correct double frames by taking average - whisker pad
    ix_to_correct = find(diff(whisk_time)==0);
    for ix = ix_to_correct
        whisk_angle(ix,1) = (whisk_angle(ix,1)+whisk_angle(ix+1,1))/2;
        whisk_angle(ix+1,:) = [];
        whisk_set_point(ix,1) = (whisk_set_point(ix,1)+whisk_set_point(ix+1,1))/2;
        whisk_set_point(ix+1,:) = [];
        whisk_amp(ix,1) = (whisk_amp(ix,1)+whisk_amp(ix+1,1))/2;
        whisk_amp(ix+1,:) = [];
        whisk_time(ix+1) = [];
    end
    
    % Smooth speed
    Speed_time = wheel_MI(:,2) / 1000;
    dt_speed = mean(diff(Speed_time));
    Speed_smooth = smoothdata(wheel_MI(:,1),'gaussian',[round(smooth_win_s/dt_speed) 0] *2);
    
    % Convert to units of standard deviations
    Speed_smooth = zscore(Speed_smooth);
    
    % Interpolate all data to functional time
    if ~isempty(time)
        Speed_smooth = interp1(Speed_time,Speed_smooth,time);
        whisk_angle = interp1(whisk_time,whisk_angle,time);
        whisk_set_point = interp1(whisk_time,whisk_set_point,time);
        whisk_amp = interp1(whisk_time,whisk_amp,time);
    end