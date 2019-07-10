% Load DLC whisking angle and extract whisking variables
%
%
% Required input:
%    dataset_ix   Dataset number (1-8)
%    time at which to interpolate
% 
% Output:
%    dlc_whisk_time      time for whisking data
%    whisk_angle_filt    filtered whisking angle
%    whisk_set_point     whisker set point
%    whisk_amp           Whisking amplitude
%    whisk_phase  

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
    
    % Correct double frames by taking average
    ix_to_correct = find(diff(wheel_MI(:,2))==0);
    for ix = ix_to_correct
        wheel_MI(ix,1) = (wheel_MI(ix,1)+wheel_MI(ix+1,1))/2;
        wheel_MI(ix+1,:) = [];
    end
    
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
    
    Speed_time = wheel_MI(:,2) / 1000;
    dt_speed = mean(diff(Speed_time));
    Speed_smooth = smoothdata(wheel_MI(:,1),'gaussian',[round(smooth_win_s/dt_speed) 0] *2);
    
    % Interpolate all data to 
    if ~isempty(time)
        Speed_smooth = interp1(Speed_time,Speed_smooth,time);
        whisk_angle = interp1(whisk_time,whisk_angle,time);
        whisk_set_point = interp1(whisk_time,whisk_set_point,time);
        whisk_amp = interp1(whisk_time,whisk_amp,time);
    end