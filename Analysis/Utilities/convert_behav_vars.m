% Load DLC whisking angle and extract whisking variables
%
%
% Required input:
%    dataset_ix     Dataset number (1-8)
%    time           Interpolation time - if empty, just loads values w/out
%                     interpolation
%    smooth_win_s   Smoothing window for speed/loco
% 
% Output:
%    whisk_angle         whisking angle
%    whisk_set_point     whisker set point
%    whisk_amp           Whisking amplitude
%    loco_smooth        Smoothed wheel MI (locomotion)
%    whisk_time          Native (non-interp) time of whisking
%    loco_time          Native (non_interp) time of locomotion 


function [whisk_angle,whisk_set_point,whisk_amp,loco_smooth,speed_smooth,whisk_time,loco_time,speed_time] = convert_behav_vars(time,smooth_win_s,whisk_time,whisk_angle,whisk_set_point,whisk_amp,wheel_MI,SpeedTimeMatrix,SpeedDataMatrix)

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
    
    % Smooth loco with a half-Gaussian
    loco_time = wheel_MI(:,2) / 1000;
    dt_loco = mean(diff(loco_time));
    loco_smooth = smoothdata(wheel_MI(:,1),'gaussian',[round(smooth_win_s/dt_loco) 0] *2);
    
    % Convert to units of standard deviations
    loco_smooth = zscore(loco_smooth);
    
    % Smooth speed with a half-Gaussian
    speed_time = SpeedTimeMatrix/1000;
    dt_speed = mean(diff(speed_time));
    speed_smooth = smoothdata(SpeedDataMatrix,'gaussian',[round(smooth_win_s/dt_speed) 0] *2);
    
    % Interpolate all data to functional time
    if ~isempty(time)
        loco_smooth = interp1(loco_time,loco_smooth,time);
        speed_smooth = interp1(speed_time,speed_smooth,time);
        whisk_angle = interp1(whisk_time,whisk_angle,time);
        whisk_set_point = interp1(whisk_time,whisk_set_point,time);
        whisk_amp = interp1(whisk_time,whisk_amp,time);
    end