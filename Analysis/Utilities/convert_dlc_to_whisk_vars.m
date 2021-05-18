% Extract and filter whisking variables from DeepLabCut extracted tracking
%
% Input:
%    dlc_whisk_time      Time returned from DLC (in ms)
%    dlc_whisk_angle     Whisker position/angle from DLC
% 
% Output:
%    dlc_whisk_time      Time returned from DLC (in s)
%    whisk_angle_filt    Filtered whisking angle
%    whisk_set_point     Whisker set point
%    whisk_amp           Whisking amplitude
%    whisk_phase         Whisking phase (not used)

function [dlc_whisk_time,whisk_angle_filt,whisk_set_point,whisk_amp,whisk_phase]  = convert_dlc_to_whisk_vars(dlc_whisk_time,dlc_whisk_angle)

    % Convert to rad and set so protraction is positive
    dlc_whisk_angle = pi - deg2rad(dlc_whisk_angle);

    % Get whisking resting point
    whisk_rest = mode(dlc_whisk_angle);
    dlc_whisk_angle = dlc_whisk_angle - whisk_rest;

    % Sampling rate
    fs = length(dlc_whisk_time)/((dlc_whisk_time(end)-dlc_whisk_time(1))/1000);

    % Filter to get denoised whisker angle
    fc = 30;
    [b,a] = butter(4,fc/(fs/2));
    whisk_angle_filt = filter(b,a,dlc_whisk_angle);

    whisk_angle_filt = filter(b,a,whisk_angle_filt(end:-1:1));

    whisk_angle_filt = whisk_angle_filt(end:-1:1);

    % To get set point, median filter by 500 ms
    whisk_set_point = smoothdata(whisk_angle_filt,'gaussian',round(0.5* fs));

    % To get whisking amplitude and phase,
    fc2 = 30; fc1 = 8;
    [b,a] = butter(4,[fc1,fc2]/(fs/2));
    whisk_angle_bp = hilbert(filter(b,a,dlc_whisk_angle));
    
    whisk_phase = angle(whisk_angle_bp);
    whisk_amp = abs(whisk_angle_bp);
    
    dlc_whisk_time = dlc_whisk_time / 1000;
    
    % Convert to degrees
    whisk_set_point = rad2deg(whisk_set_point);
    whisk_amp = rad2deg(whisk_amp);