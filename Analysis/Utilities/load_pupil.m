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

function [Pupil_smooth,Pupil_time] = load_pupil(dataset_ix,time,smooth_win_s)

    if nargin < 2 
        time = [];
    end

    if nargin < 3 || isempty(smooth_win_s)
        smooth_win_s = 0.2;
    end
    
    define_dirs;
    
    fname = datasets{dataset_ix};
    
    Pupil_time = [];
    Pupil_smooth = [];
    
    try
        load([basedir,fname,'/',fname,'_Pupil.mat']);

        Pupil_time = Pupil(:,1) / 1000;
        dt = mean(diff(Pupil_time));
        Pupil_smooth = smoothdata(Pupil(:,2),'gaussian',[round(.2/dt) 0] *2);

        % Interpolate all data to 
        if ~isempty(time)
            Pupil_smooth = interp1(Pupil_time,Pupil_smooth,time);
        end
    end