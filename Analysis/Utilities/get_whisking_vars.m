% Load whisking variables
%
% Required input:
%    dataset_ix   Dataset number (1-13)
% 
% Output:
%    dlc_whisk_time      Time returned from DLC (in s)
%    whisk_angle_filt    Filtered whisking angle
%    whisk_set_point     Whisker set point
%    whisk_amp           Whisking amplitude
%    whisk_phase         Whisking phase (not used)


function [dlc_whisk_time,whisk_angle_filt,whisk_set_point,whisk_amp,whisk_phase] = get_whisking_vars(dataset_ix)

    define_dirs;
    
    % Load DLC wihsking info
    fname = datasets{dataset_ix};
    load([basedir,fname,'/',fname,'.mat'],'dlc_whisk_angle','dlc_whisk_time')
    
    [dlc_whisk_time,whisk_angle_filt,whisk_set_point,whisk_amp,whisk_phase] ...
        = convert_dlc_to_whisk_vars(dlc_whisk_time,dlc_whisk_angle);
    
    

