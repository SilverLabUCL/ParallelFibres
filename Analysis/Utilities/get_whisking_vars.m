% Load DLC whisking angle and extract whisking variables
% Version for original datasets - don't need to manually load
%
% Required input:
%    dataset_ix   Dataset number (1-8)
% 
% Output:
%    dlc_whisk_time      time for whisking data
%    whisk_angle_filt    filtered whisking angle
%    whisk_set_point     whisker set point
%    whisk_amp           Whisking amplitude
%    whisk_phase  


function [dlc_whisk_time,whisk_angle_filt,whisk_set_point,whisk_amp,whisk_phase] = get_whisking_vars(dataset_ix)

    define_dirs;
    
    % Load DLC wihsking info
    fname = datasets{dataset_ix};
    load([basedir,fname,'/',fname,'.mat'],'dlc_whisk_angle','dlc_whisk_time')
    
    [dlc_whisk_time,whisk_angle_filt,whisk_set_point,whisk_amp,whisk_phase] ...
        = convert_dlc_to_whisk_vars(dlc_whisk_time,dlc_whisk_angle);
    
    

