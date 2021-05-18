
% This function calculates correlations between activity and behavioural
% variables with a two-tailed shuffle test
%
% Input:
%    dFF                Fluorescence data
%    whisk_set_point    Whisker set point
%    whisk_amp          Whisker amplitude
%    loco               Wheel motion index (WMI)
%    speed              Locomotion speed
%    pupil              Pupil 
%    acquisition_rate   Acquisition rate
% 
% Output:       
%    C_wsp              Correlation between rois and whisker set point
%    p_wsp              p-values for Correlation between rois and whisker set point
%    C_wamp             Correlation between rois and whisker amplitude
%    p_wamp             p-values for Correlation between rois and whisker amplitude
%    C_loco             Correlation between rois and WMI
%    p_loco             p-values for Correlation between rois and WMI
%    C_speed            Correlation between rois and locomotion speed
%    p_speed            p-values for Correlation between rois and locomotion speed
%    C_pupil            Correlation between rois and pupil
%    p_pupil            p-values for Correlation between rois and pupil

function [C_wsp,p_wsp,C_wamp,p_wamp,C_loco,p_loco,C_speed,p_speed,C_pupil,p_pupil] = corr_sig(dFF,whisk_set_point,whisk_amp,loco,speed,pupil,acquisition_rate)

    [N,T] = size(dFF);
    num_reps = 1000;

    % True correlation
    C_wsp = corr(dFF',whisk_set_point,'rows','complete');
    C_wamp = corr(dFF',whisk_amp,'rows','complete');
    C_loco = corr(dFF',loco,'rows','complete');
    C_speed = corr(dFF',speed,'rows','complete');
    if isempty(pupil)
        C_pupil = [];
        p_pupil = [];
    else
        C_pupil = corr(dFF',pupil,'rows','complete');
        C_pupil_shuff = nan(num_reps,N);
        p_pupil = nan(N,1);
    end

    % Get shuffle distribution (null distribution)
    C_wsp_shuff = nan(num_reps,N);
    C_wamp_shuff = nan(num_reps,N);
    C_loco_shuff = nan(num_reps,N);
    C_speed_shuff = nan(num_reps,N);
    for rep = 1:num_reps
        T_shuffled = block_shuffle_time(T,acquisition_rate);
        dFF_shuff = dFF(:,T_shuffled);
        % Correlations after shuffling timing
        C_wsp_shuff(rep,:) = corr(dFF_shuff',whisk_set_point,'rows','complete');
        C_wamp_shuff(rep,:) = corr(dFF_shuff',whisk_amp,'rows','complete');
        C_loco_shuff(rep,:) = corr(dFF_shuff',loco,'rows','complete');
        C_speed_shuff(rep,:) = corr(dFF_shuff',speed,'rows','complete');
        if ~isempty(pupil)
            C_pupil_shuff(rep,:) = corr(dFF_shuff',pupil,'rows','complete');
        end
    end
    
    % Calculate two-tailed p value
    p_wsp = nan(N,1);
    p_wamp = nan(N,1);
    p_loco = nan(N,1);
    p_speed = nan(N,1);
    for k = 1:N
        p_wsp(k) = sum(C_wsp_shuff(:,k) > abs(C_wsp(k)) |  C_wsp_shuff(:,k) < -abs(C_wsp(k)))/num_reps;
        p_wamp(k) = sum(C_wamp_shuff(:,k) > abs(C_wamp(k)) |  C_wamp_shuff(:,k) < -abs(C_wamp(k)))/num_reps;
        p_loco(k) = sum(C_loco_shuff(:,k) > abs(C_loco(k)) |  C_loco_shuff(:,k) < -abs(C_loco(k)))/num_reps;
        p_speed(k) = sum(C_speed_shuff(:,k) > abs(C_speed(k)) |  C_speed_shuff(:,k) < -abs(C_speed(k)))/num_reps;
        if ~isempty(pupil)
            p_pupil(k) = sum(C_pupil_shuff(:,k) > abs(C_pupil(k)) |  C_pupil_shuff(:,k) < -abs(C_pupil(k)))/num_reps;
        end
    end
    
    
    
