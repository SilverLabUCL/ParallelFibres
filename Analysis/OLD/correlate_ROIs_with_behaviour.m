% This script calculates the correlations of individual ROIs with behaviour
% using a shuffle test and saves the results

% This section calculates just 0-lag corrs, but with shuffle test

% C - correlation
% p - pvalue

clear all; clc

define_dirs;
        
[C_spd, C_ang, C_wsp, C_amp] = deal(cell(8,1)); 
[p_spd, p_ang, p_wsp, p_amp] = deal(cell(8,1)); 

tic
for dataset_ix =  [1:6,8,16,17]
    
    % Load data
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    N = size(dFF,1);
    
    % Load behavioural data
    [whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    notR = define_not_running(speed,acquisition_rate,1); 
    
    [C_spd{dataset_ix}, C_ang{dataset_ix}, C_wsp{dataset_ix}, C_amp{dataset_ix}] = deal(nan(N,1)); 
    [p_spd{dataset_ix}, p_ang{dataset_ix}, p_wsp{dataset_ix}, p_amp{dataset_ix}] = deal(nan(N,1)); 

    for n = 1:N
        
        [C,p] = corr_sig(dFF(n,:),[speed,whisk_angle,whisk_set_point,whisk_amp],acquisition_rate);
        
        C_spd{dataset_ix}(n) = C(1);
        C_ang{dataset_ix}(n) = C(2);
        C_wsp{dataset_ix}(n) = C(3);
        C_amp{dataset_ix}(n) = C(4);

        p_spd{dataset_ix}(n) = p(1);
        p_ang{dataset_ix}(n) = p(2);
        p_wsp{dataset_ix}(n) = p(3);
        p_amp{dataset_ix}(n) = p(4);
        
        [C,p] = corr_sig(dFF(n,notR),[speed(notR),whisk_angle(notR),whisk_set_point(notR),whisk_amp(notR)],acquisition_rate);
        
        C_spd_notR{dataset_ix}(n) = C(1);
        C_ang_notR{dataset_ix}(n) = C(2);
        C_wsp_notR{dataset_ix}(n) = C(3);
        C_amp_notR{dataset_ix}(n) = C(4);

        p_spd_notR{dataset_ix}(n) = p(1);
        p_ang_notR{dataset_ix}(n) = p(2);
        p_wsp_notR{dataset_ix}(n) = p(3);
        p_amp_notR{dataset_ix}(n) = p(4);
        
    end
    toc
    
end

save([basedir,'processed/corr_axons_behav.mat'],'C_spd','p_spd',...
    'C_ang','C_amp','p_ang','C_wsp','p_wsp','C_amp','p_amp','p_amp',...
    'C_spd_notR','p_spd_notR','C_ang_notR','C_amp_notR','p_ang_notR',...
    'C_wsp_notR','p_wsp_notR','C_amp_notR','p_amp_notR','p_amp_notR');

%% This section calculates lagged correlations, but without shuffle test
% Also calculates SNR 

clear all; clc

define_dirs;

[SNR, C_spd_lag, C_ang_lag, C_wsp_lag, C_amp_lag, lags] = deal(cell(8,1));
[C_spd_lag_std, C_ang_lag_std, C_wsp_lag_std, C_amp_lag_std] = deal(cell(8,1));

for dataset_ix =  1:8
    
    % Load data
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    N = size(dFF,1);
    
    % Load behavioural data
    [whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    
    [~,lags_temp] = corr_lag_trial(dFF(1,:),whisk_angle,time);
    
    lags{dataset_ix} = lags_temp / acquisition_rate;
    
    [C_spd_lag{dataset_ix}, C_ang_lag{dataset_ix}, C_wsp_lag{dataset_ix}, C_amp_lag{dataset_ix},...
        C_spd_lag_std{dataset_ix}, C_ang_lag_std{dataset_ix}, C_wsp_lag_std{dataset_ix}, C_amp_lag_std{dataset_ix}]...
        = deal(nan(N,length(lags_temp))); 
    
    % Get SNR and lagged cross-correlogram (with behaviour) for all ROIS
    SNR{dataset_ix} = nan(N,1);
    for n = 1:N
        % Calculate SNR
        SNR{dataset_ix}(n) = max(medfilt1(dFF(n,:),round(.2*acquisition_rate)))/GetSn(dFF(n,:));
        
        C_temp = corr_lag_trial(dFF(n,:),whisk_angle,time);
        C_ang_lag{dataset_ix}(n,:) = mean(C_temp,2);
        C_ang_lag_std{dataset_ix}(n,:) = std(C_temp,[],2);
        
        C_temp = corr_lag_trial(dFF(n,:),whisk_set_point,time);
        C_wsp_lag{dataset_ix}(n,:) = mean(C_temp,2);
        C_wsp_lag_std{dataset_ix}(n,:) = std(C_temp,[],2);
        
        C_temp = corr_lag_trial(dFF(n,:),whisk_amp,time);
        C_amp_lag{dataset_ix}(n,:) = mean(C_temp,2);
        C_amp_lag_std{dataset_ix}(n,:) = std(C_temp,[],2);
        
        C_temp = corr_lag_trial(dFF(n,:),speed,time);
        C_spd_lag{dataset_ix}(n,:) = mean(C_temp,2);
        C_spd_lag_std{dataset_ix}(n,:) = std(C_temp,[],2);
    end
end

save([basedir,'processed/corr_axons_behav_lags.mat'], 'SNR', 'C_spd_lag', 'C_ang_lag', 'C_wsp_lag', 'C_amp_lag', 'lags', 'C_spd_lag_std', 'C_ang_lag_std', 'C_wsp_lag_std', 'C_amp_lag_std');
