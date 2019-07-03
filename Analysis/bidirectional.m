% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;
        
[C_spd, C_ang, C_wsp, C_amp] = deal(cell(8,1)); 
[p_spd, p_ang, p_wsp, p_amp] = deal(cell(8,1)); 

tic
for dataset_ix =  1:8
    
    % Load data
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    N = size(dFF,1);
    
    % Load behavioural data
    [whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    
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
        
    end
    toc
    
end

save([basedir,'processed/corr_axons_behav.mat'],'C_spd','p_spd','C_ang','p_ang','C_wsp','p_wsp','C_amp','p_amp')

%% Plot histograms of correlations for all neurons

bins = linspace(-1,1,60);
bins_c = bins(2:end)-mean(diff(bins))/2;

% Plot whisker set point
for k = 3%1:4
    switch k
        case 1
            C_temp = C_spd; p_temp = p_spd; title_name = 'speed';
        case 2
            C_temp = C_ang; p_temp = p_ang; title_name = 'whisker angle';
        case 3
            C_temp = C_wsp; p_temp = p_wsp; title_name = 'whisker set point';
        case 4
            C_temp = C_amp; p_temp = p_amp; title_name = 'whisker amplitude';
    end
    
    C_all = vertcat(C_temp{:});
    p_all = vertcat(p_temp{:});

    C_pass_shuff = C_all(p_all < 0.05);
    C_fail_shuff = C_all(p_all >= 0.05);

    h_pass = histcounts(C_pass_shuff,bins);
    h_fail = histcounts(C_fail_shuff,bins);

    figure, bar(bins_c,h_fail,'FaceColor',[.6,.8,1],'EdgeColor','w')
    hold on, bar(bins_c,h_pass,'FaceColor','b','EdgeColor','w')
    set(gca,'FontSize',18)
    title(title_name)
end

%% Correlation - different time lags

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
    
    [C_spd_lag{dataset_ix}, C_ang_lag{dataset_ix}, C_wsp_lag{dataset_ix}, C_amp_lag{dataset_ix}...
        C_spd_lag_std{dataset_ix}, C_ang_lag_std{dataset_ix}, C_wsp_lag_std{dataset_ix}, C_amp_lag_std{dataset_ix}]...
        = deal(nan(N,length(lags_temp))); 
    
    % Get SNR for all ROIS
    SNR{dataset_ix} = nan(N,1);
    for n = 1:N
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

%% SNR distribution is the same

load([basedir,'processed/corr_axons_behav.mat'])
bins = linspace(0,150,60);
bins_c = bins(2:end)-mean(diff(bins))/2;

% Plot whisker set point
SNR_all = vertcat(SNR{:});
p_all = vertcat(p_wsp{:});

SNR_pass_shuff = SNR_all(p_all < 0.05);
SNR_fail_shuff = SNR_all(p_all >= 0.05);

h_pass = histcounts(SNR_pass_shuff,bins);
h_fail = histcounts(SNR_fail_shuff,bins);

figure, bar(bins_c,h_fail,'FaceColor',[.6,.8,1],'EdgeColor','w')
figure, bar(bins_c,h_pass,'FaceColor','b','EdgeColor','w')
set(gca,'FontSize',18)

%% Get histogram of lags at maximum correlations

[lags_up, C_max_up, lags_down, C_max_down] = deal([]);
for dataset_ix =  1:8
    for n = 1:size(C_wsp_lag{dataset_ix},1)
        % Only consider significant
        if p_wsp{dataset_ix}(n) < .05
            [amax,bmax] = max(C_wsp_lag{dataset_ix}(n,:));
            [amin,bmin] = min(C_wsp_lag{dataset_ix}(n,:));
            if abs(lags{dataset_ix}(bmax)) < abs(lags{dataset_ix}(bmin))
                lags_up = [lags_up; lags{dataset_ix}(bmax)];
                C_max_up = [C_max_up; amax];
            else
                lags_down = [lags_down; lags{dataset_ix}(bmin)];
                C_max_down = [C_max_down; amin];
            end
        end
    end
end

figure, histogram([lags_up;lags_down]*1000,-200:15:200,'Normalization','probability')
