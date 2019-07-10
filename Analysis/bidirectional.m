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
for k = 1%1:4
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
    
    C_up = C_pass_shuff(C_pass_shuff>0);
    C_down = C_pass_shuff(C_pass_shuff<0);

    h_up = histcounts(C_up,bins);
    h_down = histcounts(C_down,bins);
    h_fail = histcounts(C_fail_shuff,bins);

    figure, bar(bins_c,h_fail,'FaceColor',[.4,.8,1],'EdgeColor',[.4,.8,1])
    hold on, bar(bins_c,h_up,'FaceColor',[.1,.4,1],'EdgeColor',[.1,.4,1])
    hold on, bar(bins_c,h_down,'FaceColor',[.1,.4,1],'EdgeColor',[.1,.4,1])
    set(gca,'FontSize',18)
    title(title_name)
    set(gca,'Box','off')
    
    fracs = [sum(h_up), sum(h_down),sum(h_fail)];
    fracs = fracs/sum(fracs);
    figure, pie(fracs,{'','',''})
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

save([basedir,'processed/corr_axons_behav_lags.mat'], 'SNR', 'C_spd_lag', 'C_ang_lag', 'C_wsp_lag', 'C_amp_lag', 'lags', 'C_spd_lag_std', 'C_ang_lag_std', 'C_wsp_lag_std', 'C_amp_lag_std');


%% SNR distribution is the same

load([basedir,'processed/corr_axons_behav.mat'])
bins = linspace(0,150,30);
bins_c = bins(2:end)-mean(diff(bins))/2;

% Plot whisker set point
SNR_all = vertcat(SNR{:});
p_all = vertcat(p_wsp{:});
C_all = vertcat(C_wsp{:});

SNR_up = SNR_all(p_all < 0.05 & C_all > 0);
SNR_down = SNR_all(p_all < 0.05 & C_all < 0);
SNR_fail_shuff = SNR_all(p_all >= 0.05);

figure, histogram(SNR_fail_shuff,bins,'Normalization','pdf','FaceAlpha',1,'FaceColor',[.74,.74,.74],'EdgeColor',[.74,.74,.74]);
hold on, histogram(SNR_up,bins-10,'Normalization','pdf','FaceAlpha',0,'EdgeColor',[1,.84,0],'LineWidth',2);
hold on, histogram(SNR_down,bins+10,'Normalization','pdf','FaceAlpha',0,'EdgeColor',[.1,.4,1],'LineWidth',2);

set(gca,'FontSize',15)
set(gca,'Box','off')
xlim([0,150])
ylim([0,.05])
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


bins = -200:15:200;

h_up = histcounts(lags_up*1000,bins,'Normalization','pdf');
h_down = histcounts(lags_down*1000,bins,'Normalization','pdf');
figure,
hold on, histogram(lags_up*1000,bins-2,'Normalization','pdf','FaceAlpha',0,'EdgeColor',[1,.84,0],'LineWidth',2);
hold on, histogram(lags_down*1000,bins+2,'Normalization','pdf','FaceAlpha',0,'EdgeColor',[.1,.4,1],'LineWidth',2);

set(gca,'FontSize',15)
set(gca,'Box','off')
xlim([-200,200])

kstest2(lags_up,lags_down)

%% 

define_dirs;
dataset_ix = 2;
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[N,T] = size(dFF);
    
load([basedir,'processed/corr_axons_behav.mat'])
ix_fail = find(p_wsp{dataset_ix} >.05);
ix_up = find(p_wsp{dataset_ix} < .05 & C_wsp{dataset_ix} > 0);
ix_down = find(p_wsp{dataset_ix} < .05 & C_wsp{dataset_ix} < 0);

figure, imagesc((1:T)/acquisition_rate,1:N,dFF)
colormap(hot), caxis([0,3])
figure, imagesc((1:T)/acquisition_rate,1:N,dFF([ix_up; ix_down; ix_fail],:))
colormap(hot), caxis([0,3])
set(gca,'FontSize',15)

[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
figure, plot(time, whisk_set_point,'k','LineWidth',1.5)
set(gca,'FontSize',15)
%%

define_dirs;

load([basedir,'processed/corr_axons_behav_lags.mat'])
load([basedir,'processed/corr_axons_behav_lags.mat'])

[lags_up, C_max_up, lags_down, C_max_down, ix_up, ix_down] = deal([]);
for n = 1:size(C_wsp_lag{dataset_ix},1)
    % Only consider significant
    if p_wsp{dataset_ix}(n) < .05
        [amax,bmax] = max(C_wsp_lag{dataset_ix}(n,:));
        [amin,bmin] = min(C_wsp_lag{dataset_ix}(n,:));
        if abs(lags{dataset_ix}(bmax)) < abs(lags{dataset_ix}(bmin))
            lags_up = [lags_up; lags{dataset_ix}(bmax)];
            C_max_up = [C_max_up; amax];
            ix_up = [ix_up; n];
        else
            lags_down = [lags_down; lags{dataset_ix}(bmin)];
            C_max_down = [C_max_down; amin];
            ix_down = [ix_down; n];
        end
    end
end
[lags_up, I] = sort(lags_up);
ix_up = ix_up(I);
[lags_down, I] = sort(lags_down);
ix_down = ix_down(I);


%% Plot avg over onset times

%t_onset = [149,315,407,475,544,702,784,995,1097,1366,1534,1614,1672,1714,1906,2029,2248,2337,2473,2596,2712,2898,3021,3087,3155,3206,3272,3400,3557,3717,3978,4079,4186,4305,4407,4530,4679,4799,4969,5218,5444]; %d1
t_onset = [220,325,471,588,749,1129,1399,1552,1706,1828,1922,2204,2485,2693,2648,3142,3262,3361,3468,3673,3808,3945,4025,4210,4362,4483,4647,4866,5021,5176,5290]; % d2

t_before = round(2 * acquisition_rate);
t_after = round(2 * acquisition_rate);
dFF_avg_onset = zeros(N,t_before+t_after);
WSP_avg =  zeros(length(t_onset),t_before+t_after);
speed_avg =  zeros(length(t_onset),t_before+t_after);
for ix = 1:length(t_onset)
    t_ = (t_onset(ix)-t_after):(t_onset(ix)+t_after-1);
    dFF_avg_onset = dFF_avg_onset + dFF(:,t_);
    WSP_avg(ix,:) = whisk_set_point(t_);
    speed_avg(ix,:) = speed(t_);
end

dFF_avg_onset = dFF_avg_onset/length(t_onset);


t_ = (((1:size(WSP_avg,2))-t_before)/acquisition_rate);
figure, imagesc(t_,1:length(ix_up),zscore(dFF_avg_onset(ix_up,:)')'), caxis([-5,5]), colormap(bluewhitered)%caxis([0,2]), colormap(hot)
figure, imagesc(t_,1:length(ix_down),zscore(dFF_avg_onset(ix_down,:)')'), caxis([-5,]), colormap(bluewhitered)
figure, plot_error_snake(t_,WSP_avg,[1,0,0])
figure, plot_error_snake(t_,(speed_avg),[1,0,0])

%%
for
    