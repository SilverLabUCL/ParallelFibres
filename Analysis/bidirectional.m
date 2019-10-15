% This script plots the correlations of individual PFs with behaviour 

clear all; clc

define_dirs;

% Load pre-calculated correlations
load([basedir,'processed/corr_axons_behav.mat'])
load([basedir,'processed/corr_axons_behav_lags.mat'])

%% Histogram of correlations of each PF with behaviours 
% Breaks down into up, down and non-modulated PFs
% Note that correlations were calculated previously in script:
%       correlate_ROIs_with_behaviour.m

bins = linspace(-1,1,60);
bins_c = bins(2:end)-mean(diff(bins))/2;
 
C_all = vertcat(C_wsp{:});
p_all = vertcat(p_wsp{:});

C_pass_shuff = C_all(p_all < 0.05);
C_fail_shuff = C_all(p_all >= 0.05);

C_up = C_pass_shuff(C_pass_shuff>0);
C_down = C_pass_shuff(C_pass_shuff<0);

h_up = histcounts(C_up,bins);
h_down = histcounts(C_down,bins);
h_fail = histcounts(C_fail_shuff,bins);

figure, bar(bins_c,h_fail,'FaceColor',[.74,.74,.74],'EdgeColor',[.74,.74,.74])
hold on, bar(bins_c,h_up,'FaceColor',[1,.84,0],'EdgeColor',[1,.84,0])
hold on, bar(bins_c,h_down,'FaceColor',[.1,.4,1],'EdgeColor',[.1,.4,1])
set(gca,'FontSize',18)
title('Whisker set point')
set(gca,'Box','off')
xlabel('Correlation')
ylabel('Number')

fracs = [sum(h_up), sum(h_down),sum(h_fail)];
fracs = fracs/sum(fracs);
figure, p = pie(fracs,{'up','down','n.s.'})
set(gca,'FontSize',18)
title('Whisker set point')    
t = p(1); t.FaceColor = [1,.84,0]; t.EdgeColor = [1,.84,0];
t = p(2); t.FontSize=20;
t = p(3); t.FaceColor = [.1,.4,1]; t.EdgeColor = [.1,.4,1];
t = p(4); t.FontSize=20;
t = p(5); t.FaceColor = [.74,.74,.74]; t.EdgeColor = [.74,.74,.74];
t = p(6); t.FontSize=20;


%% SNR Distribution of SNR is similar for non-modulated as ON/OFF

load([basedir,'processed/corr_axons_behav.mat'])
bins = linspace(0,150,30);
bins_c = bins(2:end)-mean(diff(bins))/2;

SNR_all = vertcat(SNR{:});
p_all = vertcat(p_wsp{:});
C_all = vertcat(C_wsp{:});

SNR_up = SNR_all(p_all < 0.05 & C_all > 0);
SNR_down = SNR_all(p_all < 0.05 & C_all < 0);
SNR_fail_shuff = SNR_all(p_all >= 0.05);

figure, histogram(SNR_fail_shuff,bins,'Normalization','probability','FaceAlpha',1,'FaceColor',[.74,.74,.74],'EdgeColor',[.74,.74,.74]);
hold on, histogram(SNR_up,bins-10,'Normalization','probability','FaceAlpha',0,'EdgeColor',[1,.84,0],'LineWidth',2);
hold on, histogram(SNR_down,bins+10,'Normalization','probability','FaceAlpha',0,'EdgeColor',[.1,.4,1],'LineWidth',2);

set(gca,'FontSize',18)
set(gca,'Box','off')
xlim([0,150])
xlabel('SNR')
ylabel('Fraction')

%% Plot examples of high SNR non-modulated PFs

load([basedir,'processed/corr_axons_behav.mat'])

dataset_ix = 3;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[~,whisk_set_point,~,speed] = load_behav_data(dataset_ix,time);

ix = find(p_wsp{dataset_ix} < .05 & abs(C_wsp{dataset_ix}) < .2);
SNR_nonmod_lowC = SNR{dataset_ix}(ix);
[SNR_nonmod_lowC,resort_ix] = sort(SNR_nonmod_lowC);
ix = ix(resort_ix);

figure, hold on
for k = 1:5
    plot(time,dFF(ix(k),:)+1.1*k,'k','LineWidth',1), [C_wsp{dataset_ix}(ix(k)),SNR{dataset_ix}(ix(k))]
end
xlabel('Time (s)'), ylabel('DF/F')
set(gca,'FontSize',18)
xlim([100,200]), ylim([.5,7])

figure, plot(time,rad2deg(whisk_set_point),'Color',[.1,.4,1],'LineWidth',1)
ylabel('Whisker set point (deg)'), xlabel('Time (s)')
set(gca,'FontSize',18), set(gca,'Box','off')
xlim([100,200]), ylim([-40,20])

figure, plot(time,speed/std(speed),'Color',[1,.3,0],'LineWidth',1)
ylabel('Locomotion (std)'), xlabel('Time (s)')
set(gca,'FontSize',18), set(gca,'Box','off')
xlim([100,200]), ylim([0,15])

%% Get histogram of lags at maximum correlations

% Find correlations
[lags_up, lags_down, C_max_lag] = deal(cell(8,1));
for dataset_ix =  1:8
    for n = 1:size(C_wsp_lag{dataset_ix},1)
        % Only consider axons with significant correlation at 0-lag
        if p_wsp{dataset_ix}(n) < .05
            [amax,bmax] = max(C_wsp_lag{dataset_ix}(n,:));
            [amin,bmin] = min(C_wsp_lag{dataset_ix}(n,:));
            if abs(lags{dataset_ix}(bmax)) < abs(lags{dataset_ix}(bmin))
                lags_up{dataset_ix} = [lags_up{dataset_ix}; lags{dataset_ix}(bmax)];
                C_max_lag{dataset_ix} = [C_max_lag{dataset_ix}; amax];
            else
                lags_down{dataset_ix} = [lags_down{dataset_ix}; lags{dataset_ix}(bmin)];
                C_max_lag{dataset_ix} = [C_max_lag{dataset_ix}; amin];
            end
        else
            C_max_lag{dataset_ix} = [C_max_lag{dataset_ix}; nan];
        end
    end
end


bins = -200:15:200;

figure, histogram(vertcat(lags_up{:})*1000,bins-2,'Normalization','probability','FaceAlpha',0,'EdgeColor',[1,.84,0],'LineWidth',2);
hold on, histogram(vertcat(lags_down{:})*1000,bins+2,'Normalization','probability','FaceAlpha',0,'EdgeColor',[.1,.4,1],'LineWidth',2);
xlabel('Time lag (ms)')
ylabel('Fraction')

set(gca,'FontSize',15)
set(gca,'Box','off')
xlim([-200,200])

kstest2(vertcat(lags_up{:}),vertcat(lags_down{:}))

%% Show example of diversity of responses

dataset_ix = 2;
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[N,T] = size(dFF);
    
ix_fail = find(p_wsp{dataset_ix} >.05);
ix_up = find(p_wsp{dataset_ix} < .05 & C_max_lag{dataset_ix} > 0);
ix_down = find(p_wsp{dataset_ix} < .05 & C_max_lag{dataset_ix} < 0);

% CAREFUL : INVERTED COLOR SCALE HERE
figure, imagesc((1:T)/acquisition_rate,1:N,-dFF([ix_up; ix_down; ix_fail],:))
colormap(gray), caxis([-1,0])
set(gca,'FontSize',15)

[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
figure, plot(time, whisk_set_point,'k','LineWidth',1.5)
set(gca,'FontSize',15)

[~, I] = sort(lags_up{dataset_ix});
ix_up = ix_up(I);
[~, I] = sort(lags_down{dataset_ix});
ix_down = ix_down(I);

%% Plot avg over onset times

%t_onset = [220,325,471,588,749,1129,1399,1552,1706,1828,1922,2204,2485,2693,2648,3142,3262,3361,3468,3673,3808,3945,4025,4210,4362,4483,4647,4866,5021,5176,5290]; % d2

[~, t_onset] = get_onsets(whisk_amp,acquisition_rate);

%t_onset = t_onset(5:8);

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
figure, imagesc(t_,1:N,zscore(dFF_avg_onset([ix_up;ix_down;ix_fail],:)')'), caxis([-5,5]), colormap(bluewhitered)%caxis([0,2]), colormap(hot)
figure, plot_error_snake(t_,WSP_avg,[1,0,0])
figure, plot_error_snake(t_,(speed_avg),[1,0,0])

%%

    