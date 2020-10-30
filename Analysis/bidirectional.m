% This script plots the correlations of individual PFs with behaviour 

clear all; clc

define_dirs;

%% Histogram of ON/OFF GCs compared to shuffle test
% Figure 2B

change_dFF = cell(13,1);
p_val = cell(13,1);

for dataset_ix = 1:13
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,loco] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get A and QW states
    [A,QW] = define_behav_periods(whisk_amp,loco,acquisition_rate); 
    [change_dFF{dataset_ix},p_val{dataset_ix}] = change_dFF_sig(dFF,A,QW,acquisition_rate);
    
end

bins = linspace(-1.,1.,60);
bins_c = bins(2:end)-mean(diff(bins))/2;

change_dFF_all = vertcat(change_dFF{:});
p_all = vertcat(p_val{:});

C_pass_shuff = change_dFF_all(p_all < 0.05);
C_fail_shuff = change_dFF_all(p_all >= 0.05);

C_up = C_pass_shuff(C_pass_shuff>0);
C_down = C_pass_shuff(C_pass_shuff<0);

h_up = histcounts(C_up,bins);
h_down = histcounts(C_down,bins);
h_fail = histcounts(C_fail_shuff,bins);

figure, bar(bins_c,h_fail,'FaceColor',[.7,.7,.7],'EdgeColor','k','BarWidth',1)
hold on, bar(bins_c,h_up,'FaceColor',[1,.3,1],'EdgeColor','k','BarWidth',1)
hold on, bar(bins_c,h_down,'FaceColor',[.2,1,1],'EdgeColor','k','BarWidth',1)
set(gca,'FontSize',18)
set(gca,'Box','off')
xlabel('Change in dFF')
ylabel('Number')

fracs = [sum(h_up), sum(h_down),sum(h_fail)];
fracs = fracs/sum(fracs);
figure, p = pie(fracs);
set(gca,'FontSize',18)
t = p(1); t.FaceColor = [1,.3,1]; t.EdgeColor = [1,.3,1];
t = p(2); t.FontSize=20;
t = p(3); t.FaceColor = [.2,1,1]; t.EdgeColor = [.2,1,1];
t = p(4); t.FontSize=20;
t = p(5); t.FaceColor = [.7,.7,.7]; t.EdgeColor = [.7,.7,.7];
t = p(6); t.FontSize=20;
%% Plot distributions separately for all experiments
% Supplementary figure SX

for dataset_ix = 1:13
    C_pass_shuff = change_dFF{dataset_ix}(p_val{dataset_ix} < 0.05);
    C_fail_shuff = change_dFF{dataset_ix}(p_val{dataset_ix} >= 0.05);

    C_up = C_pass_shuff(C_pass_shuff>0);
    C_down = C_pass_shuff(C_pass_shuff<0);

    h_up = histcounts(C_up,bins);
    h_down = histcounts(C_down,bins);
    h_fail = histcounts(C_fail_shuff,bins);

    figure, bar(bins_c,h_fail,'FaceColor',[.7,.7,.7],'EdgeColor','k','BarWidth',1)
    hold on, bar(bins_c,h_up,'FaceColor','r','EdgeColor','k','BarWidth',1)
    hold on, bar(bins_c,h_down,'FaceColor','b','EdgeColor','k','BarWidth',1)
    set(gca,'FontSize',18)
    set(gca,'Box','off')
    xlabel('Change in dFF')
    ylabel('Number')
    
    title(datasets{dataset_ix},'Interpreter','None')

    fracs = [sum(h_up), sum(h_down),sum(h_fail)];
    fracs = fracs/sum(fracs);
    figure, p = pie(fracs);
    set(gca,'FontSize',18)
    t = p(1); t.FaceColor = 'r'; t.EdgeColor = 'r';
    t = p(2); t.FontSize=20;
    t = p(3); t.FaceColor = 'b'; t.EdgeColor = 'b';
    t = p(4); t.FontSize=20;
    t = p(5); t.FaceColor = [.7,.7,.7]; t.EdgeColor = [.7,.7,.7];
    t = p(6); t.FontSize=20;
    
    title(datasets{dataset_ix},'Interpreter','None')
    
end

%% Histogram of GC corrs with different behaviours
% Figure S3
% Warning - slow

C_wsp = cell(13,1);
p_wsp = cell(13,1);

C_wamp = cell(13,1);
p_wamp = cell(13,1);

C_loco = cell(13,1);
p_loco = cell(13,1);

C_speed = cell(13,1);
p_speed = cell(13,1);

C_pupil = cell(13,1);
p_pupil = cell(13,1);

for dataset_ix = 1:13
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,whisk_amp,loco,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get A and QW states
    [C_wsp{dataset_ix},p_wsp{dataset_ix},...
        C_wamp{dataset_ix},p_wamp{dataset_ix},...
        C_loco{dataset_ix},p_loco{dataset_ix},...
        C_speed{dataset_ix},p_speed{dataset_ix},...
        C_pupil{dataset_ix},p_pupil{dataset_ix}] ...
        = corr_sig(dFF,whisk_set_point,whisk_amp,loco,speed,pupil,acquisition_rate);
    
end
%% Generates figures for previous cell 
% Figure S3

% Replace below with whisker set point, amplitude, locomotion, or pupil
C = C_speed; p_val = p_speed;

bins = linspace(-1,1,50);
bins_c = bins(2:end)-mean(diff(bins))/2;

C_all = vertcat(C{:});
p_all = vertcat(p_val{:});

C_pass_shuff = C_all(p_all < 0.05);
C_fail_shuff = C_all(p_all >= 0.05);

C_up = C_pass_shuff(C_pass_shuff>0);
C_down = C_pass_shuff(C_pass_shuff<0);

h_up = histcounts(C_up,bins);
h_down = histcounts(C_down,bins);
h_fail = histcounts(C_fail_shuff,bins);

figure, bar(bins_c,h_fail,'FaceColor',[.7,.7,.7],'EdgeColor','k','BarWidth',1)
hold on, bar(bins_c,h_up,'FaceColor',[1,.3,1],'EdgeColor','k','BarWidth',1)
hold on, bar(bins_c,h_down,'FaceColor',[.2,1,1],'EdgeColor','k','BarWidth',1)
set(gca,'FontSize',18)
set(gca,'Box','off')
xlabel('Correlation')
ylabel('Number')

fracs = [sum(h_up), sum(h_down),sum(h_fail)];
fracs = fracs/sum(fracs);
figure, p = pie(fracs);
set(gca,'FontSize',18)
t = p(1); t.FaceColor = [1,.3,1]; t.EdgeColor = [1,.3,1];
t = p(2); t.FontSize=20;
t = p(3); t.FaceColor = [.2,1,1]; t.EdgeColor = [.2,1,1];
t = p(4); t.FontSize=20;
t = p(5); t.FaceColor = [.7,.7,.7]; t.EdgeColor = [.7,.7,.7];
t = p(6); t.FontSize=20;

%% Correlate behavioural data
% Figure S3

labels = {'wsp','wamp','loco','speed'};

C_wsp_wamp = nan(13,1);
C_wsp_loco = nan(13,1);
C_wsp_speed = nan(13,1);
C_wamp_loco = nan(13,1);
C_wamp_speed = nan(13,1);
C_loco_speed = nan(13,1);
for dataset_ix = 1:13
    [~,time,acsquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,whisk_amp,loco,speed] = load_behav_data(dataset_ix,time);
    C_wsp_wamp(dataset_ix) = corr(whisk_set_point,whisk_amp,'rows','complete');
    C_wsp_loco(dataset_ix) = corr(whisk_set_point,loco,'rows','complete');
    C_wsp_speed(dataset_ix) = corr(whisk_set_point,abs(speed),'rows','complete');
    C_wamp_loco(dataset_ix) = corr(whisk_amp,loco,'rows','complete');
    C_wamp_speed(dataset_ix) = corr(whisk_amp,abs(speed),'rows','complete');
    C_loco_speed(dataset_ix) = corr(loco,abs(speed),'rows','complete');
end

figure, plot(zeros,C_wsp_wamp,'ok')
hold on, plot(ones,C_wsp_loco,'ok')
plot(2*ones,C_wsp_speed,'ok')
plot(3*ones,C_wamp_loco,'ok')
plot(4*ones,C_wamp_speed,'ok')
plot(5*ones,C_loco_speed,'ok')
set(gca,'XTick',0:5);
set(gca,'XTickLabel',{'wsp-wamp','wsp-loco','wsp-speed','wamp-loco','wamp-speed','loco-speed'});
%% Plot examples of behavioural data
% Figure S3

dataset_ix = 13;
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,loco_smooth,speed_smooth] = load_behav_data(dataset_ix,time);

figure, plot(time,whisk_set_point,'k')
xlabel('time (s)'), ylabel('whisk set point (deg.)')

figure, plot(time,whisk_amp,'k')
xlabel('time (s)'), ylabel('whisk amplitude (deg.)')

figure, plot(time,loco_smooth,'k')
xlabel('time (s)'), ylabel('locomotion')

figure, plot(time,speed_smooth,'k')
xlabel('time (s)'), ylabel('speed')

%% Show example of population activity
% Requires running previous cell 
% Figure 1C

dataset_ix = 6;
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[N,~] = size(dFF);

ix_fail = find(p_val{dataset_ix} >.05);

ix_up = find(p_val{dataset_ix} < .05 & change_dFF{dataset_ix} > 0);
C = change_dFF{dataset_ix}(ix_up);
[~,ix] = sort(C,'descend');
C = C(ix); ix_up = ix_up(ix);

ix_down = find(p_val{dataset_ix} < .05 & change_dFF{dataset_ix} < 0);
C = change_dFF{dataset_ix}(ix_down);
[~,ix] = sort(-C,'descend');
C = C(ix); ix_down = ix_down(ix);

% CAREFUL : INVERTED COLOR SCALE HERE
figure, imagesc(time,1:N,-dFF([ix_up; ix_down; ix_fail],:))
 colormap(gray), caxis([-1.5,0])
set(gca,'FontSize',15,'YTick',[1,700])
xlabel('Time (s)'), ylabel('Axon number')

[~,whisk_set_point,whisk_amp,loco,speed] = load_behav_data(dataset_ix,time);

figure, plot(time, whisk_set_point,'k','LineWidth',1.5)
set(gca,'FontSize',15,'Box','off')
xlim([time(1),time(end)])
xlabel('Time (s)'), ylabel('WSP (deg.)')

figure, plot(time, speed,'k','LineWidth',1.5)
set(gca,'FontSize',15,'Box','off')
xlim([time(1),time(end)])
xlabel('Time (s)'), ylabel('Speed')

%% Same example as previous cell, plot example of A / QW periods
% All figure parameters optimized for given example dataset
% Figure 2A

[A,QW] = define_behav_periods(whisk_amp,loco,acquisition_rate); 

% Plot whisker set point
figure, hold on
plot(time, whisk_set_point,'k','LineWidth',1.5)
for k = 1:length(QW)
    ix = QW(k,1):QW(k,2);
    plot(time(ix), whisk_set_point(ix),'c','LineWidth',1.5)
end

for k = 1:length(A)
    ix = A(k,1):A(k,2);
    plot(time(ix), whisk_set_point(ix),'m','LineWidth',1.5)
end
set(gca,'FontSize',15,'Box','off')
axis([220,280,-20,60])
xlabel('Time (s)'), ylabel('WSP (deg.)')

% Plot speed
figure, hold on
plot(time, speed,'k','LineWidth',1.5)
for k = 1:length(QW)
    ix = QW(k,1):QW(k,2);
    plot(time(ix), speed(ix),'c','LineWidth',1.5)
end

for k = 1:length(A)
    ix = A(k,1):A(k,2);
    plot(time(ix), speed(ix),'m','LineWidth',1.5)
end
set(gca,'FontSize',15,'Box','off')
axis([220,280,-5,20])
xlabel('Time (s)'), ylabel('Speed')

% Plot example dFF
figure, hold on
for k = 1:length(QW)
    rectangle('Position',[time(QW(k,1)) -20 time(QW(k,2))-time(QW(k,1)) 80],'FaceColor',[.75 1 1],'EdgeColor','w')
end
for k = 1:length(A)
    rectangle('Position',[time(A(k,1)) -20 time(A(k,2))-time(A(k,1)) 80],'FaceColor',[1 .9 1],'EdgeColor','w')
end
ix_display_on = [42,60,115,142,536];
ix_display_off = [9,149,226];
for k = 1:5
    plot(time,dFF(ix_display_on(k),:)+k*3,'Color',[.3,.3,.3],'LineWidth',1)
end
for k=1:3
    plot(time,dFF(ix_display_off(k),:)-k*3,'Color',[.3,.3,.3],'LineWidth',1)
end
set(gca,'FontSize',15,'Box','off','YTick',[0,2])
axis([220,280,-10,20])
xlabel('Time (s)'), ylabel('\Delta F/F')

%% Show full example of up / down /non mod 
% Figure S5

figure, hold on
for k = 1:length(QW)
    rectangle('Position',[time(QW(k,1)) -20 time(QW(k,2))-time(QW(k,1)) 80],'FaceColor',[.75 1 1],'EdgeColor','w')
end
for k = 1:length(A)
    rectangle('Position',[time(A(k,1)) -20 time(A(k,2))-time(A(k,1)) 80],'FaceColor',[1 .9 1],'EdgeColor','w')
end

ix_display_fail = [4,13,65];
ix_display_on = [42,60,115,142,536];
ix_display_off = [9,149,226]; % dataset 6

for k = 1:5
    plot(time,dFF(ix_display_on(k),:)+k*3,'Color',[.3,.3,.3],'LineWidth',1)
end
for k=1:3
    plot(time,dFF(ix_display_off(k),:)-k*3-2,'Color',[.3,.3,.3],'LineWidth',1)
end
for k = 1:3
    plot(time,dFF(ix_display_fail(k),:)+k*3+18,'Color',[.3,.3,.3],'LineWidth',1)
end
set(gca,'FontSize',15,'Box','off','YTick',[0,2])
axis([0,400,-13,31])
xlabel('Time (s)'), ylabel('\Delta F/F')
%% Save avg of negatively modulated PFs
% Supplementary figure SX

figure,
for dataset_ix = 1:13
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    
    ix_NM = find(change_dFF{dataset_ix} < 0 & p_val{dataset_ix} < 0.05);
    
    dFF_NM = mean(dFF(ix_NM,:));
    
    plot(time,dFF_NM,'k','LineWidth',1.2)
    
    save(['../Data/',datasets{dataset_ix},'_avg_NM'],'time','dFF_NM')
    clear time dFF_NM dFF
end


%% Get distribution of SNRs of all axons comparing non modulated to others
% Figure S5

clear all
clc

define_dirs;

addpath('../ROI Grouping/Utilities/')
addpath('../ROI Grouping/From CNMF_E/')

SNR_up = [];
SNR_down = [];
SNR_fail = [];

for dataset_ix = 1:13
    fname = datasets{dataset_ix};
    disp(fname)

    [~,time,acquisition_rate] = load_data(dataset_ix);
    load([basedir,fname,'/',fname,'.mat'],'Numb_patches')
    load([basedir,fname,'/processed/',fname,'_GroupedData.mat'],'Ain_axons','dFF_axons')

    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate); 
    
    
    for patch_no = 1:Numb_patches
        patch_no
        
        [change_dFF,p_val] = change_dFF_sig(dFF_axons{patch_no},A,QW,acquisition_rate);
        
        ix_up = find(p_val < 0.05 & change_dFF>0);
        ix_down = find(p_val < 0.05 & change_dFF<0);
        ix_fail = find(p_val > 0.05);
        
        
        % Load data
        load([basedir,fname,'/raw/Patch',sprintf('%03d',patch_no),'.mat'])
        Y = double(Y);
        
        % Get SNR for all ROIs
        [~,SNR_all,~,~] = remove_bad_cells(Ain_axons{patch_no},Y,[d1,d2],acquisition_rate,0);
        

        SNR_up = [SNR_up; SNR_all(ix_up)];
        SNR_down = [SNR_down; SNR_all(ix_down)];
        SNR_fail = [SNR_fail; SNR_all(ix_fail)];
    end
end

figure, histogram(SNR_up,1:35,'FaceColor','r')
set(gca,'FontSize',15)
xlabel('SNR')
ylabel('Number')

figure, histogram(SNR_down,1:35,'FaceColor','b')
set(gca,'FontSize',15)
xlabel('SNR')
ylabel('Number')

figure, histogram(SNR_fail,1:35,'FaceColor','k')
set(gca,'FontSize',15)
xlabel('SNR')
ylabel('Number')




