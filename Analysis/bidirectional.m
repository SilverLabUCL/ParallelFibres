% This script plots the correlations of individual PFs with behaviour 

clear all; clc

define_dirs;

%% Histogram of ON/OFF GCs compared to shuffle test
% Figure 1E

change_dFF = cell(15,1);
p_val = cell(15,1);

for dataset_ix = 1:15
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,~,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get A and QW states
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate); 
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

%% Histogram of GC corrs with different behaviours
% Figure S3
% Warning - slow

C_wsp = cell(15,1);
p_wsp = cell(15,1);

C_wamp = cell(15,1);
p_wamp = cell(15,1);

C_speed = cell(15,1);
p_speed = cell(15,1);

C_pupil = cell(15,1);
p_pupil = cell(15,1);

for dataset_ix = 1:15
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get A and QW states
    [C_wsp{dataset_ix},p_wsp{dataset_ix},...
        C_wamp{dataset_ix},p_wamp{dataset_ix},...
        C_speed{dataset_ix},p_speed{dataset_ix},...
        C_pupil{dataset_ix},p_pupil{dataset_ix}] ...
        = corr_sig(dFF,whisk_set_point,whisk_amp,speed,pupil,acquisition_rate);
    
end
%% Generates figures for previous cell

% Replace below with whisker set point, amplitude, locomotion, or pupil
C = C_wsp; p_val = p_wsp;

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

[~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate); 

figure, plot(time, whisk_set_point,'k','LineWidth',1.5)
set(gca,'FontSize',15,'Box','off')
xlim([time(1),time(end)])
xlabel('Time (s)'), ylabel('WSP (deg.)')

figure, plot(time, speed,'k','LineWidth',1.5)
set(gca,'FontSize',15,'Box','off')
xlim([time(1),time(end)])
xlabel('Time (s)'), ylabel('Loc. (std)')

%% Same example as previous cell, plot example of A / QW periods
% All figure parameters optimized for given example dataset
% Figure 1D

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
axis([220,280,-2,12])
xlabel('Time (s)'), ylabel('Loc. (std)')

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
ix_display_off = [9,149,226];
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

%% Get distribution of SNRs of all axons comparing non modulated to others
% Figure S5
addpath('../ROI Grouping/Utilities/')
addpath('../ROI Grouping/From CNMF_E/')

SNR_up = [];
SNR_down = [];
SNR_fail = [];

for dataset_ix = 1:15
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
%% Show example of single PF tuning curves
% Figure 1C

dataset_ix = 6;
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[N,~] = size(dFF);
[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,1);

[change_dFF,p_val] = change_dFF_sig(dFF,A,QW,acquisition_rate);
ix_fail = find(p_val >.05);
ix_ON = find(p_val < .05 & change_dFF > 0);
ix_OFF = find(p_val < .05 & change_dFF < 0);

dwsp = 10;
wsp_bin = -20:dwsp:60;
tuningcurves = nan(N,size(wsp_bin,2));
tuningcurves_ste = nan(N,size(wsp_bin,2));
for k = 1:size(wsp_bin,2)
    ix = find(whisk_set_point > wsp_bin(k) - dwsp/2 & whisk_set_point < wsp_bin(k) + dwsp/2);
    tuningcurves(:,k) = mean(dFF(:,ix),2);
    tuningcurves_ste(:,k) = std(dFF(:,ix),[],2) / sqrt(numel(ix));
end

%%
n = 168;
figure, plot(wsp_bin,tuningcurves(n,:),'-sk')
 hold on,
for k = 1:size(wsp_bin,2)
    plot(wsp_bin(k)*[1,1],tuningcurves(n,k)+[-1,1]*tuningcurves_ste(n,k),'k')
end
