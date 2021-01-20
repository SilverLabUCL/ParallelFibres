% Bidirectional for somatic data 

clear all; clc

basedir = '~/Documents/ParallelFibres/Data_Soma/';
addpath('Utilities/')

% Main datasets
datasets = {'FL65_170628_12_13_56_soma',... 1
            'FL_S_170905_11_40_44_soma',... 2
            'FL75_171218_11_42_30_soma',...	3
            'FL_S_171109_14_41_35_soma',... 4
            'FL75_180112_10_45_22_soma',...	5
            'FL92_180228_10_26_01_soma'};  % 6

         
%% Histogram of ON/OFF GCs compared to shuffle test
% Figure 2B

change_dFF = cell(6,1);
p_val = cell(6,1);

for dataset_ix = 1:6
    fname = datasets{dataset_ix};
    load([basedir,fname,'.mat'])
    
    time = mean(timeseries,2)/1000;
    
    % Whisking variables
    [whisk_time,whisk_angle,whisk_set_point,whisk_amp,whisk_phase] ...
        = convert_dlc_to_whisk_vars(dlc_whisk_time,dlc_whisk_angle);
    
    % Convert to proper time
    [~,~,whisk_amp,loco,speed] = ...
        convert_behav_vars(time,0.2,whisk_time,whisk_angle,whisk_set_point,whisk_amp,wheel_MI,SpeedTimeMatrix,SpeedDataMatrix);
    
    %[dFF,time,acquisition_rate] = load_data(dataset_ix);
    
    % Get A and QW states
    [A,QW] = define_behav_periods(whisk_amp,loco,acq_rate,1); 
    pause
    [change_dFF{dataset_ix},p_val{dataset_ix}] = change_dFF_sig(dFF_keep,A,QW,acq_rate);
    
end
%%
N = size(change_dFF{dataset_ix},1);
ix_ON = find(change_dFF{dataset_ix} > 0 & p_val{dataset_ix} < .05);
ix_OFF = find(change_dFF{dataset_ix} < 0 & p_val{dataset_ix} < .05);

numel(ix_ON)/N
numel(ix_OFF)/N

%%

C_wsp = cell(6,1);
p_wsp = cell(6,1);

C_wamp = cell(6,1);
p_wamp = cell(6,1);

C_loco = cell(6,1);
p_loco = cell(6,1);

C_speed = cell(6,1);
p_speed = cell(6,1);

C_pupil = cell(6,1);
p_pupil = cell(6,1);

for dataset_ix = 1:6
    dataset_ix
    fname = datasets{dataset_ix};
    load([basedir,fname,'.mat'])
    
    time = mean(timeseries,2)/1000;
    
    smooth_win_s = .5;
    smooth_bins = round(acq_rate * smooth_win_s);
    
    dFF = dFF_keep;
    for k = 1:size(dFF,1)
         dFF(k,:) = smoothdata(dFF(k,:),'gaussian', [2*smooth_bins 0]);
    end
    
    % Whisking variables
    [whisk_time,whisk_angle,whisk_set_point,whisk_amp,whisk_phase] ...
        = convert_dlc_to_whisk_vars(dlc_whisk_time,dlc_whisk_angle);
    
    % Convert to proper time
    [~,whisk_set_point,whisk_amp,loco,speed] = ...
        convert_behav_vars(time,0.2,whisk_time,whisk_angle,whisk_set_point,whisk_amp,wheel_MI,SpeedTimeMatrix,SpeedDataMatrix);

    % Get A and QW states
    [C_wsp{dataset_ix},p_wsp{dataset_ix},...
        C_wamp{dataset_ix},p_wamp{dataset_ix},...
        C_loco{dataset_ix},p_loco{dataset_ix},...
        C_speed{dataset_ix},p_speed{dataset_ix},...
        C_pupil{dataset_ix},p_pupil{dataset_ix}] ...
        = corr_sig(dFF,whisk_set_point,whisk_amp,loco,speed,[],acq_rate);
end
%%


% Replace below with whisker set point, amplitude, locomotion, or pupil
C = C_speed; p_val = p_speed;
%C = C_loco; p_val = p_loco;
%C = C_wsp; p_val = p_wsp;
%C = C_wamp; p_val = p_wamp;

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
