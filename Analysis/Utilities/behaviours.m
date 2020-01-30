% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;

%% Example behaviours
% Figure S3B

dataset_ix = 1;
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
pupil = load_pupil(dataset_ix,time);

figure, plot(time,whisk_set_point,'k'), axis([0,170,-10,40]), xlabel('Time (s)'), ylabel('Whisker set point')
figure, plot(time,whisk_amp,'k'), axis([0,170,-2,25]), xlabel('Time (s)'), ylabel('Whisker amplitude')
figure, plot(time,speed,'k'), axis([0,170,-2,10]), xlabel('Time (s)'), ylabel('Locomotion')
figure, plot(time,pupil,'k'), axis([0,170,2000,2850]), xlabel('Time (s)'), ylabel('Pupil')


%% Correlation between behaviours
% Figure S3C

C_wsp_wa = nan(15,1);
p_wsp_wa = nan(15,1);
C_wsp_wa_shuff = nan(15,1);

C_wsp_l = nan(15,1);
p_wsp_l = nan(15,1);
C_wsp_l_shuff = nan(15,1);

C_wsp_p = nan(15,1);
p_wsp_p = nan(15,1);
C_wsp_p_shuff = nan(15,1);

C_wa_l = nan(15,1);
p_wa_l = nan(15,1);
C_wa_l_shuff = nan(15,1);

C_wa_p = nan(15,1);
p_wa_p = nan(15,1);
C_wa_p_shuff = nan(15,1);

C_l_p = nan(15,1);
p_l_p = nan(15,1);
C_l_p_shuff = nan(15,1);

for dataset_ix = 1:15
    dataset_ix
    [whisk_angle,whisk_set_point,whisk_amp,speed,whisk_time,speed_time]...
        = load_behav_data(dataset_ix);
    [pupil, pupil_time] = load_pupil(dataset_ix);

    acquisition_rate = 100;
    
    time_0 = max([whisk_time(1),speed_time(1)]);
    time_end = min([whisk_time(end),speed_time(end)]);
    if ~isempty(pupil)
        time_0 = max([time_0,pupil_time(1)]);
        time_end = min([time_end,pupil_time(end)]);
    end
    time = (time_0:(1/acquisition_rate):time_end)';
    
    whisk_set_point = interp1(whisk_time,whisk_set_point,time);
    whisk_amp = interp1(whisk_time,whisk_amp,time);
    speed = interp1(speed_time,speed,time);
    if ~isempty(pupil)
        pupil = interp1(pupil_time,pupil,time);
    end
    
    [C_wsp_wa(dataset_ix),p_wsp_wa(dataset_ix),C_wsp_wa_shuff(dataset_ix)] = corr_behav_sig(whisk_set_point,whisk_amp,acquisition_rate);
    [C_wsp_l(dataset_ix),p_wsp_l(dataset_ix),C_wsp_l_shuff(dataset_ix)] = corr_behav_sig(whisk_set_point,speed,acquisition_rate);
    
    [C_wa_l(dataset_ix),p_wa_l(dataset_ix),C_wa_l_shuff(dataset_ix)] = corr_behav_sig(whisk_amp,speed,acquisition_rate);
    
    if ~isempty(pupil)
        [C_wsp_p(dataset_ix),p_wsp_p(dataset_ix),C_wsp_p_shuff(dataset_ix)] = corr_behav_sig(whisk_set_point,pupil,acquisition_rate);
        [C_wa_p(dataset_ix),p_wa_p(dataset_ix),C_wa_p_shuff(dataset_ix)] = corr_behav_sig(whisk_amp,pupil,acquisition_rate);
        [C_l_p(dataset_ix),p_l_p(dataset_ix),C_l_p_shuff(dataset_ix)] = corr_behav_sig(speed,pupil,acquisition_rate);
    end
end

%% Plot behavioural correlations 
% Figure S3D 

c = [.5,.5,.5];

figure, hold on
plot(zeros(15,1),C_wsp_wa,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones(15,1),C_wsp_wa_shuff,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(0+.2*[-1,1],mean(C_wsp_wa)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],mean(C_wsp_wa_shuff)*[1,1],'k','LineWidth',3)

plot(2.5*ones(15,1),C_wsp_l,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(3.5*ones(15,1),C_wsp_l_shuff,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(2.5+.2*[-1,1],mean(C_wsp_l)*[1,1],'k','LineWidth',3)
plot(3.5+.2*[-1,1],mean(C_wsp_l_shuff)*[1,1],'k','LineWidth',3)

plot(5*ones(15,1),C_wa_l,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(6*ones(15,1),C_wa_l_shuff,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(5+.2*[-1,1],mean(C_wa_l)*[1,1],'k','LineWidth',3)
plot(6+.2*[-1,1],mean(C_wa_l_shuff)*[1,1],'k','LineWidth',3)

xlim([-.5,6.5])
set(gca,'FontSize',15,'XTick',[0,1,2.5,3.5,5,6],'XTickLabel',[])
ylabel('Correlation')

signrank(C_wsp_wa,C_wsp_wa_shuff)
signrank(C_wsp_l,C_wsp_l_shuff)
signrank(C_wa_l,C_wa_l_shuff)
%% Correlations with pupil
% Figure S3D 
figure, hold on
plot(zeros(15,1),C_wsp_p,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones(15,1),C_wsp_p_shuff,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(0+.2*[-1,1],nanmean(C_wsp_p)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],nanmean(C_wsp_p_shuff)*[1,1],'k','LineWidth',3)

plot(2.5*ones(15,1),C_wa_p,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(3.5*ones(15,1),C_wa_p_shuff,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(2.5+.2*[-1,1],nanmean(C_wa_p)*[1,1],'k','LineWidth',3)
plot(3.5+.2*[-1,1],nanmean(C_wa_p_shuff)*[1,1],'k','LineWidth',3)

plot(5*ones(15,1),C_l_p,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(6*ones(15,1),C_l_p_shuff,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(5+.2*[-1,1],nanmean(C_l_p)*[1,1],'k','LineWidth',3)
plot(6+.2*[-1,1],nanmean(C_l_p_shuff)*[1,1],'k','LineWidth',3)

xlim([-.5,6.5])
set(gca,'FontSize',15,'XTick',[0,1,2.5,3.5,5,6],'XTickLabel',[])
ylabel('Correlation')