clear all; clc

datasets_Crus1_POIs = {...
'Crus1_POIs_FL77_180207_10_12_41.mat',... % 117*b
'Crus1_POIs_FL77_180207_10_42_44.mat',... % 40 *
'Crus1_POIs_FL77_180207_11_01_27.mat',... % 23 *b
'Crus1_POIs_FL87_180220_10_52_03.mat',... % 35
'Crus1_POIs_FL87_180413_11_40_37.mat'};   % 45

datasets_Crus1_patches = {...
'Crus1_patches_FL77_180213_10_46_41',... % 158 *b
'Crus1_patches_FL_S_170905_10_40_52',... % 120 *
'Crus1_patches_FL87_180220_10_38_55',... % 77 // *b
'Crus1_patches_FL87_180501_10_36_14',... % 150 *
'Crus1_patches_FL87_180501_10_47_25',... % 150 *
'Crus1_patches_FL87_180501_11_03_09',...   % 169 *
'Crus1_patches_FL_S_170906_11_26_25'};

datasets_Lob4_5_POIs = {...
'Lob4_5_POIs_FL42_161109_10_38_39.mat',... % 13 *
'Lob4_5_POIs_FL93_180226_11_08_23.mat'};   % 55

datasets_Lob4_5_patches = {...
'Lob4_5_patches_FL93_180424_10_45_32.mat',... % 40
'Lob4_5_patches_FL93_180502_11_13_55.mat',... % 38
'Lob4_5_patches_FL93_180724_10_51_11.mat',... % 21
'Lob4_5_patches_FL42_161109_10_47_57.mat',... % 33 *
'Lob4_5_patches_FL99_180523_10_18_32.mat',... % 108 *W
'Lob4_5_patches_FL93_180424_10_30_35.mat'};   % 45

xcorr_speed_all = nan(7,1);
xcorr_wmi_all = nan(7,1);
xcorr_pupil_all = nan(7,1);

xcorr_speed_active = nan(7,1);
xcorr_wmi_active = nan(7,1);
xcorr_pupil_active = nan(7,1);

for dataset =  1:7
    %%
    
    dataset
    
    load(strcat('for_Alex/',datasets_Crus1_patches{dataset},'.mat')) % 1,2,4,5,6
    load(strcat('for_Alex/',datasets_Crus1_patches{dataset},'_MIwheel.mat')) % 1,2,4,5,6

    N = size(Axon_dFF,2);
     
    if dataset ~=2
        load(strcat('correct_baseline_',int2str(dataset),'.mat'))


        Axon_dFF_ = zeros(size(Axon_dFF));
        Axon_dFF_s_ = zeros(size(Axon_dFF));
        for k = 1:size(Axon_dFF,2)
            t= TimeAxon(:,k);
            Axon_dFF_(:,k) = Axon_dFF(:,k)-(f.a*exp(f.b*t)+f.c*exp(f.d*t));
            Axon_dFF_s_(:,k) = smoothdata(Axon_dFF_(:,k),'gaussian',round(200*acquisition_rate/1000));
        end
        Axon_dFF =Axon_dFF_;
        Axon_dFF_smooth =Axon_dFF_s_;
        clear Axon_dFF_
        clear Axon_dFF_s_
    end

    [coeff, score, ~, ~, explained] = pca(Axon_dFF_smooth);
    
    % Resample behaviour
    TimeAxon_avg =  mean(TimeAxon,2);

    dt_MI = mean(diff(MI_facepad(:,2)));
    MI_smooth = smoothdata(MI_facepad(:,1),'gaussian',round(200/dt_MI));
    MI_interp = interp1(MI_facepad(:,2),MI_smooth,TimeAxon_avg);
    
    C = corrcoef(score(:,1),MI_interp);
    if C(1,2)<0
        score(:,1) = -score(:,1);
    end

    % wheel mi
    ix_to_correct = find(diff(wheel_MI(:,2))==0);
    for ix = ix_to_correct
        wheel_MI(ix,1) = (wheel_MI(ix,1)+wheel_MI(ix+1,1))/2;
        wheel_MI(ix+1,:) = [];
    end
    dt_speed = mean(diff(wheel_MI(:,2)));
    Speed_smooth = smoothdata(wheel_MI(:,1),'gaussian',round(200/dt_speed));
    Speed_interp = interp1(wheel_MI(:,2),Speed_smooth,TimeAxon(:,end));

    % encoder
%     dt_speed = mean(diff(SpeedTimeMatrix));
%     Speed_smooth = smoothdata(SpeedDataMatrix,'gaussian',round(200/dt_speed));
%     Speed_interp = interp1(SpeedTimeMatrix,Speed_smooth,TimeAxon_avg);

    try
        dt_pupil = mean(diff(Pupil(:,1)));
        Pupil_smooth = smoothdata(Pupil(:,2),'gaussian',round(200/dt_pupil));
        Pupil_interp = interp1(Pupil(:,1),Pupil_smooth,TimeAxon_avg);
    end
    find_arousal_periods_pc;
    
    %% Get times of trials

    ix = find(abs(diff(TimeAxon_avg)-1000/acquisition_rate)>.001);
    if numel(ix)~=Numb_trials-1
        error
    end

    trial_times = zeros(Numb_trials,2);
    trial_times(1,1) = TimeAxon_avg(1);
    trial_times(2:end,1) = TimeAxon_avg(ix+1);
    trial_times(1:end-1,2) = TimeAxon_avg(ix-1);
    trial_times(end,end) = TimeAxon_avg(end);

    %% Get correlations by trial

    xcorr_mi = [];
    xcorr_speed = [];
    xcorr_pupil = [];
    lags = [];

    pc_ix=1;

    min_trial_length = floor(min(trial_times(:,2) - trial_times(:,1))*acquisition_rate/1000);

    for trial = 1:Numb_trials
        ix = find(TimeAxon_avg >= trial_times(trial,1) & TimeAxon_avg < trial_times(trial,2));
        ix = ix(2:min_trial_length);
        x = score(ix,pc_ix);
        y = MI_interp(ix);

        C = corrcoef(x,y);
        if C(1,2)<0
            x = -x;
        end

        [xcorr_temp, lags_temp] = xcov(x,y,'coeff');
        xcorr_mi = [xcorr_mi, xcorr_temp];
        lags = [lags, lags_temp'];

        y = Speed_interp(ix);
        [xcorr_temp, lags_temp] = xcov(x,y,'coeff');
        xcorr_speed = [xcorr_speed, xcorr_temp];
        lags = [lags, lags_temp'];
        
        try
            y = Pupil_interp(ix);
            [xcorr_temp, lags_temp] = xcov(x,y,'coeff');
            xcorr_pupil = [xcorr_pupil, xcorr_temp];
            lags = [lags, lags_temp'];
        end

    end
    
    if mean(var(lags,[],2))~=0
        error
    else 
        lags = lags(:,1)/acquisition_rate;
    end

    %lag_choice_ix = find(lags==0);
    lag_choice_ix = find(mean(xcorr_mi,2) == max(mean(xcorr_mi,2)));
    lag_choice = lags(lag_choice_ix);
    
    xcorr_speed_all(dataset) = nanmean(xcorr_speed(lag_choice_ix,:));
    xcorr_wmi_all(dataset) =  nanmean(xcorr_mi(lag_choice_ix,:));
    
    try
       xcorr_pupil_all(dataset) = nanmean(xcorr_pupil(lag_choice_ix,:));
    end
   
    %% Get correlations w/in arousal period for each trial

    Time_aroused = [];
    Axon_dFF_smooth_aroused = [];
    MI_interp_a = [];
    Pupil_interp_a = [];
    Speed_interp_a = [];
    score_orig_a = [];
    for k = 1:length(aroused)

        min_length = zeros(N,1);
        for j = 1:N
            min_length(j) = length(find(TimeAxon(:,j) >= aroused(k,1) & TimeAxon(:,j) < aroused(k,2)));
        end
        min_length = min(min_length);

        temp = zeros(min_length,N);
        ix = find(TimeAxon_avg >= aroused(k,1) & TimeAxon_avg < aroused(k,2));
        temp = Axon_dFF_smooth(ix(1:min_length),:);

        Axon_dFF_smooth_aroused = [Axon_dFF_smooth_aroused; temp];

        ix = find(TimeAxon_avg >= aroused(k,1) & TimeAxon_avg < aroused(k,2));
        Time_aroused = [Time_aroused; TimeAxon_avg(ix(1:min_length))];
        MI_interp_a = [MI_interp_a; MI_interp(ix)];
        Speed_interp_a = [Speed_interp_a; Speed_interp(ix)];
        try
            Pupil_interp_a = [Pupil_interp_a; Pupil_interp(ix)];
        end
        score_orig_a = [score_orig_a;score(ix)];
    end
    
    %% Get correlations by trial

    xcorr_mi = [];
    xcorr_speed = [];
    xcorr_pupil = [];
    lags = [];

    pc_ix=1;

    for trial = 1:Numb_trials

        ix = find(Time_aroused >= trial_times(trial,1) & Time_aroused < trial_times(trial,2));
        x = score_orig_a(ix,pc_ix);
            
        if numel(ix) > 50
            
            y = Speed_interp_a(ix);
            [xcorr_temp, lags_temp] = xcov(x,y,100,'coeff');
            xcorr_mi = [xcorr_mi, xcorr_temp];
            lags = [lags, lags_temp'];

            y = MI_interp_a(ix); 
            [xcorr_temp, lags_temp] = xcov(x,y,100,'coeff');
            xcorr_speed = [xcorr_speed, xcorr_temp];
            lags = [lags, lags_temp'];
        
            try
                y = Pupil_interp_a(ix);
                [xcorr_temp, lags_temp] = xcov(x,y,100,'coeff');
                xcorr_pupil = [xcorr_pupil, xcorr_temp];
                lags = [lags, lags_temp'];
            end

        end

    end

    if mean(var(lags,[],2))~=0
        error
    else 
        lags = lags(:,1)/acquisition_rate;
    end
    
    lag_choice_ix = find(lags>=lag_choice);
    lag_choice_ix = lag_choice_ix(1);
    
    xcorr_speed_active(dataset) =  nanmean(xcorr_speed(lag_choice_ix,:));
    xcorr_wmi_active(dataset) = nanmean(xcorr_mi(lag_choice_ix,:));
    try
        xcorr_pupil_active(dataset) = nanmean(xcorr_pupil(lag_choice_ix,:));
    end
end

%% Plot distributions of correlations

figure,  hold on
plot(ones,xcorr_wmi_all,'ob')
plot(1+.1*[1,-1],[1,1]*nanmean(xcorr_wmi_all),'-k','LineWidth',2)
plot(2*ones,xcorr_speed_all,'or')
plot(2+.1*[1,-1],[1,1]*nanmean(xcorr_speed_all),'-k','LineWidth',2)
plot(3*ones,xcorr_pupil_all,'o','Color',[.8,.65,0])
plot(3+.1*[1,-1],[1,1]*nanmean(xcorr_pupil_all),'-k','LineWidth',2)
set(gca,'Xtick',[1,2,3])
set(gca,'XtickLabel',{'WMI','Speed','Pupil'})
set(gca,'FontSize',15)
set(ylabel('Correlation coefficient'))
xlim([0,4])
ylim([-.5,1])

%%
figure,  hold on
plot(zeros,xcorr_wmi_all,'ob','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(0+.3*[1,-1],[1,1]*nanmean(xcorr_wmi_all),'-k','LineWidth',2)
plot(ones,xcorr_wmi_active,'om','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(1+.3*[1,-1],[1,1]*nanmean(xcorr_wmi_active),'-k','LineWidth',2)
set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Cross correlation'))
xlim([-.5,1.5])
ylim([-.5,1])

%%
figure,  hold on
plot(zeros,xcorr_speed_all,'or','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(0+.3*[1,-1],[1,1]*nanmean(xcorr_speed_all),'-k','LineWidth',2)
plot(ones,xcorr_speed_active,'om','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(1+.3*[1,-1],[1,1]*nanmean(xcorr_speed_active),'-k','LineWidth',2)
set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Cross correlation'))
xlim([-.5,1.5])
ylim([-.5,1])


%%
figure,  hold on
plot(zeros,xcorr_pupil_all,'o','Color',[.8,.65,0],'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(0+.3*[1,-1],[1,1]*nanmean(xcorr_pupil_all),'-k','LineWidth',2)
plot(ones,xcorr_pupil_active,'om','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(1+.3*[1,-1],[1,1]*nanmean(xcorr_pupil_active),'-k','LineWidth',2)
set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{'All','Active'})
set(gca,'FontSize',15)
set(ylabel('Cross correlation'))
xlim([-.5,1.5])
ylim([-.5,1])



