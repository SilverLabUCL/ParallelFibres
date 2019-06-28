% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

datasets = {'FL87_180501_11_03_09',...  1
            'FL87_180501_10_47_25',...  2
            'FL87_180501_10_36_14',...  3
            'FL87_180220_10_38_55',...  4
            'FL77_180213_10_46_41',...  5
            'FL_S_170906_11_26_25',...  6
            'FL_S_170905_10_40_52',...  7
            'FL45_170125_14_47_04'}; %  8
        
C_all_mi = [];C_pass_shuff_mi=[];C_fail_shuff_mi=[];
C_all_speed = [];C_pass_shuff_speed=[];C_fail_shuff_speed=[];

for dataset_ix =  1:8
    
    % Load data
    [dFF,acquisition_rate] = load_data(dataset_ix);
    
    % Load behavioural data
    = load_behav(datset_ix,acquisition_rate);

    
    dataset_ix
    load(strcat('for_Alex/',dataset_ixs_Crus1_patches{dataset_ix},'.mat')) % 1,2,4,5,6
    load(strcat('for_Alex/',dataset_ixs_Crus1_patches{dataset_ix},'_MIwheel.mat')) % 1,2,4,5,6

    if dataset_ix ~=2
        load(strcat('correct_baseline_',int2str(dataset_ix),'.mat'))
    end

    [T,N] = size(Axon_dFF);
    
    Axon_dFF_ = zeros(size(Axon_dFF));
    Axon_dFF_s_ = zeros(size(Axon_dFF));
    for k = 1:size(Axon_dFF,2)
        t= TimeAxon(:,k);
        if dataset_ix ~=2
            Axon_dFF_(:,k) = Axon_dFF(:,k)-(f.a*exp(f.b*t)+f.c*exp(f.d*t));
        else
            Axon_dFF_(:,k) = Axon_dFF(:,k);
        end
        Axon_dFF_s_(:,k) = smoothdata(Axon_dFF_(:,k),'gaussian',round(200*acquisition_rate/1000));
    end
    Axon_dFF =Axon_dFF_;
    Axon_dFF_smooth =Axon_dFF_s_;
    clear Axon_dFF_
    clear Axon_dFF_s_

    dt_MI = mean(diff(MI_facepad(:,2)));
    MI_smooth = smoothdata(MI_facepad(:,1),'gaussian',round(200/dt_MI));
    MI_interp = interp1(MI_facepad(:,2),MI_smooth,TimeAxon(:,end));
    
    % wheel mi
    ix_to_correct = find(diff(wheel_MI(:,2))==0);
    for ix = ix_to_correct
        wheel_MI(ix,1) = (wheel_MI(ix,1)+wheel_MI(ix+1,1))/2;
        wheel_MI(ix+1,:) = [];
    end
    dt_speed = mean(diff(wheel_MI(:,2)));
    Speed_smooth = smoothdata(wheel_MI(:,1),'gaussian',round(200/dt_speed));
    Speed_interp = interp1(wheel_MI(:,2),Speed_smooth,TimeAxon(:,end));


%     encoder
%     dt_speed = mean(diff(SpeedTimeMatrix));
%     Speed_smooth = smoothdata(SpeedDataMatrix,'gaussian',round(200/dt_speed));
%     Speed_interp = interp1(SpeedTimeMatrix,Speed_smooth,TimeAxon(:,end));


    C_mi = zeros(N,1); C_speed = zeros(N,1);
    for n = 1:N
        C_ = corrcoef(MI_interp,Axon_dFF_smooth(:,n));
        C_mi(n) = C_(1,2);
        
        C_ = corrcoef(Speed_interp,Axon_dFF_smooth(:,n));
        C_speed(n) = C_(1,2);
        
        num_reps = 1000;
        C_sh_mi = zeros(num_reps,1);
        C_sh_speed = zeros(num_reps,1);
        for rep = 1:num_reps
            shift = floor(rand*(T-2*acquisition_rate)+2*acquisition_rate);
            dFF_temp = zeros(1,T);
            dFF_temp(shift:T) = Axon_dFF_smooth(1:T-shift+1,n);
            dFF_temp(1:shift-1) = Axon_dFF_smooth(T-shift+2:end,n);  
            C_ = corrcoef(MI_interp,dFF_temp);
            C_sh_mi(rep)= C_(1,2);
            C_ = corrcoef(Speed_interp,dFF_temp);
            C_sh_speed(rep)= C_(1,2);
        end
        p = sum(C_sh_mi > abs(C_mi(n)) |  C_sh_mi < -abs(C_mi(n)))/num_reps;
        if p<.05
            C_pass_shuff_mi = [C_pass_shuff_mi; C_mi(n)];
            
            
%             min_trial_length = floor(min(trial_times(:,2) - trial_times(:,1))*acquisition_rate/1000);
% 
%             xcorr_temp
%             for trial = 1:Numb_trials
%                 ix = find(TimeAxon_avg >= trial_times(trial,1) & TimeAxon_avg < trial_times(trial,2));
%                 ix = ix(2:min_trial_length);
%                 x = score(ix,pc_ix);
%                 y = MI_interp(ix);
% 
%                 C = corrcoef(x,y);
%                 if C(1,2)<0
%                     x = -x;
%                 end
% 
%                 [xcorr_temp, lags_temp] = xcov(x,y,'coeff');
%                 lags = [lags, lags_temp'];
% 
%             end
            
        else
            C_fail_shuff_mi = [C_fail_shuff_mi; C_mi(n)];
        end
        p = sum(C_sh_speed > abs(C_speed(n)) |  C_sh_speed < -abs(C_speed(n)))/num_reps;
        if p<.05
            C_pass_shuff_speed = [C_pass_shuff_speed; C_speed(n)];
        else
            C_fail_shuff_speed = [C_fail_shuff_speed; C_speed(n)];
        end
    end
    
    C_all_speed = [C_all_speed;C_speed];
    C_all_mi = [C_all_mi;C_mi];
    
    
end

%figure, histogram(C_all)
%%
bins = linspace(-1,1,30);
bins_c = bins(2:end)-mean(diff(bins))/2;
h = histogram(C_pass_shuff_speed,bins);
h_pass_speed = h.Values;
h = histogram(C_fail_shuff_speed,bins);
h_fail_speed = h.Values;

h = histogram(C_pass_shuff_mi,bins);
h_pass_mi = h.Values;
h = histogram(C_fail_shuff_mi,bins);
h_fail_mi = h.Values;

figure, bar(bins_c,h_fail_mi,'FaceColor',[.6,.8,1],'EdgeColor','w')
hold on, bar(bins_c,h_pass_mi,'FaceColor','b','EdgeColor','w')

figure, bar(bins_c,h_fail_speed,'FaceColor',[1,.8,.8],'EdgeColor','w')
hold on, bar(bins_c,h_pass_speed,'FaceColor','r','EdgeColor','w')

%% Plot dFF

figure, hold on
for k = 1:Nf
    plot(TimeAxon(:,k),Axon_dFF_smooth(:,k)+k,'k')
end

%% With lag, no shuffle

C_mi_max = [];

for dataset_ix =  1:7
    dataset_ix
    load(strcat('for_Alex/',dataset_ixs_Crus1_patches{dataset_ix})) % 1,2,4,5,6

    if dataset_ix ~=2
        load(strcat('correct_baseline_',int2str(dataset_ix),'.mat'))
    end

    [T,N] = size(Axon_dFF);
    
    Axon_dFF_ = zeros(size(Axon_dFF));
    Axon_dFF_s_ = zeros(size(Axon_dFF));
    for k = 1:size(Axon_dFF,2)
        t= TimeAxon(:,k);
        if dataset_ix ~=2
            Axon_dFF_(:,k) = Axon_dFF(:,k)-(f.a*exp(f.b*t)+f.c*exp(f.d*t));
        else
            Axon_dFF_(:,k) = Axon_dFF(:,k);
        end
        Axon_dFF_s_(:,k) = smoothdata(Axon_dFF_(:,k),'gaussian',round(200*acquisition_rate/1000));
    end
    Axon_dFF =Axon_dFF_;
    Axon_dFF_smooth =Axon_dFF_s_;
    clear Axon_dFF_
    clear Axon_dFF_s_

    dt_MI = mean(diff(MI_facepad(:,2)));
    MI_smooth = smoothdata(MI_facepad(:,1),'gaussian',round(200/dt_MI));
    MI_interp = interp1(MI_facepad(:,2),MI_smooth,TimeAxon(:,end));

    dt_speed = mean(diff(SpeedTimeMatrix));
    Speed_smooth = smoothdata(SpeedDataMatrix,'gaussian',round(200/dt_speed));
    Speed_interp = interp1(SpeedTimeMatrix,Speed_smooth,TimeAxon(:,end));


    x = diff(mean(TimeAxon,2)); trial_switch = find(abs(x-median(x))>.01);
    if numel(trial_switch)~=29
        error
    end
    
    
    
    C_mi = zeros(N, 2*(T/30)-1); C_speed = zeros(N, 2*(T/30)-1);
    
    for n = 1:N
        temp_mi = zeros(30,2*(T/30)-1);
        temp_speed = zeros(30,2*(T/30)-1);
        
        for tr = 1:30
            if tr==1
                ix_trial = 1:trial_switch(1);
            elseif tr == 30
                ix_trial = trial_switch(end)+1:T;
            else
                ix_trial = trial_switch(tr-1)+1:trial_switch(tr);
            end
            [temp_mi(tr,:), lag] = xcov(MI_interp(ix_trial),Axon_dFF_smooth(ix_trial,n),'coeff');
            [temp_speed(tr,:), lag_speed] = xcov(Speed_interp(ix_trial),Axon_dFF_smooth(ix_trial,n),'coeff');
            if norm(lag-lag_speed)~=0
                error
            end
        end

        C_mi(n,:) = mean(temp_mi,1);

        C_speed(n,:) = mean(temp_speed,1);
        
        [max_val,b] = max(abs(C_mi(n,:)));
        if max_val > 3*std(C_mi(n,:)) %abs(lag(b))/acquisition_rate < 1 & max_val > 3*std(C_mi(n,:))
            C_mi_max = [C_mi_max; (lag(b) + rand - 0.5 )/acquisition_rate, C_mi(n,b)];
            if max_val > 4*std(C_mi(n,:)) & abs(lag(b))/acquisition_rate > .1 & max_val > .4%.8 & C_mi(n,b)<0
                error
            end
        end
    end
    
    %figure, plot(lag/acquisition_rate,C_mi)
    %xlabel('Lag (s)'); ylabel('Correlation with WMI')
    %title(dataset_ix)
    
    
end

%figure, histogram(C_all)
