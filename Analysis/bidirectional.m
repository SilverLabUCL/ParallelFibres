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

bins = linspace(-1,1,30);
bins_c = bins(2:end)-mean(diff(bins))/2;

% Plot whisker set point
for k = 1:4
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

%% Plot dFF

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
