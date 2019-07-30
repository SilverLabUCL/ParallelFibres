% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;

%%

[C_wsp,C_amp,C_spd,C_pupil] = deal(nan(7,1));

[C_wsp_Ac,C_amp_Ac,C_spd_Ac,C_pupil_Ac] = deal(nan(7,1));

[C_wsp_Q,C_amp_Q,C_spd_Q,C_pupil_Q] = deal(nan(7,1));

for dataset_ix = 1:7
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get score
    [coeff, score] = pca(dFF');
    
    % Correlation with pc1
    pc1 = score(:,1);
    temp = corr(whisk_set_point,pc1);
    
    if temp < 0
        pc1 = -pc1;
        C_wsp(dataset_ix) = -temp;
    else
        C_wsp(dataset_ix) = temp;
    end
    
    get_behavioural_periods
    
    C_amp(dataset_ix) = corr(whisk_amp,pc1);
    C_spd(dataset_ix) = corr(speed,pc1);
    if ~isempty(pupil)
        C_pupil(dataset_ix) = corr(pupil,pc1);
    end
    
    
    get_active_per_hack
    %pause
    
    [pc1_Ac,whisk_set_point_Ac,whisk_amp_Ac,pupil_Ac,speed_Ac] = deal([]);
    for k = 1:length(aroused)
        ix = find(time >= aroused(k,1) & time <= aroused(k,2));
        pc1_Ac = [pc1_Ac; pc1(ix)];
        whisk_set_point_Ac = [whisk_set_point_Ac; whisk_set_point(ix)];
        whisk_amp_Ac = [whisk_amp_Ac; whisk_amp(ix)];
        speed_Ac = [speed_Ac; speed(ix)];
        if ~isempty(pupil)
            pupil_Ac = [pupil_Ac; pupil(ix)];
        end
    end
        
    C_wsp_Ac(dataset_ix) = corr(whisk_set_point_Ac,pc1_Ac);
    C_amp_Ac(dataset_ix) = corr(whisk_amp_Ac,pc1_Ac);
    C_spd_Ac(dataset_ix) = corr(speed_Ac,pc1_Ac);
    if ~isempty(pupil)
        C_pupil_Ac(dataset_ix) = corr(pupil_Ac,pc1_Ac);
    end
    
    [pc1_Q,whisk_set_point_Q,whisk_amp_Q,pupil_Q,speed_Q] = deal([]);
    for k = 1:length(quiescent)
        ix = find(time >= quiescent(k,1) & time <= quiescent(k,2));
        pc1_Q = [pc1_Q; pc1(ix)];
        whisk_set_point_Q = [whisk_set_point_Q; whisk_set_point(ix)];
        whisk_amp_Q = [whisk_amp_Q; whisk_amp(ix)];
        speed_Q = [speed_Q; speed(ix)];
        if ~isempty(pupil)
            pupil_Q = [pupil_Q; pupil(ix)];
        end
    end
        
    C_wsp_Q(dataset_ix) = corr(whisk_set_point_Q,pc1_Q);
    C_amp_Q(dataset_ix) = corr(whisk_amp_Q,pc1_Q);
    C_spd_Q(dataset_ix) = corr(speed_Q,pc1_Q);
    if ~isempty(pupil)
        C_pupil_Q(dataset_ix) = corr(pupil_Q,pc1_Q);
    end
    
    

end
%%
figure,  hold on
plot(zeros,C_wsp,'ob','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(ones,C_wsp_Q,'oc','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(2*ones,C_wsp_Ac,'om','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
set(gca,'Xtick',[0,1,2])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Cross correlation'))
xlim([-.5,2.5])
ylim([-.5,1])
%%


figure,  hold on
plot(zeros,C_amp,'ob','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(ones,C_amp_Q,'oc','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(2*ones,C_amp_Ac,'om','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
set(gca,'Xtick',[0,1,2])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Cross correlation'))
xlim([-.5,2.5])
ylim([-.5,1])







%% ALL FOLLOWING CODE FOR variability
%%
[offset_indices, onset_indices] = get_onsets(whisk_amp,acquisition_rate);

figure, x = whisk_amp;
plot(x,'k'), hold on
for k = 1:size(onset_indices,2)
    t = offset_indices(k):onset_indices(k);
    plot(t,x(t),'c')
    plot(onset_indices(k),x(onset_indices(k)),'ok','MarkerFaceColor','c')
end
size(onset_indices,2)
%%
[N,T] = size(dFF);

% [U,S,V] = svd(dFF');
% S(1,1)=0;
% dFF = (U*S*V)';

[coeff, score] = pca(dFF');
dFF = score';

num_bins = round(.2 * acquisition_rate);

dFF_pre_onset = zeros(N,1);

for t = 1:size(onset_indices,2)
    t_ = ((onset_indices(t)-num_bins)-1) : (onset_indices(t)-1);
    dFF_pre_onset = [dFF_pre_onset,dFF(:,t_)];
end

E = eig(cov(dFF_pre_onset')); L = sqrt(sum(E.^2))


num_reps = 1000;

L_shuff = zeros(num_reps,1);

for i = 1:num_reps
    dFF_shuff = zeros(N,1);
    
    for t = 1:size(onset_indices,2)
        t_start = randsample(offset_indices(t): ((onset_indices(t)-8*num_bins)-1),1);
        t_ = t_start:(t_start+num_bins);
        dFF_shuff = [dFF_shuff,dFF(:,t_)];
    end
    
    E = eig(cov(dFF_shuff'));
    L_shuff(i) = sqrt(sum(E.^2));
end

figure, histogram(L_shuff,'Normalization','probability')
hold on, plot([L,L],ylim,':k','LineWidth',2)
set(gca,'FontSize',15)
xlabel('Total variance')
ylabel('Probability')

p = sum(L_shuff < L)/num_reps

%% Plot 3d manifold
[coeff, score] = pca(dFF');

zsp = zscore(whisk_set_point);

figure, hold on
for k = 1:length(score)-1
    c = (zsp(k)+1)/2;
    if c < 0
        c = 0;
    elseif c > 1
        c = 1;
    end
    c=c*[1,0,1]+(1-c)*[0,1,1];
   plot3(score(k:k+1,1),score(k:k+1,2),score(k:k+1,3),'-','Color',c)
end
view(-53,-52)

plot3(score(onset_indices,1),score(onset_indices,2),score(onset_indices,3),'or')
plot3(score(offset_indices,1),score(offset_indices,2),score(offset_indices,3),'ok')

%%

figure, plot(zscore(whisk_angle))
hold on, plot(zscore(score(:,1))+3)
hold on, plot(zscore(score(:,2))+6)
hold on, plot(zscore(score(:,3))+9)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot PC 1 vs coding dimension - summary plots

ang_cd_pc1 = zeros(7,1);
rho_cd_pc1 = zeros(7,1);

for dataset_ix = 1:7
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get score
    [coeff, score] = pca(dFF');

    % Swap sign if negatively correlated with set point
    if corr(whisk_set_point,score(:,1)) < 0
        coeff(:,1) = - coeff(:,1);
    end
    
    pc1 = coeff(:,1);
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
    cd = get_coding_dimension(dFF,A,QW);
    
    if round(norm(coeff(:,1)),10) ~= 1 || round(norm(w),10) ~=1
        error
    end
    
    ang_cd_pc1(dataset_ix) = subspace(pc1,cd);
    rho_cd_pc1(dataset_ix) = corr(pc1,cd);
end

figure, plot(ones(7,1),rho_cd_pc1,'ok','MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.5)
ylabel('Correlation'), ylim([0,1])
set(gca,'FontSize',15,'XTick',[],'Box','off')
signrank(rho_cd_pc1 - 1.5)

figure, plot(ones(7,1),ang_cd_pc1,'ok','MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.5)
ylabel('Angle (rad.)'), ylim([0,1.5])
set(gca,'FontSize',15,'XTick',[],'Box','off')
signrank(ang_cd_pc1 - 1.5)

%% Plot figures of PC 1 vs coding dimension - example dataset

dataset_ix = 6; % set to 'worst' example (lowest corr)

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
pupil = load_pupil(dataset_ix,time);

% Get score
[coeff, score] = pca(dFF');

% Swap sign if negatively correlated with set point
if corr(whisk_set_point,score(:,1)) < 0
    coeff(:,1) = - coeff(:,1);
else
end

pc1 = coeff(:,1);
[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,1);
cd = get_coding_dimension(dFF,A,QW);

figure, plot(cd,'r','LineWidth',2)
hold on, plot(pc1,':b','LineWidth',2)
xlim([1,size(dFF,1)])
set(gca,'FontSize',15,'Box','off')
xlabel('Axon Number')
ylabel('Coefficient')

corr(cd,pc1)