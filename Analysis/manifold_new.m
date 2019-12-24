% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;

%% Show that distance within is 

angle_A_QW = nan(17,1);
angle_shuff = nan(17,1);

dist = @(x,y) sqrt(sum((x-y).^2));

pval_1 = nan(17,1);
pval_2 = nan(17,1);
pval_3 = nan(17,1);

norm_dist_in_A = nan(17,1);
norm_dist_in_QW = nan(17,1);
norm_dist_A_QW = nan(17,1);

figure, hold on
for dataset_ix = [1:6,9,16,17]
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get score
    [~,T] = size(dFF);
    
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
    
    dFF_A = [];
    for k = 1:length(A)
        ix = A(k,1):A(k,2);
        dFF_A = [dFF_A, dFF(:,ix)];
    end
    
    dFF_QW = [];
    for k = 1:length(QW)
        ix = QW(k,1):QW(k,2);
        dFF_QW = [dFF_QW, dFF(:,ix)];
    end
    
    T_A = size(dFF_A,2);
    dist_in_A = nan(T_A,1); ix = 1;
    for t1 = 1:T_A
        for t2 = t1:T_A
            dist_in_A(ix) = dist(dFF_A(:,t1),dFF_A(:,t2));
            ix = ix+1;
        end
    end
    
    T_QW = size(dFF_QW,2);
    dist_in_QW = nan(T_QW,1); ix = 1;
    for t1 = 1:T_QW
        for t2 = t1:T_QW
            dist_in_QW(ix) = dist(dFF_QW(:,t1),dFF_QW(:,t2));
            ix = ix+1;
        end
    end

    
    dist_A_QW = nan(T_QW,1); ix = 1;
    for t1 = 1:T_A
        for t2 = 1:T_QW
            dist_A_QW(ix) = dist(dFF_A(:,t1),dFF_QW(:,t2));
            ix = ix+1;
        end
    end
    
    pval_3(dataset_ix) = ranksum((dist_in_A),(dist_in_QW));
    pval_1(dataset_ix) = ranksum((dist_in_A),(dist_A_QW));
    pval_2(dataset_ix) = ranksum((dist_in_QW),(dist_A_QW));
    
    
    norm_dist_in_A(dataset_ix) = nanmean(dist_in_A)/nanmean(dist_in_QW);
    norm_dist_in_QW(dataset_ix) = nanmean(dist_in_QW)/nanmean(dist_in_QW);
    norm_dist_A_QW(dataset_ix) = nanmean(dist_A_QW)/nanmean(dist_in_QW);
    
    plot([0,1,2],[norm_dist_in_QW(dataset_ix),norm_dist_in_A(dataset_ix),norm_dist_A_QW(dataset_ix)],'Color',[.8,.8,.8],'LineWidth',2)
    
    plot(0,norm_dist_in_QW(dataset_ix),'oc','MarkerFaceColor','w','LineWidth',2,'MarkerSize',8)
    %plot([0,0],nanmean(dist_in_A)+[-1,1]*std(dist_in_A),'-m')
    
    plot(1,norm_dist_in_A(dataset_ix),'om','MarkerFaceColor','w','LineWidth',2,'MarkerSize',8)
    %plot([1,1],nanmean(dist_in_QW)+[-1,1]*std(dist_in_QW),'-c')
    
    plot(2,norm_dist_A_QW(dataset_ix),'ok','MarkerFaceColor','w','LineWidth',2,'MarkerSize',8)
    %plot([2,2],nanmean(dist_A_QW)+[-1,1]*std(dist_A_QW),'-k')
    pause(.1)
    
end


%% Show that coding dimension captures transitions between states

dataset_ix = 1;
[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
pupil = load_pupil(dataset_ix,time);

[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,1);
cd = get_coding_dimension(dFF,A,QW);

A_or_QW = cd'*(dFF - mean(dFF,2));

figure, 
hold on, plot(cd'*(dFF - mean(dFF,2)))



%%
figure, plot(ones(7,1),rho_cd_pc1,'ok','MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.5)
ylabel('Correlation'), ylim([0,1])
set(gca,'FontSize',15,'XTick',[],'Box','off')
signrank(rho_cd_pc1 - 1.5)

figure, plot(ones(7,1),ang_cd_pc1,'ok','MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.5)
ylabel('Angle (rad.)'), ylim([0,1.5])
set(gca,'FontSize',15,'XTick',[],'Box','off')
signrank(ang_cd_pc1 - 1.5)

%% Plot figures of PC 1 vs coding dimension - example dataset
% Fig 3A,C

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

%% Figure 3G

[C_wsp,C_amp,C_spd,C_pupil] = deal(nan(7,1));

[C_wsp_A,C_amp_A,C_spd_A,C_pupil_A] = deal(nan(7,1));

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
    
    %[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,1);
    buffer = round(acquisition_rate * .5);
    [A,QW] = get_A_QW_periods(pc1,buffer,1);
    pause

    C_amp(dataset_ix) = corr(whisk_amp,pc1);
    C_spd(dataset_ix) = corr(speed,pc1);
    if ~isempty(pupil)
        C_pupil(dataset_ix) = corr(pupil,pc1);
    end
    
    
    [pc1_A,whisk_set_point_A,whisk_amp_A,pupil_A,speed_A] = deal([]);
    for k = 1:length(A)
        ix = A(k,1):A(k,2);
        pc1_A = [pc1_A; pc1(ix)];
        whisk_set_point_A = [whisk_set_point_A; whisk_set_point(ix)];
        whisk_amp_A = [whisk_amp_A; whisk_amp(ix)];
        speed_A = [speed_A; speed(ix)];
        if ~isempty(pupil)
            pupil_A = [pupil_A; pupil(ix)];
        end
    end
        
    C_wsp_A(dataset_ix) = corr(whisk_set_point_A,pc1_A);
    C_amp_A(dataset_ix) = corr(whisk_amp_A,pc1_A);
    C_spd_A(dataset_ix) = corr(speed_A,pc1_A);
    if ~isempty(pupil)
        C_pupil_A(dataset_ix) = corr(pupil_A,pc1_A);
    end
    
    [pc1_Q,whisk_set_point_Q,whisk_amp_Q,pupil_Q,speed_Q] = deal([]);
    for k = 1:length(QW)
        ix = QW(k,1):QW(k,2);
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

figure,  hold on
plot(zeros,C_wsp,'ob','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(ones,C_wsp_Q,'oc','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(2*ones,C_wsp_A,'om','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
set(gca,'Xtick',[0,1,2])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Cross correlation'))
xlim([-.5,2.5])
ylim([-.5,1])

signrank(C_wsp,C_wsp_Q)
signrank(C_wsp_Q,C_wsp_A)
signrank(C_wsp,C_wsp_A)

%% Figure 3H


%% A and QW subspaces are orthogonal

angle_A_QW = nan(17,1);
angle_shuff = nan(17,1);

num_PCs = 3;

for dataset_ix = [1:6,9,16]
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get score
    [~,T] = size(dFF);
    
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
    
    dFF_A = [];
    for k = 1:length(A)
        ix = A(k,1):A(k,2);
        dFF_A = [dFF_A, dFF(:,ix)];
    end
    coeff_A = pca(dFF_A');
    
    dFF_QW = [];
    for k = 1:length(QW)
        ix = QW(k,1):QW(k,2);
        dFF_QW = [dFF_QW, dFF(:,ix)];
    end
    coeff_QW = pca(dFF_QW');
    
    angle_A_QW(dataset_ix) = subspace(coeff_QW(:,1:num_PCs),coeff_A(:,1:num_PCs));
    
    train_ixs = block_shuffle_time(T,acquisition_rate);
    test_ixs = train_ixs(1:round(T/2));
    train_ixs = setdiff(train_ixs,test_ixs); 

    coeff_1 = pca(dFF(:,test_ixs)');
    coeff_2 = pca(dFF(:,train_ixs)');
    
    angle_shuff(dataset_ix) = subspace(coeff_1(:,1:num_PCs),coeff_2(:,1:num_PCs));
    
end

figure,  hold on
plot(zeros,angle_A_QW,'ok','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(ones,angle_shuff,'ok','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Angle'))
xlim([-.5,1.5])
ylim([0,1.6])

signrank(angle_A_QW,angle_shuff)


%% ALL FOLLOWING CODE FOR variability
%
clear all; clc

define_dirs

dataset_ix = 1;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
pupil = load_pupil(dataset_ix,time);



[offset_indices, onset_indices] = get_onsets(whisk_amp,acquisition_rate);

figure, x = whisk_amp;
plot(time,x,'k'), hold on
for k = 1:size(onset_indices,2)
    t = offset_indices(k):onset_indices(k);
    plot(time(t),x(t),'c')
    plot(time(onset_indices(k)),x(onset_indices(k)),'ok','MarkerFaceColor','c')
end
size(onset_indices,2)
%%
[N,T] = size(dFF);

num_bins = round(.2 * acquisition_rate);

num_bins_gap = 2*round(.2 * acquisition_rate);

dFF_pre_onset = zeros(N,1);

for t = 1:size(onset_indices,2)
    t_ = ((onset_indices(t)-num_bins)-1) : (onset_indices(t)-1);
    dFF_pre_onset = [dFF_pre_onset,dFF(:,t_)];
end

E = eig(cov(dFF_pre_onset')); L = sqrt(sum(E.^2));


num_reps = 1000;

L_shuff = zeros(num_reps,1);

for i = 1:num_reps
    dFF_shuff = zeros(N,1);
    
    for t = 1:size(onset_indices,2)
        t_start = randsample(offset_indices(t): ((onset_indices(t)-num_bins_gap)-1),1);
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

dataset_ix = 17;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,1);
cd = get_coding_dimension(dFF,A,QW);
A_or_QW = cd'*(dFF - mean(dFF,2));

[coeff, score] = pca(dFF');

zsp = zscore(A_or_QW);

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


%% Plot variance unexplained vs # PCs

dataset_ix = 6;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[~,whisk_set_point,~,~] = load_behav_data(dataset_ix,time);

[~, score] = pca(dFF');

[N,T] = size(dFF);
num_its = 20;

err = nan(N,num_its);
for n_PCs = 1:N

    reg = [score(:,1:n_PCs),ones(T,1)];
    
    for it_ix = 1:num_its

        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 

        b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

        mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
        err(n_PCs,it_ix) = mse / var(whisk_set_point(test_ixs));
    end
end
%% Give example of linear regression fitting to set pt

dataset_ix = 6;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[coeff, score] = pca(dFF');

T = size(dFF,2);

test_ixs = 1:floor(T * .2); 
train_ixs = setdiff(1:T,test_ixs);

reg = score(:,1:10);%(:,1:50);
reg = [reg,ones(T,1)];

b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

figure, plot(time(test_ixs),whisk_set_point(test_ixs),'k','LineWidth',1)
hold on, plot(time(test_ixs),reg(test_ixs,:)*b,'Color',[.72,.27,1],'LineWidth',1.5)
ylim([-.5,1.5])
%%

num_its = 100;
err_res = cell(17,1);
err_all = cell(17,1);

for dataset_ix = [1:6,9,17]

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,~,~] = load_behav_data(dataset_ix,time);

    [N,T] = size(dFF);
    [~, score] = pca(dFF');

    % Get residual of whisker set point, ie the component of WSP that 
    % is not well described by the first PC
    reg = [score(:,1),ones(T,1)];
    b = (reg'*reg) \ reg' * whisk_set_point;
    residual = whisk_set_point - reg*b;

    err_all{dataset_ix} = nan(N,num_its);
    err_res{dataset_ix} = nan(N,num_its);
    for n_PCs = 1:N

        reg = [score(:,1:n_PCs),ones(T,1)];
        for it_ix = 1:num_its
            train_ixs = block_shuffle_time(T,acquisition_rate);
            test_ixs = train_ixs(1:round(T * 0.2));
            train_ixs = setdiff(train_ixs,test_ixs); 

            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * residual(train_ixs);

            mse = mean((residual(test_ixs) - reg(test_ixs,:)*b).^2);
            err_res{dataset_ix}(n_PCs,it_ix) = mse / var(residual(test_ixs));
            
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_all{dataset_ix}(n_PCs,it_ix) = mse / var(whisk_set_point(test_ixs));
        end
    end
end
%%

figure, plot_error_snake(1:703,err_res{6}',[.5,.5,.5])

set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'FontSize',15)
xlabel('Number of PCs')
ylabel('Unexplained variance') 

%%

x = nanmean(err_res{5},2);
y = smoothdata(x,'gaussian',10);
figure, loglog(x,'k','LineWidth',1)
hold on, loglog(y,'r','LineWidth',1)
[a,b] = min(y)
%%

%% Figure 3G

coeff_1 = [];
coeff_2 = [];
coeff_3 = [];

for dataset_ix = [1:6,9,17]
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
end