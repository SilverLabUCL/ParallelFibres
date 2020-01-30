% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;

%% Example
% Figure 3A

dataset_ix = 2;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
cd = get_coding_dimension(dFF,A,QW);
A_or_QW = cd'*(dFF - mean(dFF,2));

[~, score] = pca(dFF');

figure, plot(time,whisk_set_point,'k','LineWidth',1.5)
xlim([285,350])
set(gca,'FontSize',15)
xlabel('Time (s)')
ylabel('Whisker set point')

figure, plot(time,score(:,1),'k','LineWidth',1.5)
xlim([285,350])
set(gca,'FontSize',15)
xlabel('Time (s)')
ylabel('PC 1')
%%
zsp = zscore(A_or_QW);

figure, hold on
for k = 1:length(score)
    c = (zsp(k)+1)/2;
    if c < 0
        c = 0;
    elseif c > 1
        c = 1;
    end
    c=c*[1,0,1]+(1-c)*[0,1,1];
    
    plot(score(k,1),whisk_set_point,'.','Color',c)
end
set(gca,'FontSize',15)
xlabel('PC 1')
ylabel('Whisker set point')


%% Figure 3C

[C_wsp,C_amp,C_spd,C_pupil] = deal(nan(15,1));

[C_wsp_A,C_amp_A,C_spd_A,C_pupil_A] = deal(nan(15,1));

[C_wsp_Q,C_amp_Q,C_spd_Q,C_pupil_Q] = deal(nan(15,1));

C_state = nan(17,1);

for dataset_ix = 1:15
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
    
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,1);
    %buffer = round(acquisition_rate * .5);
    %[A,QW] = get_A_QW_periods(pc1,buffer,1);
    pause

    C_amp(dataset_ix) = corr(whisk_amp,pc1);
    C_spd(dataset_ix) = corr(speed,pc1);
    if ~isempty(pupil)
        C_pupil(dataset_ix) = corr(pupil,pc1);
    end
    
    
    [pc1_A,whisk_set_point_A,whisk_amp_A,pupil_A,speed_A] = deal([]);
    ix_A = [];
    for k = 1:length(A)
        ix_A = [ix_A,A(k,1):A(k,2)];
    end
    pc1_A = pc1(ix_A);
    whisk_set_point_A = whisk_set_point(ix_A);
    whisk_amp_A = whisk_amp(ix_A);
    speed_A = speed(ix_A);
    if ~isempty(pupil)
        pupil_A = pupil(ix_A);
    end
        
    C_wsp_A(dataset_ix) = corr(whisk_set_point_A,pc1_A);
    C_amp_A(dataset_ix) = corr(whisk_amp_A,pc1_A);
    C_spd_A(dataset_ix) = corr(speed_A,pc1_A);
    if ~isempty(pupil)
        C_pupil_A(dataset_ix) = corr(pupil_A,pc1_A);
    end
    
    [pc1_QW,whisk_set_point_QW,whisk_amp_QW,pupil_QW,speed_QW] = deal([]);
    
    ix_QW = [];
    for k = 1:length(QW)
        ix_QW = [ix_QW,QW(k,1):QW(k,2)];
    end
    pc1_QW = pc1(ix_QW);
    whisk_set_point_QW = whisk_set_point(ix_QW);
    whisk_amp_QW = whisk_amp(ix_QW);
    speed_QW = speed(ix_QW);
    if ~isempty(pupil)
        pupil_QW = pupil(ix_QW);
    end
        
    C_wsp_QW(dataset_ix) = corr(whisk_set_point_QW,pc1_QW);
    C_amp_QW(dataset_ix) = corr(whisk_amp_QW,pc1_QW);
    C_spd_QW(dataset_ix) = corr(speed_QW,pc1_QW);
    if ~isempty(pupil)
        C_pupil_QW(dataset_ix) = corr(pupil_QW,pc1_QW);
    end
    
    C_state(dataset_ix) = corr([pc1_A;pc1_QW],[ones(size(pc1_A));zeros(size(pc1_QW))]);
end

figure,  hold on
plot(zeros-1,C_state,'ok','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(zeros,C_wsp,'ob','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(ones,C_wsp_QW,'oc','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
plot(2*ones,C_wsp_A,'om','MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',12)
set(gca,'Xtick',[-1,0,1,2])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Cross correlation'))
xlim([-1.5,2.5])
ylim([-.5,1])

signrank(C_wsp,C_wsp_QW)
signrank(C_wsp_QW,C_wsp_A)
signrank(C_wsp,C_wsp_A)
signrank(C_state,C_wsp)


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

dataset_ix = 15;

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

dataset_ix = 2; %6

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[coeff, score] = pca(dFF');

T = size(dFF,2);

test_ixs = (1:floor(T * .2)) ; 
train_ixs = setdiff(1:T,test_ixs);

for num_PCs = [1,10,100]

    reg = score(:,1:num_PCs);
    reg = [reg,ones(T,1)];

    b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

    figure, plot(time(test_ixs),whisk_set_point(test_ixs),'k','LineWidth',1)
    hold on, plot(time(test_ixs),reg(test_ixs,:)*b,'Color',[.72,.27,1],'LineWidth',1.5)
    ylim([-.5,1.5])
    mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
    title(mse)
end

% now compare to best fitting 
mse = zeros(size(dFF,1),1);
for n = 1:size(dFF,1)

    reg = dFF(n,:)';
    reg = [reg,ones(T,1)];

    b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);
    mse(n) = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
end

[~,n_min] = min(mse);
reg = dFF(n_min,:)';
reg = [reg,ones(T,1)];

b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

figure, plot(time(test_ixs),whisk_set_point(test_ixs),'k','LineWidth',1)
hold on, plot(time(test_ixs),reg(test_ixs,:)*b,'Color',[.52,.52,.52],'LineWidth',1.5)
ylim([-.5,1.5])
title(mse(n_min))

%%

num_its = 3;
err_res = cell(17,1);
err_all = cell(17,1);

for dataset_ix = [1:6,9,17]

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,~,~] = load_behav_data(dataset_ix,time);

    [N,T] = size(dFF);
    [~, score] = pca(dFF');

    % Get residual of whisker set point, ie the component of WSP that 
    % is not well described by the first PC
    %reg = [score(:,1),ones(T,1)];
    %b = (reg'*reg) \ reg' * whisk_set_point;
    %residual = whisk_set_point - reg*b;

    err_all{dataset_ix} = nan(N,num_its);
    err_res{dataset_ix} = nan(N,num_its);
    for it_ix = 1:num_its
        
        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 
        
        for n_PCs = 1:N

            reg = [score(:,1:n_PCs),ones(T,1)];
        
            %b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * residual(train_ixs);

            %mse = mean((residual(test_ixs) - reg(test_ixs,:)*b).^2);
            %err_res{dataset_ix}(n_PCs,it_ix) = mse / var(residual(test_ixs));
            
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_all{dataset_ix}(n_PCs,it_ix) = mse / var(whisk_set_point(test_ixs));
        end
    end
end

%% Randomly choose QW period

num_its = 10;
err_res = cell(17,1);
err_all = cell(17,1);

for dataset_ix = [1:6,9,17]

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,0);
    % Only keep QW periods of at least 500 ms
    
    
    [N,T] = size(dFF);
    [~, score] = pca(dFF');

    % Get residual of whisker set point, ie the component of WSP that 
    % is not well described by the first PC
    reg = [score(:,1),ones(T,1)];
    b = (reg'*reg) \ reg' * whisk_set_point;
    residual = whisk_set_point - reg*b;

    err_all{dataset_ix} = nan(N,num_its);
    err_res{dataset_ix} = nan(N,num_its);
    
    for it_ix = 1:num_its
            
        QW_resort = randsample(size(QW,1),size(QW,1));
        QW = QW(QW_resort,:);
        QW_ixes = [];
        for k = 1:size(QW,1)
            QW_ixes = [QW_ixes, QW(k,1):QW(k,2)];
        end

        test_ixs = QW_ixes(1:round(T * 0.1));
        train_ixs = setdiff(1:T,test_ixs); 

        for n_PCs = 1:N
            reg = [score(:,1:n_PCs),ones(T,1)];

            %b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * residual(train_ixs);

            %mse = mean((residual(test_ixs) - reg(test_ixs,:)*b).^2);
            %err_res{dataset_ix}(n_PCs,it_ix) = mse / var(residual(test_ixs));
            
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_all{dataset_ix}(n_PCs,it_ix) = mse / var(whisk_set_point(test_ixs));
        end
    end
end

%%
dataset_ix=1;
figure, plot(1:size(err_all{dataset_ix},1),err_all{dataset_ix}')
%plot_error_snake(1:size(err_all{dataset_ix},1),err_all{dataset_ix}',[.5,.5,.5])
[m,ix] = min((err_all{dataset_ix}));
mean(ix)

set(gca,'YScale','log')
set(gca,'XScale','log')
set(gca,'FontSize',15)
xlabel('Number of PCs')
ylabel('Unexplained variance') 
%%
opt_num_PCs = nan(17,1);
for dataset_ix = [1:6,9,17]
    [m,ix] = min((err_all{dataset_ix}));
    opt_num_PCs(dataset_ix) = mean(ix);
end
figure, plot(zeros(17,1),opt_num_PCs,'ok','MarkerSize',8,'Linewidth',1.5)
set(gca,'Box','off','XTick',[],'FontSize',15,'Yscale','log')
ylabel('Optimal number PCs')
ylim([1,500])
%%

x = nanmean(err_res{5},2);
y = smoothdata(x,'gaussian',10);
figure, loglog(x,'k','LineWidth',1)
hold on, loglog(y,'r','LineWidth',1)
[a,b] = min(y)
%%