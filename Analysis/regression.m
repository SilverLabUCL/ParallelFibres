% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;

%% Example PC1 vs WSP
% Figure 4A,B

dataset_ix = 2;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
cd = get_coding_dimension(dFF,A,QW);
A_or_QW = cd'*(dFF - mean(dFF,2));

[~, score] = pca(dFF');

pc1 = score(:,1);
if corr(pc1,whisk_set_point) < 0
    pc1 = -pc1;
end

figure, plot(time,whisk_set_point,'k','LineWidth',1.5)
xlim([285,350])
set(gca,'FontSize',15)
xlabel('Time (s)')
ylabel('Whisker set point')

figure, plot(time,pc1,'k','LineWidth',1.5)
xlim([285,350])
set(gca,'FontSize',15)
xlabel('Time (s)')
ylabel('PC 1')

zsp = zscore(A_or_QW);

figure, hold on
for k = 1:length(pc1)
    c = (zsp(k)+1)/2;
    if c < 0
        c = 0;
    elseif c > 1
        c = 1;
    end
    c=c*[1,0,1]+(1-c)*[0,1,1];
    
    plot(pc1(k),whisk_set_point(k),'.','Color',c,'MarkerSize',10)
end
set(gca,'FontSize',15)
xlabel('PC 1')
ylabel('Whisker set point')


%% PC1 correlations with state / kinematics
% Figure 4C

[C_wsp,C_amp,C_spd,C_pupil] = deal(nan(13,1));

[C_wsp_A,C_amp_A,C_spd_A,C_pupil_A] = deal(nan(13,1));

[C_wsp_QW,C_amp_QW,C_spd_QW,C_pupil_QW] = deal(nan(13,1));

C_state = nan(13,1);

for dataset_ix = 1:13
    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    pupil = load_pupil(dataset_ix,time);

    % Get score
    [coeff, score] = pca(dFF');
    pc1 = score(:,1);
    
    % Remove NaNs
    ix = find(~isnan(whisk_set_point) & ~isnan(speed))';
    whisk_set_point = whisk_set_point(ix);
    whisk_amp = whisk_amp(ix);
    speed = speed(ix);
    pc1 = pc1(ix);
    
    % Correlation with pc1
    temp = corr(whisk_set_point,pc1);
    
    if temp < 0
        pc1 = -pc1;
        C_wsp(dataset_ix) = -temp;
    else
        C_wsp(dataset_ix) = temp;
    end
    
    
    [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
    
    C_amp(dataset_ix) = corr(whisk_amp,pc1);
    C_spd(dataset_ix) = corr(speed,pc1);
    if ~isempty(pupil)
        C_pupil(dataset_ix) = corr(pupil,pc1);
    end
    
    
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

c = [.5,.5,.5];

figure,  hold on
plot(zeros-1,C_state,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(zeros,C_wsp,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones,C_wsp_QW,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(2*ones,C_wsp_A,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)

plot(-1+.2*[-1,1],mean(C_state)*[1,1],'k','LineWidth',3)
plot(0+.2*[-1,1],mean(C_wsp)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],mean(C_wsp_QW)*[1,1],'k','LineWidth',3)
plot(2+.2*[-1,1],mean(C_wsp_A)*[1,1],'k','LineWidth',3)

set(gca,'Xtick',[-1,0,1,2])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Correlation'))
xlim([-1.5,2.5])
ylim([-.5,1])

signrank(C_state,C_wsp)
signrank(C_state,C_wsp_QW)
signrank(C_state,C_wsp_A)

%% Give example of linear regression fitting to set pt
% Figure 5A

dataset_ix = 2;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);

[coeff, score] = pca(dFF');

T = size(dFF,2);

test_ixs = (1:floor(T * .2)) ; 
train_ixs = setdiff(1:T,test_ixs);

figure, plot(time(test_ixs),whisk_set_point(test_ixs),'k','LineWidth',1)
xlabel('Time (s)'), ylabel('WSP')
set(gca,'FontSize',15)
title('WSP')
ylim([-40,50])

% Lin reg with increasing numbers of PCs
for num_PCs = [1,10,100]

    reg = score(:,1:num_PCs);
    reg = [reg,ones(T,1)];

    b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

    figure, plot(time(test_ixs),reg(test_ixs,:)*b,'Color',[.72,.27,1],'LineWidth',1.5)
    xlabel('Time (s)'), ylabel('Predicted WSP')

    mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
    title(['#PCs = ',num2str(num_PCs),', MSE = ',num2str(mse)])
    set(gca,'FontSize',15)
    ylim([-40,50])
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

figure, plot(time(test_ixs),reg(test_ixs,:)*b,'Color',[.52,.52,.52],'LineWidth',1.5)
xlabel('Time (s)'), ylabel('Predicted WSP')
title(['Best PF, MSE = ',num2str(mse(n_min))])
set(gca,'FontSize',15)
ylim([-40,50])

%% Regression - lasso, PCR, individual PFs
% takes forever

num_its = 10;
err_lasso = cell(13,1);
frac_nonzero_b = cell(13,1);
err_PFs = cell(13,1);
err_PCs = cell(13,1);
b_lasso = cell(13,1);
b_PCs = cell(13,1);

lambda = [0,logspace(-3,0,15)];

% Set to 0 if regressing against WSP
%        1 if regressing against speed
regress_speed = 0;

tic
for dataset_ix = 1:13   
    toc, tic

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,~,~,speed] = load_behav_data(dataset_ix,time);
    
    % overwrite 
    if regress_speed == 1
        % remove nans if exist
        if numel(isnan(speed))>0
            ixnan = find(isnan(speed));
            speed(ixnan)=[];
            time(ixnan)=[];
            dFF(:,ixnan)=[];
        end
        whisk_set_point = speed;
    end
    
    [N,T] = size(dFF);
    
    [~, score] = pca(dFF');
    
    err_lasso{dataset_ix} = nan(length(lambda),num_its);
    b_lasso{dataset_ix} = nan(length(lambda),num_its,N);
    frac_nonzero_b{dataset_ix} = nan(length(lambda),num_its);
    
    err_PCs{dataset_ix} = nan(N,num_its);
    b_PCs{dataset_ix} = nan(N,num_its);
    
    err_PFs{dataset_ix} = nan(N,num_its);
    
    % Random iterations
    for it_ix = 1:num_its
        disp(it_ix)
        
        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 
                
        for lambda_ix = 1:length(lambda)
    
            reg = dFF';
            
            [b,fitinfo] = lasso(reg(train_ixs,:),whisk_set_point(train_ixs),'Lambda',lambda(lambda_ix));
            mse = mean((whisk_set_point(test_ixs) - ( reg(test_ixs,:)*b+ fitinfo.Intercept ) ).^2);
            
            err_lasso{dataset_ix}(lambda_ix,it_ix) = mse / var(whisk_set_point(test_ixs));
            b_lasso{dataset_ix}(lambda_ix,it_ix,:) = b;
            frac_nonzero_b{dataset_ix}(lambda_ix,it_ix) = sum(b~=0) / N;
        end
        
        mse_PFs = nan(N,1);
        for n = 1:N
            % different numbers of PCs
            reg = [score(:,1:n),ones(T,1)];
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_PCs{dataset_ix}(n,it_ix) = mse / var(whisk_set_point(test_ixs));
            b_PCs{dataset_ix}(lambda_ix,it_ix,1:n) = b(1:end-1);
            
            % each PF
            reg = [dFF(n,:)',ones(T,1)];
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse_PFs(n) = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_PFs{dataset_ix}(n,it_ix) = mse_PFs(n) / var(whisk_set_point(test_ixs));
        end
        
    end
end

if regress_speed == 0
    filename = [basedir,'processed/regression_results_wsp'];
    disp(filename)
    save(filename,'err_PCs','err_PFs','err_lasso','frac_nonzero_b','lambda','b_PCs','b_lasso')
elseif regress_speed == 1
    filename = [basedir,'processed/regression_results_speed'];
    disp(filename)
    save(filename,'err_PCs','err_PFs','err_lasso','frac_nonzero_b','lambda','b_PCs','b_lasso')
end

%% Plot example variance unexplained vs. #PCs
% Figure 5B

load([basedir,'processed/regression_results_wsp'])
%load([basedir,'processed/regression_results_speed'])

%
dataset_ix=2;
figure, plot_error_snake(1:size(err_PCs{dataset_ix},1),err_PCs{dataset_ix}','k')

set(gca,'XScale','log')
set(gca,'FontSize',15)
xlabel('Number of PCs')
ylabel('Unexplained variance')

%% Compare 10 PCs vs optimal vs best PF
% Figure 5C

err_1PC = nan(13,1);
err_10PCs = nan(13,1);
err_best_PCs = nan(13,1);
num_best_PCs = nan(13,1);

for dataset_ix = 1:13
    err_1PC(dataset_ix) = min(mean(err_PFs{dataset_ix},2));
    err_10PCs(dataset_ix) = mean(err_PCs{dataset_ix}(10,:));
    [err_best_PCs(dataset_ix), num_best_PCs(dataset_ix)] = min(mean(err_PCs{dataset_ix},2));
end

c = [.5,.5,.5];

% Remove datasts 4 and 10 ONLY if doing analysis on speed!
% err_1PC([4,10]) = nan;
% err_10PCs([4,10]) = nan;
% err_best_PCs([4,10]) = nan;
% num_best_PCs([4,10]) = nan;

figure,  hold on
for dataset_ix = 1:13
    plot(0:2,[err_1PC(dataset_ix),err_10PCs(dataset_ix),err_best_PCs(dataset_ix)],'k','LineWidth',1)
end
plot(zeros,err_1PC,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones,err_10PCs,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(2*ones,err_best_PCs,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(0+.2*[-1,1],nanmean(err_1PC)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],nanmean(err_10PCs)*[1,1],'k','LineWidth',3)
plot(2+.2*[-1,1],nanmean(err_best_PCs)*[1,1],'k','LineWidth',3)
set(gca,'Xtick',[0,1,2])
set(gca,'XtickLabel',{'1 PC','10 PCs', 'Optimal'})
set(gca,'FontSize',15)
set(ylabel('Unexplained variance'))
xlim([-.5,2.5])
ylim([0,1])

signrank(err_best_PCs,err_10PCs)
signrank(err_best_PCs,err_1PC)

%% Plot number of PCs in optimal case
% Figure 5D

figure,  hold on
plot(zeros,num_best_PCs,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
set(gca,'FontSize',15)
set(ylabel('Optimal # PCs'))
ylim([0,200])

%% LASSO example

dataset_ix = 1;
figure, plot_error_snake(lambda,err_lasso{dataset_ix}','k')

set(gca,'FontSize',15)
xlabel('LASSO penalty')
ylabel('Unexplained variance') 

%% Find minimum and plot lasso vs lambda
% Figure 5E

err_best_PF = nan(13,1);
err_lasso_min = nan(13,1);
num_PFs_lasso = nan(13,1);

for dataset_ix = 1:13

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [N,T] = size(dFF);

    err_best_PF(dataset_ix) = min(mean(err_PFs{dataset_ix},2));
    [err_lasso_min(dataset_ix),ix] = min(mean(err_lasso{dataset_ix},2));
    num_PFs_lasso(dataset_ix) = N * frac_nonzero_b{dataset_ix}(ix);

end

% Remove datasts 4 and 10 ONLY if doing analysis on speed!
% err_best_PF([4,10]) = nan;
% err_lasso_min([4,10]) = nan;
% num_PFs_lasso([4,10]) = nan;

c = [.5,.5,.5];

figure,  hold on
for dataset_ix = 1:13
    plot(0:1,[err_best_PF(dataset_ix),err_lasso_min(dataset_ix)],'k','LineWidth',1)
end
plot(zeros,err_best_PF,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones,err_lasso_min,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(0+.2*[-1,1],nanmean(err_best_PF)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],nanmean(err_lasso_min)*[1,1],'k','LineWidth',3)
set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{'Best PF','Optimal # PFs'})
set(gca,'FontSize',15)
set(ylabel('Unexplained variance'))
xlim([-.5,1.5])
ylim([0,1])

signrank(err_best_PF,err_lasso_min)

%% Plot optimal number PFs in LASSO regression
% Figure 5F

figure,  hold on
plot(zeros,num_PFs_lasso,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
set(gca,'FontSize',15)
set(ylabel('Optimal # PFs'))
ylim([0,500])

%% Plot regression coefficients for speed vs wsp

load([basedir,'processed/regression_results_wsp'])
b_lasso_wsp = b_lasso;
err_lasso_wsp = err_lasso;

load([basedir,'processed/regression_results_speed'])
b_lasso_speed = b_lasso;
err_lasso_speed = err_lasso;

b_corr = nan(13,1);
err = nan(13,1);
for dataset_ix = 1:13
    if ~(dataset_ix == 4 || dataset_ix == 10)
        
        [err_wsp,ix] = min(mean(err_lasso_wsp{dataset_ix},2));
        b_mean_wsp = (mean(b_lasso_wsp{dataset_ix}(ix,:,:),2));
        
        [err_speed,ix] = min(mean(err_lasso_speed{dataset_ix},2));
        b_mean_speed = (mean(b_lasso_speed{dataset_ix}(ix,:,:),2));
        
        b_corr(dataset_ix) = corr(b_mean_wsp(:),b_mean_speed(:));
        err(dataset_ix) = 1/2*(err_wsp + err_speed);
    end
end

signrank(b_corr)

figure, plot(err,b_corr,'ok','MarkerFaceColor','k')
set(gca,'FontSize',15)
xlabel('Regression error')
ylabel('Correlation between coefficients')
axis([0,.8,0,.5])

%% OLD version
% Calculate variance unexplained vs # PCs 
% ONLY during QW
% Takes forever

num_its = 10;
err_PCs = cell(13,1);
err_PCs_sh = cell(13,1);

for dataset_ix = 1:13

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    [~,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
    
    % Shuffle time points
    T_sh = block_shuffle_time(size(dFF,2),acquisition_rate);
    dFF_sh = dFF(:,T_sh);
    whisk_set_point_sh = whisk_set_point(T_sh);
        
    % Only do regression on QW times
    dFF_QW = [];
    dFF_QW_sh = [];
    whisk_set_point_QW = [];
    whisk_set_point_QW_sh = [];
    for k = 1:length(QW)
        ix = QW(k,1):QW(k,2);
        dFF_QW = [dFF_QW, dFF(:,ix)];
        dFF_QW_sh = [dFF_QW_sh, dFF_sh(:,ix)];
        whisk_set_point_QW = [whisk_set_point_QW; whisk_set_point(ix)];
        whisk_set_point_QW_sh = [whisk_set_point_QW_sh; whisk_set_point_sh(ix)];
    end
    
    clear dFF dFF_sh whisk_set_point whisk_set_point_sh
    
    [N,T] = size(dFF_QW);
    [~, score] = pca(dFF_QW');
    [~, score_sh] = pca(dFF_QW_sh');
    
    err_PCs{dataset_ix} = nan(N,num_its);
    err_PCs_sh{dataset_ix} = nan(N,num_its);
    % Random iterations
    for it_ix = 1:num_its
        disp(it_ix)
        
        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 
        
        for n = 1:N

            % different numbers of PCs
            reg = [score(:,1:n),ones(T,1)];
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point_QW(train_ixs);

            mse = mean((whisk_set_point_QW(test_ixs) - reg(test_ixs,:)*b).^2);
            err_PCs{dataset_ix}(n,it_ix) = mse / var(whisk_set_point_QW(test_ixs));
            
            % Shuffled case
            reg = [score_sh(:,1:n),ones(T,1)];
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point_QW_sh(train_ixs);

            mse = mean((whisk_set_point_QW_sh(test_ixs) - reg(test_ixs,:)*b).^2);
            err_PCs_sh{dataset_ix}(n,it_ix) = mse / var(whisk_set_point_QW_sh(test_ixs));
            
        end        
        
    end
end

%save([basedir,'processed/regression_results_QW_only'],'err_PCs','err_PCs_sh')

%% State-dependent regression
% ONLY during QW
% Takes forever 

num_its = 10;
err_PCs_QW = cell(13,1);
err_PCs_sh = cell(13,1);
b_QW = cell(13,1);

% Set to 0 if regressing during QW
%        1 if regressing during AS
regress_AS = 1;

tic
for dataset_ix = [1:4,6:9,11:13]%1:13
    toc, tic

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    [AS,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
    
    if regress_AS == 1
        QW = AS;
    end
    
    % Only do regression on QW times
    ix_QW = [];
    for k = 1:length(QW)
        ix_QW = [ix_QW, QW(k,1):QW(k,2)];
    end
    
    dFF_QW = dFF(:,ix_QW);
    whisk_set_point_QW = whisk_set_point(ix_QW);
    
    [N,T] = size(dFF_QW);
    [~, score_QW] = pca(dFF_QW'); 
    [~, score] = pca(dFF');   
    
    err_PCs_QW{dataset_ix} = nan(N,num_its);
    err_PCs_sh{dataset_ix} = nan(N,num_its);
    b_QW{dataset_ix} = nan(N,num_its);
    % Random iterations
    for it_ix = 1:num_its
        disp(it_ix)
        
        % Get training and testing indices
        % Indexed in terms of ix_QW
        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 
        
        % Testing data
        whisk_set_point_QW_test = whisk_set_point_QW(test_ixs);
        score_QW_test = score_QW(test_ixs,:);
        score_test = score(ix_QW(test_ixs),:);
        
        % Training data for QW-trained decoder
        whisk_set_point_QW_train = whisk_set_point_QW(train_ixs);
        score_QW_train = score_QW(train_ixs,:);
        
        % To get available indices for training control decoder
        % First remove all test indices
        train_ixs_sh = setdiff(1:size(dFF,2),ix_QW(test_ixs));
        % Then shuffle and take the same number as train_ixs
        ix_sh = block_shuffle_time(numel(train_ixs_sh),acquisition_rate);
        train_ixs_sh = train_ixs_sh(ix_sh);
        % Finally take the same number as train_ixs
        train_ixs_sh = train_ixs_sh(1:numel(train_ixs));

        % Training data for QW-trained decoder
        whisk_set_point_train = whisk_set_point(train_ixs_sh);
        score_train = score(train_ixs_sh,:);
        
        for n = 1:N

            % different numbers of PCs
            reg_train = [score_QW_train(:,1:n),ones(numel(train_ixs),1)];
            reg_test = [score_QW_test(:,1:n),ones(numel(test_ixs),1)];
            
            b = (reg_train'*reg_train) \ reg_train' * whisk_set_point_QW_train;

            mse = mean((whisk_set_point_QW_test - reg_test*b).^2);
            err_PCs_QW{dataset_ix}(n,it_ix) = mse / var(whisk_set_point_QW_test);
            
            b_QW{dataset_ix}(1:n,it_ix) = b(1:end-1);
            
            % Shuffled case
            reg_train = [score_train(:,1:n),ones(numel(train_ixs),1)];
            reg_test = [score_test(:,1:n),ones(numel(test_ixs),1)];
            
            b = (reg_train'*reg_train) \ reg_train' * whisk_set_point_train;

            mse = mean((whisk_set_point_QW_test - reg_test*b).^2);
            err_PCs_sh{dataset_ix}(n,it_ix) = mse / var(whisk_set_point_QW_test);
            
        end        
        
    end
end

if regress_AS == 0
    filename = [basedir,'processed/regression_results_state_dependent_decoding_QW'];
    disp(filename);
    save(filename,'err_PCs_QW','err_PCs_sh','b_QW')
elseif regress_AS == 1
    filename = [basedir,'processed/regression_results_state_dependent_decoding_AS'];
    disp(filename);
    err_PCs_AS=err_PCs_QW; b_AS = b_QW;
    save(filename,'err_PCs_AS','err_PCs_sh','b_AS')
end

%% Find minimum and plot vs shuffled

err_best_PCs = nan(13,1);
err_best_PCs_sh = nan(13,1);

for dataset_ix = 1:13

    err_best_PCs(dataset_ix) = min(mean(err_PCs_QW{dataset_ix},2));
    err_best_PCs_sh(dataset_ix) = min(mean(err_PCs_sh{dataset_ix},2));
    
end

c = [.5,.5,.5];

figure,  hold on
for dataset_ix = 1:13
    plot(0:1,[err_best_PCs(dataset_ix),err_best_PCs_sh(dataset_ix)],'k','LineWidth',1)
end
plot(zeros,err_best_PCs,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones,err_best_PCs_sh,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(0+.2*[-1,1],nanmean(err_best_PCs)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],nanmean(err_best_PCs_sh)*[1,1],'k','LineWidth',3)
set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{'QW only','Shuffled times'})
set(gca,'FontSize',15)
set(ylabel('Unexplained variance'))
xlim([-.5,1.5])
ylim([0,1])

signrank(err_best_PCs,err_best_PCs_sh)