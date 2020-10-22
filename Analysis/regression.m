% This script plots the correlations of individual PFs with behaviour to
% show bidirectionality of responses

clear all; clc

define_dirs;

%% Example PC1 vs WSP
% Figure 3A,B

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
% Figure 3C

[C_wsp,C_amp,C_spd,C_pupil] = deal(nan(15,1));

[C_wsp_A,C_amp_A,C_spd_A,C_pupil_A] = deal(nan(15,1));

[C_wsp_QW,C_amp_QW,C_spd_QW,C_pupil_QW] = deal(nan(15,1));

C_state = nan(15,1);

for dataset_ix = 1:15
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
% Figure 3D

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

%% try lasso regression for PFs
% takes forever

num_its = 1;%0;
err_lasso = cell(15,1);
frac_nonzero_b = cell(15,1);
lambda = [0,logspace(-3,0,15)];
err_PFs = cell(15,1);

for dataset_ix = 1:15

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,~,~] = load_behav_data(dataset_ix,time);

    [N,T] = size(dFF);
    
    err_lasso{dataset_ix} = nan(length(lambda),num_its);
    % Random iterations
    for it_ix = 1:num_its
        
        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 
                
        for lambda_ix = 1:length(lambda)
    
            reg = dFF';
            
            [b,fitinfo] = lasso(reg(train_ixs,:),whisk_set_point(train_ixs),'Lambda',lambda(lambda_ix));
            mse = mean((whisk_set_point(test_ixs) - ( reg(test_ixs,:)*b+ fitinfo.Intercept ) ).^2);
            
            err_lasso{dataset_ix}(lambda_ix,it_ix) = mse / var(whisk_set_point(test_ixs));
            frac_nonzero_b{dataset_ix}(lambda_ix,it_ix) = sum(b~=0) / N;
        end
        
        mse_PFs = nan(N,1);
        for n = 1:N
            % each PF
            reg = [dFF(n,:)',ones(T,1)];
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse_PFs(n) = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_PFs{dataset_ix}(n,it_ix) = mse_PFs(n) / var(whisk_set_point(test_ixs));
        end
    end
end
%% plot lasso

figure, hold on

err_PF_min = nan(15,1);
err_lasso_min = nan(15,1);
num_PFs_lasso = nan(15,1);

for dataset_ix = 1:15

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[N,T] = size(dFF);

err_PF_min(dataset_ix) = min(err_PFs{dataset_ix});
[err_lasso_min(dataset_ix),ix] = min(err_lasso{dataset_ix});
num_PFs_lasso(dataset_ix) = N * frac_nonzero_b{dataset_ix}(ix);
% figure, plot(lambda,err_lasso{dataset_ix},'o-r','LineWidth',2,'MarkerFaceColor','r')
% hold on, plot([0,1], err_PF_min*[1,1],':k','LineWidth',2)
% [~,ix] = min(err_lasso{dataset_ix});
% title(['Optimal number of PFs: ',num2str(N * frac_nonzero_b{dataset_ix}(ix))])
% set(gca,'FontSize',20)
% xlabel('lambda')
% ylabel('Unexplained variance')
end
%%

num_PFs_lasso(8)=nan;
err_PF_min(8)=nan;
err_lasso_min(8)=nan;

figure, hold on
for dataset_ix = 1:15
	plot([0,1], [err_PF_min,err_lasso_min],'o-','Color',[.7,.7,.7],'LineWidth',2)
end
xlim([-.5,1.5])
set(gca,'FontSize',20,'XTick',[0,1],'XTickLabel',{'Best PF','LASSO'})
ylabel('Unexplained Variance')

figure, plot(num_PFs_lasso*0,num_PFs_lasso,'o','Color',[.7,.7,.7],'LineWidth',2)
ylabel('Optimal #PFs')
set(gca,'FontSize',20,'XTick',0,'XTickLabel',{''})

%% Calculate variance unexplained vs # PCs 
% Takes forever

num_its = 10;
err_PFs = cell(15,1);
err_PCs = cell(15,1);

for dataset_ix = 1:15

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,~,~] = load_behav_data(dataset_ix,time);

    [N,T] = size(dFF);
    [~, score] = pca(dFF');
    
    err_PCs{dataset_ix} = nan(N,num_its);
    err_best_PF{dataset_ix} = nan(N,num_its);
    err_PFs{dataset_ix} = nan(N,num_its);
    % Random iterations
    for it_ix = 1:num_its
        
        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 
        
        mse_PFs = nan(N,1);
        for n = 1:N

            % different numbers of PCs
            reg = [score(:,1:n),ones(T,1)];
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_PCs{dataset_ix}(n,it_ix) = mse / var(whisk_set_point(test_ixs));
            
            % each PF
            reg = [dFF(n,:)',ones(T,1)];
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse_PFs(n) = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_PFs{dataset_ix}(n,it_ix) = mse_PFs(n) / var(whisk_set_point(test_ixs));
        end
        
        
    end
end

%% Plot example variance unexplained vs. #PCs
% Figure 3E

dataset_ix=5;
figure, plot_error_snake(1:size(err_PCs{dataset_ix},1),err_PCs{dataset_ix}','k')

set(gca,'XScale','log')
set(gca,'FontSize',15)
xlabel('Number of PCs')
ylabel('Unexplained variance') 
%% Compare 10 PCs vs optimal vs best PF
% Figure 3F
err_best_PF = nan(15,1);
err_10PCs = nan(15,1);
err_best_PCs = nan(15,1);
num_best_PCs = nan(15,1);

for dataset_ix = 1:15
    err_best_PF(dataset_ix) = min(mean(err_PFs{dataset_ix},2));
    err_10PCs(dataset_ix) = mean(err_PCs{dataset_ix}(10,:));
    [err_best_PCs(dataset_ix), num_best_PCs(dataset_ix)] = min(mean(err_PCs{dataset_ix},2));
end

c = [.5,.5,.5];

figure,  hold on
plot(zeros,err_best_PF,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(ones,err_10Pcs,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(2*ones,err_best_PCs,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
plot(0+.2*[-1,1],nanmean(err_best_PF)*[1,1],'k','LineWidth',3)
plot(1+.2*[-1,1],nanmean(err_10PCs)*[1,1],'k','LineWidth',3)
plot(2+.2*[-1,1],nanmean(err_best_PCs)*[1,1],'k','LineWidth',3)
set(gca,'Xtick',[0,1])
set(gca,'XtickLabel',{})
set(gca,'FontSize',15)
set(ylabel('Unexplained variance'))
xlim([-.5,2.5])
ylim([0,1])

signrank(err_best_PCs,err_10PCs)
signrank(err_best_PCs,err_best_PF)

% Plot number of PCs in optimal case
figure,  hold on
plot(zeros,num_best_PCs,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
set(gca,'FontSize',15)
set(ylabel('Best # PCs'))
ylim([0,200])

figure,  hold on
plot(zeros,num_best_PCs,'o','MarkerFaceColor','w','Color',c,'MarkerSize',8)
set(gca,'FontSize',15)
set(ylabel('Best # PCs'))
ylim([600,800])
