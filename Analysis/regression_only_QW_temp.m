
dataset_ix = 2;

[dFF,time,acquisition_rate] = load_data(dataset_ix);
[whisk_angle,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
[~,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
[~,T] = size(dFF);
T_shuff = block_shuffle_time(T,acquisition_rate);
dFF = dFF(:,T_shuff);
whisk_set_point = whisk_set_point(T_shuff);

dFF_QW = [];
whisk_set_point_QW = [];
for k = 1:length(QW)
    ix = QW(k,1):QW(k,2);
    dFF_QW = [dFF_QW, dFF(:,ix)];
    whisk_set_point_QW = [whisk_set_point_QW; whisk_set_point(ix)];
end
dFF = dFF_QW;
whisk_set_point = whisk_set_point_QW;

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


%% Calculate variance unexplained vs # PCs 
% ONLY during QW
% Takes forever

num_its = 10;
err_PCs = cell(15,1);

for dataset_ix = 1:15

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    [~,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);
    
    % Only do regression on QW times
    dFF_QW = [];
    whisk_set_point_QW = [];
    for k = 1:length(QW)
        ix = QW(k,1):QW(k,2);
        dFF_QW = [dFF_QW, dFF(:,ix)];
        whisk_set_point_QW = [whisk_set_point_QW; whisk_set_point(ix)];
    end
    dFF = dFF_QW;
    whisk_set_point = whisk_set_point_QW;

    [N,T] = size(dFF);
    [~, score] = pca(dFF');
    
    err_PCs{dataset_ix} = nan(N,num_its);
    % Random iterations
    for it_ix = 1:num_its
        
        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 
        
        for n = 1:N

            % different numbers of PCs
            reg = [score(:,1:n),ones(T,1)];
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_PCs{dataset_ix}(n,it_ix) = mse / var(whisk_set_point(test_ixs));
            
        end        
        
    end
end


%% Calculate variance unexplained vs # PCs 
% SHUFFLE to compare to ONLY QW
% Takes forever

num_its = 10;
err_PCs_shuff = cell(15,1);

for dataset_ix = 1:15

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,whisk_amp,speed] = load_behav_data(dataset_ix,time);
    [~,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate);

    [N,T] = size(dFF);
    
    % SHUFFLE DATA
    T_shuff = block_shuffle_time(T,acquisition_rate);
    dFF = dFF(:,T_shuff);
    whisk_set_point = whisk_set_point(T_shuff);

    % Only do regression on QW times
    dFF_QW = [];
    whisk_set_point_QW = [];
    for k = 1:length(QW)
        ix = QW(k,1):QW(k,2);
        dFF_QW = [dFF_QW, dFF(:,ix)];
        whisk_set_point_QW = [whisk_set_point_QW; whisk_set_point(ix)];
    end
    dFF = dFF_QW;
    whisk_set_point = whisk_set_point_QW;

    [N,T] = size(dFF);
    [~, score] = pca(dFF');
    
    err_PCs_shuff{dataset_ix} = nan(N,num_its);
    % Random iterations
    for it_ix = 1:num_its
        
        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 
        
        for n = 1:N

            % different numbers of PCs
            reg = [score(:,1:n),ones(T,1)];
            
            b = (reg(train_ixs,:)'*reg(train_ixs,:)) \ reg(train_ixs,:)' * whisk_set_point(train_ixs);

            mse = mean((whisk_set_point(test_ixs) - reg(test_ixs,:)*b).^2);
            err_PCs_shuff{dataset_ix}(n,it_ix) = mse / var(whisk_set_point(test_ixs));
            
        end
        
    end
end

%% Plot example variance unexplained vs. #PCs
% Figure 3E

dataset_ix=11;
figure, plot_error_snake(1:size(err_PCs{dataset_ix},1),err_PCs{dataset_ix}','k')
hold on, plot_error_snake(1:size(err_PCs_shuff{dataset_ix},1),err_PCs_shuff{dataset_ix}','r')

set(gca,'XScale','log')
ylim([0,1])
set(gca,'FontSize',15)
xlabel('Number of PCs')
ylabel('Unexplained variance') 

%% Compare 10 PCs vs optimal vs best PF
% Figure 3F
err_PC_min = nan(15,1);
err_shuff_min = nan(15,1);

for dataset_ix = 1:15
    err_PC_min(dataset_ix) = min(mean(err_PCs{dataset_ix},2));
    err_shuff_min(dataset_ix) = min(mean(err_PCs_shuff{dataset_ix},2));
end

figure, hold on
for dataset_ix = 1:15
	plot([0,1], [err_PC_min,err_shuff_min],'o-','Color',[.7,.7,.7],'LineWidth',2)
end

plot([0,1], [nanmean(err_PC_min),nanmean(err_shuff_min)],'s','Color','k','MarkerSize',10,'MarkerFaceColor','k')

xlim([-.5,1.5])
set(gca,'FontSize',20,'XTick',[0,1],'XTickLabel',{'QW','Shuff'})
ylabel('Unexplained Variance')


