clear all

define_dirs

load([basedir,'processed/regression_results'],'err_lasso','lambda')

lambda_wsp = nan(13,1);
for dataset_ix = 1:13
    [err_lasso_min(dataset_ix),ix] = min(mean(err_lasso{dataset_ix},2));
    lambda_wsp(dataset_ix) = lambda(ix);
end

clear err_lasso lambda
load([basedir,'processed/regression_results_speed'],'err_lasso','lambda')

lambda_speed = nan(13,1);
for dataset_ix = 1:13
    [err_lasso_min(dataset_ix),ix] = min(mean(err_lasso{dataset_ix},2));
    lambda_speed(dataset_ix) = lambda(ix);
end
clear err_lasso lambda
%% redo regression

num_its = 10;
b_wsp = cell(13,1);
b_speed = cell(13,1);
b_corr = nan(13,num_its);

for dataset_ix = 1:13   

    [dFF,time,acquisition_rate] = load_data(dataset_ix);
    [~,whisk_set_point,~,~,speed] = load_behav_data(dataset_ix,time);
    
    [N,T] = size(dFF);
    
    b_wsp{dataset_ix} = nan(N,num_its);
    b_speed{dataset_ix} = nan(N,num_its);
    
    % Random iterations
    for it_ix = 1:num_its
        disp(it_ix)
        
        train_ixs = block_shuffle_time(T,acquisition_rate);
        test_ixs = train_ixs(1:round(T * 0.2));
        train_ixs = setdiff(train_ixs,test_ixs); 
        
        reg = dFF';
        
        % WSP
        b_wsp{dataset_ix}(:,it_ix) = lasso(reg(train_ixs,:),whisk_set_point(train_ixs),'Lambda',lambda_wsp(dataset_ix));

        % Speed
        b_speed{dataset_ix}(:,it_ix) = lasso(reg(train_ixs,:),speed(train_ixs),'Lambda',lambda_speed(dataset_ix));

        b_corr(dataset_ix,it_ix) = corr(b_wsp{dataset_ix}(:,it_ix),b_speed{dataset_ix}(:,it_ix));
        
    end
end
