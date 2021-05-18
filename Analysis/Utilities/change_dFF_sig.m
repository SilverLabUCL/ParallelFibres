
% This function determines whether dFF is higher for A or QW state
%
% with two tailed significance test
% Input:
%    dFF                  Population activity
%    A                    Active state periods
%    QW                   Quiet wakefulness periods
%    acquisition_rate     Acquisition rate
% 
% Output:
%    change_dFF           Vector of difference between dFF during A vs QW
%                           for each ROI
%    p                    Vector of p-values for each ROI


function [change_dFF,p] = change_dFF_sig(dFF,A,QW,acquisition_rate)

    [N,T] = size(dFF);

    % Get indices during A and QW states
    ix_A = [];
    for k = 1:length(A)
        ix_A = [ix_A,A(k,1):A(k,2)];
    end
    
    ix_QW = [];
    for k = 1:length(QW)
        ix_QW = [ix_QW,QW(k,1):QW(k,2)];
    end
    
    change_dFF = nanmean(dFF(:,ix_A),2)-nanmean(dFF(:,ix_QW),2);

    % Get shuffle distribution (null distribution)
    num_reps = 1000;
    change_dFF_shuff = nan(num_reps,N);
    for rep = 1:num_reps
        T_shuffled = block_shuffle_time(T,acquisition_rate);
        dFF_shuff = dFF(:,T_shuffled);
        change_dFF_shuff(rep,:) = nanmean(dFF_shuff(:,ix_A),2)-nanmean(dFF_shuff(:,ix_QW),2);
    end
    
    % Calculate two-tailed p value 
    p = nan(N,1);
    for k = 1:N
        p(k) = sum(change_dFF_shuff(:,k) > abs(change_dFF(k)) |  change_dFF_shuff(:,k) < -abs(change_dFF(k)))/num_reps;
    end