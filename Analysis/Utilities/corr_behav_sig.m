% Correlation with significance test (two tailed)
% Inputs must be columns

function [C,p,C_shuff] = corr_behav_sig(x,y,acquisition_rate)

    num_reps = 1000;

    T = size(x,1);
    
    ix = find(~isnan(x) & ~isnan(y))';
    x = x(ix); y = y(ix);
    
    % True correlation
    C = corr(x,y);

    % Get shuffle distribution (null distribution)
    C_shuff = nan(num_reps,1);
    for rep = 1:num_reps
        T_shuffled = block_shuffle_time(T,acquisition_rate);
        C_shuff(rep) = corr(x,y(T_shuffled));
    end
    
    % Calculate two-tailed p value
    p = sum(C_shuff > abs(C) |  C_shuff < -abs(C))/num_reps;
    
    % Return mean of shuffle distribution
    C_shuff = mean(C_shuff);

    
    
