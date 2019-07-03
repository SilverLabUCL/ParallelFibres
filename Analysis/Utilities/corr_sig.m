% Correlation with significance test (two tailed)

function [C,p] = corr_sig(x,y,acquisition_rate)

    % Convert to columns
    if isrow(x)
        x = x';
    end
    if size(x,1) ~= size(y,1)
        y = y';
    end

    % Remove nans
    ix = find(~isnan(x+sum(y,2)));
    x = x(ix); y = y(ix,:);
    
    T = size(x,1);
    
    % True correlation
    C = corr(x,y);

    % Get shuffle distribution (null distribution)
    num_reps = 1000;
    C_shuff = nan(num_reps,size(y,2));
    for rep = 1:num_reps
        T_shuffled = block_shuffle_time(T,acquisition_rate);
        C_shuff(rep,:) = corr(x,y(T_shuffled,:));
    end
    
    % Calculate two-tailed p value
    p = sum(C_shuff > abs(C) |  C_shuff < -abs(C))/num_reps;