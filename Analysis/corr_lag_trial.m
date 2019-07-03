% Get time lagged correlations with trial structure

function [C,lags] = corr_lag_trial(x,y,time)

    % Convert to columns
    if isrow(x)
        x = x';
    end
    if size(x,1) ~= size(y,1)
        y = y';
    end
    if isrow(time)
        time = time';
    end
    
    T = size(x,1);
    acq_rate_est = length(time) / (time(end)-time(1));
    
    trial_change_ix = find(diff(time) > (mode(diff(time))+2*std(diff(time))));
    trial_length = trial_change_ix(1);
    num_trials = length(trial_change_ix)+1;
    
    % True correlation
    C = nan(2*trial_length-1,num_trials);
    for trial = 1:num_trials
        if trial == 1
            trial_start = 1;
        else
            trial_start = trial_change_ix(trial-1)+1;
        end
        if trial == num_trials
            trial_end = T;
        else
            trial_end = trial_change_ix(trial);
        end
        x_ = x(trial_start:trial_end);
        y_ = y(trial_start:trial_end);
        [C(:,trial),lags] = xcov(x_,y_,'coeff');
        if sum(isnan(x_+y_)) > 0
            C(:,trial) = [];
        end
    end
    
    lags = lags / acq_rate_est;
    