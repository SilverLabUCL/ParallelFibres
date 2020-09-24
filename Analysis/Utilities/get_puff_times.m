%  This function calculates the timestamps (all time in ms)
%  Inputs:
%    dataset_ix         Dataset_number
%    plot_me            Boolean, if you want to plot dFF during puff (1) or not (0)
%
%  Ouput:
%    time_ix                   Vector of time indices 

function [time_ix,time_ix_shuff] = get_puff_times(dataset_ix,plot_me)
    
    if nargin < 2 || isempty(plot_me)
        plot_me = 0;
    end
    
    define_dirs;
    fname = datasets{dataset_ix};
    load([basedir,fname,'/',fname,'_time.mat'],'baseline')
    load([basedir,fname,'/',fname,'.mat'],'acquisition_rate','Numb_trials');
    load([basedir,fname,'/',fname,'_time.mat'],'Numb_cycle');

    if plot_me
        [dFF,~,~] = load_data(dataset_ix);
        figure(1)
        figure(2)
    end
    
    time_ix = [];
    time_ix_shuff = [];
    for tr = 1:Numb_trials
        trial_start = Numb_cycle*(tr-1)+1;
        trial_end = Numb_cycle*tr;
        
        num_pre = round(0 * acquisition_rate);
        num_post = round(.3 * acquisition_rate);
        
        puff_ix = round(trial_start + baseline / 1000 * acquisition_rate);
        window_start = puff_ix - num_pre;
        window_end = puff_ix + num_post;
        
        time_ix = [time_ix, window_start:window_end];
        
        shuffstart_ix = randsample((trial_start+num_pre):(trial_end-num_post),1);
        window_start_shuff = shuffstart_ix - num_pre;
        window_end_shuff = shuffstart_ix + num_post;
        time_ix_shuff = [time_ix_shuff, window_start_shuff:window_end_shuff];
        
        if plot_me
            figure(1), hold on, plot(mean(dFF(:,window_start:window_end),1),'LineWidth',2)
            figure(2), hold on, plot(mean(dFF(:,window_start_shuff:window_end_shuff),1),'k','LineWidth',2)
        end

    end
    