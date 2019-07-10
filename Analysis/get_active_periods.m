function get_active_periods(pc1,time)
    
    %% Remove inactive state + transitions and recalculate
    active_signal = pc1.*(pc1>0);
    
    isDone = 0;
    active_end = 0;
    
    active = [];
    while isDone == 0
        % Beginning of whisking period defined by rapid whisker movement
        active_start = find(time>=active_end & active_signal>0);

        if numel(active_start) == 0
            isDone = 1; % no whisking periods found      
        else
            active_start = time(active_start(1));

            active_end = find(time>active_start & active_signal==0);

            if numel(active_end) == 0
                active = [active;active_start,time(end)];
                isDone = 1;
            else
                active_end = time(active_end(1));
                active = [active;active_start,active_end];
            end
        end
    end