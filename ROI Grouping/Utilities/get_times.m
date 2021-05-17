%  This function calculates the timestamps (all time in ms)
%
%  Inputs:
%    offset_initial         Time it takes to reach that pixel in the cycle
%                                   (usually, in MatrixTime_patch)
%    Cycle_time             Time of a full cycle (usually,
%                                   MatrixTime(end,end)
%    Numb_cyles             Number of cycles
%    Time_between_trial     Vector of intertrial times (ignore last element)
%    Numb_trials            Number of trials in the experiment
%
%  Ouput:
%    time                   Vector of times

function time = get_times(offset_initial,Cycle_time,Numb_cycles,Time_between_trial,Numb_trials)

    % Time it takes to finish the cycle from selected pixel to the last pixel
    end_of_cycle = Cycle_time - offset_initial;

    % First trial: add initial offset to cycle time
    time = nan(1,Numb_trials*Numb_cycles);
    time(1:Numb_cycles) = offset_initial + (0:(Numb_cycles-1)) * Cycle_time;
    
    % For subsequent trials
    for tr = 2:Numb_trials
        trial_start = Numb_cycles*(tr-1)+1;
        trial_end = Numb_cycles*tr;
        % Add (1) the time to finish the last cycle of the previous trial,
        % (2) the intertrial interval
        offset = end_of_cycle +  Time_between_trial(tr-1);
        time(trial_start:trial_end) = time(1:Numb_cycles) + time(trial_start-1) + offset;
    end
    