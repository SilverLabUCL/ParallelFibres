
function [aroused,quiescent] = find_active_periods(time,speed,varargin)
    

    if ~isempty(varargin)
        MI_threshold = varargin{1};
    else
        MI_threshold = 0.2;
    end
%%
    MI_baseline = median(speed(:));
    MI_std = std(speed(:));
    start_speed_threshold = .003;
    end_speed_threshold = .002;
    end_speed_time = 150; %ms
    end_speed_bins = round(end_speed_time/mean(mean(diff(time))));


    isDone = 0;

    % First find all whisking periods
    aroused = [];
    x = speed(:);
    v = diff(x); % First derivative -- velocity
%%
    while isDone == 0


        % Beginning of whisking period defined by rapid whisker movement
        whisk_start = find(v > start_speed_threshold | x(1:end-1) > MI_threshold);

        if numel(whisk_start) == 0
            isDone = 1; % no whi,sking periods found            
        else
            whisk_start = whisk_start(1);
            % End of whisking period defined by slow movement
            % Find (1) whisker speed slower than threshold, (2) whisker retracting,
            % (3) after whisking start period, and (4) whisker position is near
            % baseline. This gives potential starting pts of whisker slowing
            slowing_start = find(v>-end_speed_threshold & x(1:end-1)<MI_baseline+0.1*MI_std & v<0 & (1:length(v))'>=whisk_start);

            if numel(slowing_start) == 0 
                whisk_end = size(time,1);
                if max(x(whisk_start:whisk_end))>MI_threshold
                    aroused = [aroused; time(whisk_start), time(whisk_end)];
                end
                isDone = 1;
            else
                whisk_start = whisk_start(1);
                % Now look for when absolute whisker speed is slower than threshold
                % over potential whisker slowing pts
                slowing = (abs(v)<end_speed_threshold & (1:length(v))'>=slowing_start(1));
                % Look for when whisker speed has been slow for sustained period of
                % time (ie end_speed_time)
                temp = slowing(1:(end-end_speed_bins+1));
                for k = 2:end_speed_bins
                    temp = temp.*slowing(k:(end-end_speed_bins+k));
                end
                % End of whisking period
                whisk_end = find(temp);
                whisk_end = whisk_end(1);

                % Only add to whisking period if peak is above threshold
                if max(x(whisk_start:whisk_end))>MI_threshold
                    aroused = [aroused; time(whisk_start), time(whisk_end)];
                end

                % Remove all activity from before that period
                x(1:whisk_end) = 0;
                v(1:whisk_end) = 0;

                if (v > start_speed_threshold) == 0
                    isDone = 1;
                end
            end
        end

    end

    quiescent = [];
    if numel(aroused)==0
        quiescent = [time(1),time(end)];
    else            
        if aroused(1,1) > time(1)
            quiescent = [quiescent; time(1), aroused(1,1)];
        end
        for k = 1:size(aroused,1)-1
            quiescent = [quiescent; aroused(k,2), aroused(k+1,1)];
        end 
        if aroused(end,2) < time(end)
            quiescent = [quiescent; aroused(end,2), time(end)];
        end
    end


%% test
%figure, plot(time,speed,'k')
