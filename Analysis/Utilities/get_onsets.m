


function onset_indices = get_onsets(speed,acquisition_rate,plot_me)

    if nargin < 3 || isempty(plot_me)
        plot_me = 0;
    end
    
    buffer = round(acquisition_rate * .5);
    
    % Find speed > 2
    ix_running = find(speed > 1.5);
    
    if ~isempty(ix_running)
        % Find indices that have a 1 s gap between running timepoints
        ix_start = ix_running(find(diff(ix_running) > acquisition_rate) + 1);
        if ix_running(1) > buffer
            ix_start = [ix_running(1);ix_start];
        end
        % Find onset 
        onset_indices = [];
        for k = 1:length(ix_start)
            onset_indices = [onset_indices; ix_start(k)-buffer, ix_start(k)+buffer];
        end
        if onset_indices(end,end) > length(speed)
            onset_indices(end,:) = [];
        end
    
    else
        onset_indices = [];
    end

    

    if plot_me
        figure, hold on, plot(speed,':k','LineWidth',1)
        for k = 1:size(onset_indices,1)
            ix = onset_indices(k,1):onset_indices(k,2);
            plot(ix,speed(ix),'r','LineWidth',1)
        end
    end
    
    
    