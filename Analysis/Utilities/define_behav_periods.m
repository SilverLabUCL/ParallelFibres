
% This function defines QW and A periods based on behavioural data
%
% Input:
%    whisk_amp            Whisking amplitude
%    speed                Wheel MI (locomotion)
%    acquisition_rate     Acquisition rate (Hz)
%    plot_me              Optional flag to plot results from definition of A vs QW
%
% Output:
%    A        Matrix defining indices of active running + whisking
%             1st column start times, 2nd column end times (indices, not
%             real time)
%    QW       Corresponding indices for QW periods

function [A,QW] = define_behav_periods(whisk_amp,speed,acquisition_rate,plot_me)
    
    if nargin < 4 || isempty(plot_me)
        plot_me = 0;
    end

    % Long smoothing window (500 ms)
    smooth_win = round(acquisition_rate * .5);

    % First normalize behavioural data (sortof zscore w/mode instead of
    % median bc mode typically corresponds to QW periods)
    % Whisking information based on whisker amplitude
    x = whisk_amp;
    x = smoothdata(x,'movmedian',smooth_win);
    x = x - mode(round(x,4));
    x = x / nanstd(x);

    % Running information based on wheel MI
    y = abs(diff(speed));
    y = smoothdata(y,'movmedian',smooth_win);
    y = y - mode(round(y,4));
    y = y / nanstd(y);

    % Make vectors equal length since speed is based on diff
    x = x(1:end-1);

    % Find Active periods (R+W)
    % Criteria 1: both R and W are above .1 STD
    RandW = (x >= .1 & y >= .1);
    ix_start = find(diff(RandW)>0);
    ix_end = find(diff(RandW)<0);

    % Correct if starts with running period
    if ix_start(1) > ix_end(1)
        ix_start = [1; ix_start];
    end

    % Correct if ends with running period
    if ix_start(end) > ix_end(end)
        ix_end = [ix_end;length(x)];
    end

    % Criteria 2: Delete any that are < 3 s long
    A = [ix_start, ix_end];
    for k = 1:size(A,1)
        if A(k,2) - A(k,1) < 3* acquisition_rate
            A(k,:) = [nan,nan];
        end
    end
    A = A(~isnan(sum(A,2)),:);

    % Find Quiet Wakefulness periods (noW,noR)
    % Criteria 1: both R and W are below .1 STD
    notRW = (x < .1 & y < .1);
    ix_start = find(diff(notRW)>0);
    ix_end = find(diff(notRW)<0);

    % Correct if starts with SW period
    if ix_start(1) > ix_end(1)
        ix_start = [1; ix_start];
    end

    % Correct if ends with SW period
    if ix_start(end) > ix_end(end)
        ix_end = [ix_end;length(x)];
    end

    QW = [ix_start, ix_end];
    QW = QW(~isnan(sum(QW,2)),:);

    % Plot to manually check behavioural labels
    if plot_me
        figure,
        subplot(2,1,1), plot((1:length(x))/acquisition_rate,x,'k','LineWidth',1.5), 
        hold on, ylabel('Whisk'), set(gca,'FontSize',15)
        xlim([0,length(speed)/acquisition_rate])
        subplot(2,1,2), plot((1:length(y))/acquisition_rate,y,'k','LineWidth',1.5), hold on
        hold on, ylabel('Loc'), set(gca,'FontSize',15)
        xlim([0,length(speed)/acquisition_rate])
        xlabel('Time (s)')

        for k = 1:size(A,1)
            ix = A(k,1):A(k,2);
            subplot(2,1,1), plot(ix/acquisition_rate,x(ix),'Color',[1,0,1],'LineWidth',1.5)
            subplot(2,1,2), plot(ix/acquisition_rate,y(ix),'Color',[1,0,1],'LineWidth',1.5)
        end
        for k = 1:size(QW,1)
            ix = QW(k,1):QW(k,2);
            subplot(2,1,1), plot(ix/acquisition_rate,x(ix),'Color',[0,1,1],'LineWidth',1.5)
            subplot(2,1,2), plot(ix/acquisition_rate,y(ix),'Color',[0,1,1],'LineWidth',1.5)
        end
    end
