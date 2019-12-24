
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

function notR = define_not_running(speed,acquisition_rate,plot_me)
    
    if nargin < 3 || isempty(plot_me)
        plot_me = 0;
    end

    % Long smoothing window (500 ms)
    smooth_win = round(acquisition_rate * .5);

    % Running information based on wheel MI
    y = abs(diff(speed));
    y = smoothdata(y,'movmedian',smooth_win);
    y = y - mode(round(y,4));
    y = y / nanstd(y);

    % Find periods when not R below .1 STD
    notR = (y < .1);

    % Plot to manually check behavioural labels
    if plot_me
        t = (1:length(y))/acquisition_rate;
        figure,
        plot(t,y,'k','LineWidth',1.5), hold on
        hold on, plot(t(notR),y(notR),'ob','MarkerFaceColor','b')
        ylabel('Loc'), set(gca,'FontSize',15)
        xlim([0,length(speed)/acquisition_rate])
        xlabel('Time (s)')

    end
