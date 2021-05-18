%
% This function plots an error snake (SEM) around mean values
% Input:
%    x                    vector of x values
%    y                    matrix of y values 
%    face_color           colour of error snake
% 

function plot_error_snake(x,y,face_color)

    N_trials = size(y,1);
    y_mean = nanmean(y,1);
    y_top = y_mean+nanstd(y,[],1)/sqrt(N_trials);
    y_bot = y_mean-nanstd(y,[],1)/sqrt(N_trials);
    
    ix = find(~isnan(y_top) & ~isnan(y_bot));
    y_top = y_top(ix);
    y_bot = y_bot(ix);
    x = x(ix);
    y_mean = y_mean(ix);

    fill([x, fliplr(x)],[y_top, fliplr(y_bot)],'k','FaceColor',face_color,'LineStyle','none','FaceAlpha',.3)
    hold on, plot(x,y_mean,'k')
    
    set(gca,'Box','off')
    set(gca,'FontSize',18)