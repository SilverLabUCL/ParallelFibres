



function plot_error_snake(x,y,face_color)

    N_trials = size(y,1);
    y_top = nanmean(y,1)+std(y,[],1)/sqrt(N_trials);
    y_bot = nanmean(y,1)-std(y,[],1)/sqrt(N_trials);

    fill([x, fliplr(x)],[y_top, fliplr(y_bot)],'k','FaceColor',face_color,'LineStyle','none','FaceAlpha',.3)
    hold on, plot(x,mean(y,1),'k')
    
    set(gca,'Box','off')
    set(gca,'FontSize',18)