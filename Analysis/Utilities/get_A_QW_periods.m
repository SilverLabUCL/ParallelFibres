% This function separates into periods in A or QW state based on the first
% principle component being above / below its median value, plus a buffer 

function [A,QW] = get_A_QW_periods(pc1,buffer,manual_check)

    if nargin < 3 || isempty(manual_check)
        manual_check = 0;
    end
    
    signal = pc1;
    signal(signal<=0) = 0;

    A = [];
    time = (1:length(pc1))';
    
    isDone = 0;
    A_end = 0;
    while isDone == 0

        % Beginning of whisking period defined by pc1
        A_start = find(time>=A_end & signal>0);

        if numel(A_start) == 0
            isDone = 1; % no whisking periods found      
        else
            A_start = time(A_start(1));

            A_end = find(time>A_start & signal==0);

            if numel(A_end) == 0
                A = [A;A_start,time(end)];
                isDone = 1;
            else
                A_end = time(A_end(1));
                A = [A;A_start,A_end];
            end
        end
    end

    QW = [];
    if numel(A)==0
        QW = [time(1),time(end)];
    else            
        if A(1,1) > time(1)
            QW = [QW; time(1), A(1,1)];
        end
        for k = 1:size(A,1)-1
            QW = [QW; A(k,2), A(k+1,1)];
        end 
        if A(end,2) < time(end)
            QW = [QW; A(end,2), time(end)];
        end
    end


    for k = 1:length(A)
        if A(k,1) <= time(1)
        else
            A(k,1) = A(k,1) + buffer;
            if A(k,1) > A(k,2)
                A(k,:) = NaN;
            end
        end
    end
    A = A(~isnan(A(:,1)),:);
    for k = 1:length(A)
        if A(k,2) >= time(end)
        else
            A(k,2) = A(k,2) - buffer;%1;
            if A(k,1) > A(k,2)
                A(k,:) = NaN;
            end
        end
    end
    A = A(~isnan(A(:,1)),:);

    for k = 1:length(QW)
        if QW(k,1) <= QW(1)
        else
            QW(k,1) = QW(k,1) + buffer;%1;
            if QW(k,1) > QW(k,2)
                QW(k,:) = NaN;
            end
        end
    end
    QW = QW(~isnan(QW(:,1)),:);
    for k = 1:length(QW)
        if QW(k,2) >= time(end)
        else
            QW(k,2) = QW(k,2) - buffer;
            if QW(k,1) > QW(k,2)
                QW(k,:) = NaN;
            end
        end
    end
    QW = QW(~isnan(QW(:,1)),:);

    if manual_check
        figure, subplot(2,1,1)
        for k = 1:length(QW)
            x = QW(k,1);
            y = -3;
            width = QW(k,2) - QW(k,1);
            height = 6;
            hold on, rectangle('Position',[x y width height],'FaceColor',[0 1 1 ],'EdgeColor','w')
        end
        plot(time,pc1,'k')  
        set(gca,'FontSize',15)
        subplot(2,1,2)
        for k = 1:length(A)
            x = A(k,1);
            y = -3;
            width = A(k,2) - A(k,1);
            height = 6;
            hold on, rectangle('Position',[x y width height],'FaceColor',[1 .5 1],'EdgeColor','w')
        end
        plot(time,pc1,'k')
        set(gca,'FontSize',15)
    end