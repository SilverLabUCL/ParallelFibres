%
% This function allows the user to delete ROIs by eye, specifically useful
% for bright debris that is mistakenly classified to an ROI because of
% bleaching
%
% Input:
%    Ain              Spatial filters matrix (num pixels x num ROIs)
%    dFF              Matrix of dFF (num ROIs x num timepoints)
% 
% Output:
%    Ain_new          Ain after removing offending ROIs
%    dFF_new          dFF after removing offending ROIs

function [Ain_new,dFF_new] = delete_ROIs_manually(Ain,dFF,acquisition_rate)
    
close all
    [N,T] = size(dFF);
    
    ROIs_to_delete = [];
    
    figure
    y_pos = 0;
    for it = 0:floor(N/10)-1
        ixs = 10*it+1 : 10*(it+1);
        hold off
        for roi = ixs
            plot((1:size(dFF,2))/acquisition_rate,dFF(roi,:)+y_pos); hold on
            text(-15,y_pos +.5,num2str(roi),'color','k','fontsize',10,'fontweight','bold');
            y_pos = y_pos + mean(dFF(roi,:)) * 3 + .5;
        end
        xlim([-20,T/acquisition_rate+20])
        
        stop = 0;
        while stop == 0
            in = input('Which ROIs should be deleted? Enter A to advance: ','s');
            if strcmp(in,'A') || strcmp(in,'a')
                    stop = 1;
            elseif all(ismember(in, '0123456789')) 
                if strcmp(in,'')
                    disp('Enter something!')
                elseif str2num(in) > ixs(end) || str2num(in) < ixs(1)
                    disp('Out of bounds');
                else
                    ROIs_to_delete = [ROIs_to_delete; str2num(in)];
                end
            else
                disp('That is not an ROI.')
            end
        end
    end
    
    hold off
    ixs = floor(N/10)*10+1:N;
    for roi = ixs
        plot((1:size(dFF,2))/acquisition_rate,dFF(roi,:)+y_pos); hold on
        text(-15,y_pos +.5,num2str(roi),'color','k','fontsize',10,'fontweight','bold');
        y_pos = y_pos + mean(dFF(roi,:)) * 2 + .5;
    end
    xlim([-20,T/acquisition_rate+20])
        
    stop = 0;
    while stop == 0
        in = input('Which ROIs should be deleted? Enter A to advance: ','s');
        if  strcmp(in,'A') || strcmp(in,'a')
                stop = 1;
        elseif all(ismember(in, '0123456789')) 
            if strcmp(in,'')
                    disp('Enter something!')
            elseif str2num(in) > ixs(end) || str2num(in) < ixs(1)
                disp('Out of bounds');
            else
                ROIs_to_delete = [ROIs_to_delete; str2num(in)];
            end
        else
            disp('That is not an ROI.')
        end
    end
    
    % Delete bad ROIs
    
    Ain_new = Ain;
    Ain_new(:,ROIs_to_delete) = [];
    
    dFF_new = dFF;
    dFF_new(ROIs_to_delete,:) = [];

    
    
    