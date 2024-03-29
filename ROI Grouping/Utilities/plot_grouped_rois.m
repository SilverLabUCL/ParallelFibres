% Plots ROIs that have been grouped into putative axons
% Originally from CNMF_E code
% Updated by NACG (2019) to add variable colours

function plot_grouped_rois(Ain_axons,Cn,dFF,ix_axons_to_rois,acquisition_rate,thr,display_numbers)

    if nargin < 6 || isempty(thr)
        thr = 0.95;
    end
    
    if nargin < 7 || isempty(display_numbers)
        display_numbers = 1;
    end

    [d1,d2] = size(Cn);
    T = size(dFF,2);
    
    close all
    
    % Color map + centers
    cmap = lines(size(Ain_axons,2));
    cm = com(Ain_axons(:,1:end),d1,d2);
    
    % Spatial filters of grouped ROIs
    figure(1), hold on
    imagesc(Cn,[min(Cn(:)),max(Cn(:))]);
    colormap(gray)
    axis tight; axis equal;
    set(gca,'XTick',[],'YTick',[]);
    
    % dFF of grouped ROIs
    figure(2), hold on
    count_fig2 = 0; 
    
    % Spatial filters of loner ROIs
    figure(3), hold on
    imagesc(Cn,[min(Cn(:)),max(Cn(:))]);
    colormap(gray)
    axis tight; axis equal;
    set(gca,'XTick',[],'YTick',[]);
    
    % dFF of grouped ROIs
    figure(4), hold on
    count_fig4 = 0;
    
    for ax = 1:size(Ain_axons,2)
        % Get spatial filter
        A_temp = full(reshape(Ain_axons(:,ax),d1,d2));
        A_temp = medfilt2(A_temp,[3,3]);
        A_temp = A_temp(:);
        [temp,ind] = sort(A_temp(:).^2,'ascend');
        temp =  cumsum(temp);
        ff = find(temp > (1-thr)*temp(end),1,'first');
        
        % If spatial filter is too small
        if isempty(ff)
            error
        end
        
        % If grouped ROIs, 
        if numel(ix_axons_to_rois{ax})>1
            % plot in color in fig 1, gray in fig 3
            figure(1), contour(reshape(A_temp,d1,d2),A_temp(ind(ff)),'LineColor',cmap(ax,:), 'linewidth', 2);
            figure(3), contour(reshape(A_temp,d1,d2),A_temp(ind(ff)),':','color',[.5,.5,.5],'linewidth', 2);
            
            % plot dFF in fig 2
            figure(2), 
            for roi = ix_axons_to_rois{ax}
                plot((1:size(dFF,2))/acquisition_rate,dFF(roi,:)+count_fig2,'Color',cmap(ax,:));
                count_fig2 = count_fig2 + mean(dFF(roi,:)) * 3 + .5;
            end
            if display_numbers
                text(-10,count_fig2 -.5,num2str(ax),'color','k','fontsize',10,'fontweight','bold');
            end

            % If loners, plot in gray in fig 1, color in fig 3
        else 
            figure(1), contour(reshape(A_temp,d1,d2),A_temp(ind(ff)),':','color',[.5,.5,.5], 'linewidth', 2);
            figure(3), contour(reshape(A_temp,d1,d2),A_temp(ind(ff)),'color',cmap(ax,:),'linewidth', 2);
            
            % plot dFF in fig 2
            roi = ix_axons_to_rois{ax};
            figure(4), plot((1:size(dFF,2))/acquisition_rate,dFF(roi,:)+count_fig4,'Color',cmap(ax,:));
            count_fig4 = count_fig4 + mean(dFF(roi,:)) * 3 + .5;
            if display_numbers
                text(-10,count_fig4 -.5,num2str(ax),'color','k','fontsize',10,'fontweight','bold');
            end

        end
    end
    
    % Write numbers on top of ROIs
    if display_numbers
        for ax = 1:size(Ain_axons,2)
            if numel(ix_axons_to_rois{ax})>1
                figure(1), text((cm(ax,2)),(cm(ax,1)),num2str(ax),'color','w','fontsize',16,'fontweight','bold');
            else
                figure(3), text((cm(ax,2)),(cm(ax,1)),num2str(ax),'color','w','fontsize',16,'fontweight','bold');
            end
        end
    end
    
    figure(1), set(gca,'YDir','reverse')
    title('Map of grouped ROIs')
    
    figure(2),
    axis([-20,(T/acquisition_rate)+10,-1,count_fig2+1])
    set(gca,'YTick',0:1); set(gca,'YTickLabel',[])
    set(gca,'FontSize',14)
    title('DFF (grouped ROIs)')
    
    figure(3), set(gca,'YDir','reverse')
    title('Map of loner ROIs')
    
    figure(4), 
    axis([-20,(T/acquisition_rate)+10,-1,count_fig4+1])
    set(gca,'YTick',0:1); set(gca,'YTickLabel',[])
    set(gca,'FontSize',14)
    title('DFF (loner ROIs)')
    
end



