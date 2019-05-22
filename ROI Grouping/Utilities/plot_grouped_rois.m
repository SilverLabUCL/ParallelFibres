% Originally from CNMF_E code
% Updated by NACG (2019) to add variable colours

function plot_grouped_rois(Aor,Cn,ix_axons_to_rois,thr,display_numbers)

    [d1,d2] = size(Cn);
    imagesc(Cn,[min(Cn(:)),max(Cn(:))]); hold on
    colormap(gray)
    axis tight; axis equal;
    set(gca,'XTick',[],'YTick',[]);
    posA = get(gca,'position');
    set(gca,'position',posA);

    cmap = lines(size(Aor,2));
    cm = com(Aor(:,1:end),d1,d2);
    i = 1;
    for ax = 1:size(Aor,2)
        A_temp = full(reshape(Aor(:,ax),d1,d2));
        A_temp = medfilt2(A_temp,[3,3]);
        A_temp = A_temp(:);
        [temp,ind] = sort(A_temp(:).^2,'ascend');
        temp =  cumsum(temp);
        ff = find(temp > (1-thr)*temp(end),1,'first');
        if ~isempty(ff)
            if numel(ix_axons_to_rois{ax})>1
                contour(reshape(A_temp,d1,d2),A_temp(ind(ff)),'LineColor',cmap(i,:), 'linewidth', 2);
                i = i+1;
            else 
                contour(reshape(A_temp,d1,d2),A_temp(ind(ff)),':','color',[.5,.5,.5], 'linewidth', 2);
            end
        end
        hold on;
    end
    
    if display_numbers
        for ax = 1:size(Aor,2)
            if numel(ix_axons_to_rois{ax})>1
            	text((cm(ax,2)),(cm(ax,1)),num2str(ax),'color','w','fontsize',16,'fontweight','bold');
            end
        end
    end
    
    
end



