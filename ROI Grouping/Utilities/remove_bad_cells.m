% This function removes detected ROIs with low SNR
%
% Input:
%    Ain              *Thresholded* spatial filters matrix (num pixels x num ROIs)
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    dims             Vector of pixel dims of patch: [d1, d2]
%    min_pnr          Minimum SNR
%    acq_rate         Acquisition rate (Hz) (for plotting /manual mode)
%                     (optional - leave empty for fully automated mode)
% 
% Output:
%    A_good           New spatial filters matrix of good ROIs
%    SNR_good         Vector of SNR for good ROIs
%    A_bad            New spatial filters matrix of BAD ROIs
%    SNR_bad          Vector of SNR for BAD ROIs

function [A_good,SNR_good,A_bad,SNR_bad] = remove_bad_cells(Ain,Y,dims,acq_rate,min_pnr,manual)

if nargin < 5 || isempty(min_pnr)
    min_pnr = 6;
end
if nargin < 6 || isempty(manual)
    manual = 0;
end

d1 = dims(1);
d2 = dims(2);

%% Ensure that input has been thresholded
if sum(sum(Ain ~=1 & Ain ~= 0))
    error('Spatial matrices must be threshoded to remove bad cells');
end

%% Calculate dFF
dFF = get_dFF(Ain,Y);

%% Correlation image on raw data (unfiltered)
Cn = correlation_image(Y, [1,2], d1,d2);

%% 
SNR = zeros(size(Ain,2),1); 
goodrois = zeros(size(Ain,2),1); 

if manual == 1
    figure, hold on
end

for k = 1:size(Ain,2)
    
    dFF_ = dFF(k,:);
    
    % Estimate SNR as peak dFF over noise measured from power spectrum
    SNR(k) = max(medfilt1(dFF_,round(.2*acq_rate)))/GetSn(dFF_);
    
    % If SNR is less than min_pnr, default is to remove that ROI
    if SNR(k) < min_pnr
        remove_me = 1;
        col_plot = 'r'; % plot in red if bad ROI
    else
        remove_me = 0;
        col_plot = 'b'; % plot in blue if good ROI
    end
    
    % In case of manual mode
    if manual
        % Plot spatial filter on top of correlation image
        subplot(2,1,1), imagesc(Cn);
        caxis(median(Cn(:))+[-1,5]*std(Cn(:)))
        colormap(gray), hold on
        B = bwboundaries(reshape(Ain(:,k),d1,d2));
        for kk = 1:length(B)
           boundary = B{kk};
           plot(boundary(:,2), boundary(:,1),'Color', col_plot, 'LineWidth', 2)
        end
        hold off
        title(['ROI ',num2str(k),' of ',num2str(size(Ain,2))])
        set(gca,'FontSize',15), set(gca,'XTick',{}), set(gca,'YTick',{})
        
        % Plot dFF
        subplot(2,1,2), plot((1:length(dFF_))/acq_rate,dFF_,'Color',col_plot)
        set(gca,'FontSize',15), axis([0,length(dFF_)/acq_rate,-.5,2.5])
        title(['SNR = ',num2str(SNR(k))])
        set(gca,'FontSize',15)
        
        % Manual override
        if remove_me == 0
            u = input('This ROI will be kept. Enter R to override: ', 's');
            if strcmp(u,'R') || strcmp(u,'r')
                remove_me = 1;
                 disp(['ROI ',num2str(k),' has been removed.'])
            end
        else
            u = input('This ROI will be removed. Enter K to override: ', 's');
            if strcmp(u,'K') || strcmp(u,'k')
                remove_me = 0;
                disp(['ROI ',num2str(k),' has been kept.'])
            end
        end
        
    end
    
    % Keep list of bad ROIs
    if ~remove_me
        goodrois(k) = 1;
        %A_new(:,k) = nan(d1*d2,1);
        %SNR_good(k) = nan;
    end
end

%% Remove bad ROIs

A_good = Ain(:,goodrois==1);
SNR_good = SNR(goodrois==1);

A_bad = Ain(:,goodrois==0);
SNR_bad = SNR(goodrois==0);
 
% %% Plot results 
% 
% figure, imagesc(reshape(sum(A_good,2),d1,d2))
% set(gca, 'xtick', []); set(gca, 'ytick', []); 
% axis tight; axis equal;
% caxis([0,3])
