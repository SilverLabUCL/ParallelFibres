% This function calculates distribution of correlations by bootstrapping
% random contiguous intervals throughout - used for standard deviation of
% correlations between dFF of two ROIs to determine whether to group them
%
% Input:
%    xx             dFF of ROI 1
%    yy             dFF of ROI 2
%    acq_rate       Acquisition rate (Hz)
%    time_interval  Interval length (optional)
% 
% Output:
%    C              Vector of bootstrapped sample correlations

function C = get_corr_dist(xx,yy,acq_rate,time_interval)
    
    if nargin < 4 || isempty(time_interval)
        time_interval = 20;
    end
    
    T = length(xx);
    K = round(time_interval * acq_rate);
    
    %% Get noise & smoothed traces to later measure SNR
    % Smoothed versions
    xx_smooth = medfilt1(xx,round(.2*acq_rate));
    yy_smooth = medfilt1(yy,round(.2*acq_rate));
        
    % Noise levels
    std_xx = GetSn(xx);
    std_yy = GetSn(yy);
    
    % Global SNR
    xx_snr = max(xx_smooth)/std_xx;
    yy_snr = max(yy_smooth)/std_yy;
    
    %% Bootstrap samples and calculate correlations
    num_samples = 1000;
    C = nan(1,num_samples);
    %figure
    for t = 1:num_samples
        ix_start = randsample(1:(T-K-1),1);
        ix = ix_start : (ix_start+K-1);
        
        % Baseline for this snippet
        xx_baseline = prctile(xx_smooth(ix),10);
        yy_baseline = prctile(yy_smooth(ix),10);
        
        % SNR for this snippet
        xx_snr_ix = max(xx_smooth(ix)-xx_baseline)/std_xx;
        yy_snr_ix = max(yy_smooth(ix)-yy_baseline)/std_yy;

        % If SNR for this snippet is 25% of total SNR, calculate
        % correlations - otherwise, there may be no event. This improves
        % calculation of correlations for ROIs that may be sparsely active
        if xx_snr_ix > xx_snr/3 || yy_snr_ix > yy_snr/3
            temp = [xx(ix);  yy(ix)];
            C_temp = corrcoef(temp');
            C_temp = C_temp(1,2);
            C(t) = C_temp;
            % Following code is for plotting intervals 
%             hold off, plot(xx), hold on, plot(yy+1)
%             plot([ix(1),ix(1)],[-.5,3],'k','LineWidth',3)
%             plot([ix(end),ix(end)],[-.5,3],'k','LineWidth',3),
%             C_temp
%             pause
%         else
%             hold off, plot(xx), hold on, plot(yy+1)
%             plot([ix(1),ix(1)],[-.5,3],'r','LineWidth',3)
%             plot([ix(end),ix(end)],[-.5,3],'r','LineWidth',3),
%             C_temp
%             pause
        end
    end
    
    %%
    


    
    
    
    