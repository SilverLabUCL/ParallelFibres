% This function calculates the dFF for a hand-drawn region of the patch
%
% Input:
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    dims             Vector of pixel dims of patch: [d1, d2]
%    acq_rate         Acquisition rate (Hz)
% 
% Output:
%    A_neuropil       Spatial mask for hand drawn neuropil region
%    dFF_neuropil     Neuropil dFF

function [A_neuropil,dFF] = get_neuropil_byhand(Ain,Y,dims,acq_rate,smooth_win_s)


if nargin < 5 || isempty(smooth_win_s)
    smooth_win_s = 0.2;
end

d1 = dims(1);
d2 = dims(2);

% Correlation image
Cn = correlation_image(Y, [1,2], d1,d2);

% Plot correlation image of raw data
figure, plot_contours(Ain,Cn,.9);
colormap(gray); caxis(median(Cn(:))+[-1,1]*std(Cn(:)))
set(gca, 'xtick', []); set(gca, 'ytick', []);
axis tight; axis equal;

% Query user to hand draw ROI
fh = imfreehand();
drawnow; pause(0.1);
x = fh.getPosition();
A_neuropil = poly2mask(x(:,1),x(:,2),d1,d2);
A_neuropil = A_neuropil(:);

% Get dFF of neuropil (unsmoothed)
dFF = get_dFF(A_neuropil,Y,smooth_win_s);

% Estimate SNR
SNR =  max(medfilt1(dFF,round(.2*acq_rate)))/GetSn(dFF);
    

% Plot neuropil region
figure, subplot(2,1,1), imagesc(Cn);
caxis(median(Cn(:))+[-1,5]*std(Cn(:)))
colormap(gray), hold on
B = bwboundaries(reshape(A_neuropil,d1,d2));
for kk = 1:length(B)
   boundary = B{kk};
   plot(boundary(:,2), boundary(:,1),'r', 'LineWidth', 2)
end
set(gca,'FontSize',15), set(gca,'XTick',{}), set(gca,'YTick',{})

% Plot dFF (blue) and smoothed dFF (red) 
subplot(2,1,2), plot((1:length(dFF))/acq_rate,dFF,'r')
set(gca,'FontSize',15)
title(['SNR = ',num2str(SNR)])
set(gca,'FontSize',15)