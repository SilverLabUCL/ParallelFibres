% This function removes detected ROIs with low SNR
%
% Input:
%    Ain              *Thresholded* spatial filters matrix (num pixels x num ROIs)
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    acq_rate         Acquisition rate (Hz) (optional)
%    smooth_win_s     Smoothing window (s) (optional)
%
% Output:
%    dFF              Matrix of dFF (num ROIs x num timepoints)


function dFF = get_dFF(Ain,Y,acq_rate,smooth_win_s)

% Default smoothing for gcamp6f timescale
if nargin < 4 || isempty(smooth_win_s)
    smooth_win_s = .2; % 200 ms
end

% If smooth_win_s is set to 0 (i.e., no smoothing)
% Then replace with [] (otherwise smooth function will throw an error)
if smooth_win_s == 0
    smooth_win_s = [];
end

% If acq_rate is not given, do not smooth
if nargin < 3
    smooth_win_s = []; 
    acq_rate = [];
end

%% Parameters

N_ROIs = size(Ain,2);
T = size(Y,2);
smooth_bins = round(acq_rate * smooth_win_s);

%% Ensure that input has been thresholded
if sum(sum(Ain ~=1 & Ain ~= 0))
    error('Spatial matrices must be threshoded to calculate dFF');
end

%% Calculate dFF

dFF = zeros(N_ROIs,T);

for k = 1:N_ROIs
    
        % Average F over all pixels within that ROI
        F = mean(Y(Ain(:,k)==1,:),1);

        % Then calculate dFF
        F0 = prctile(F,10);
        dFF(k,:) = (F - F0)/F0;
        
    
    % Smooth dFF 
    if ~isempty(smooth_win_s)
        % Asymmetric Gaussian filter
        % Have to double smoothing window because asymmetric Gaussian only
        % uses one side
        dFF(k,:) = smoothdata(dFF(k,:),'gaussian', [2*smooth_bins 0]);
    end
    
end