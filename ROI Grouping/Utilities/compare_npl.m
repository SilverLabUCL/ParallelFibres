% This function compares to neuropil
%
% Input:
%    Ain              *Thresholded* spatial filters matrix (num pixels x num ROIs)
%    Ain_npl          Thresholded spatial filters for neuropil subtraction
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    acq_rate         Acquisition rate (Hz) (optional)
%    smooth_win_s     Smoothing window (s) (optional)
%
% Output:
%    F_raw            Matrix of raw fluorescence (optional)
%    F_npl            Matrix of neuropil fluorescence (optional)


function [F_raw,F_npl] = compare_npl(Ain,Ain_npl,Y,acq_rate,smooth_win_s)

% Default no smoothing 
if nargin < 5 || isempty(smooth_win_s) || smooth_win_s == 0
    smooth_win_s = [];
end

% If acq_rate is not given, do not smooth
if nargin < 4 || isempty(acq_rate)
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

F_raw = zeros(N_ROIs,T);
F_npl = zeros(N_ROIs,T);

for k = 1:N_ROIs
    
        % Average F over all pixels within masks 
        F_raw(k,:) = mean(Y(Ain(:,k)==1,:),1);
        F_npl(k,:) = mean(Y(Ain_npl(:,k)==1,:),1);
    
    % Smooth dFF 
    if ~isempty(smooth_win_s)
        % Asymmetric Gaussian filter
        % Have to double smoothing window because asymmetric Gaussian only
        % uses one side
        F_raw(k,:) = smoothdata(F_raw(k,:),'gaussian', [2*smooth_bins 0]);
        F_npl(k,:) = smoothdata(F_npl(k,:),'gaussian', [2*smooth_bins 0]);
    end
    
end