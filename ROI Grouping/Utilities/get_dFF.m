% This function removes detected ROIs with low SNR
%
% Input:
%    Ain              *Thresholded* spatial filters matrix (num pixels x num ROIs)
%    Ain_npl          Thresholded spatial filters for neuropil subtraction
%    Y                Raw fluorescence matrix (num pixels x num timepoints)
%    acq_rate         Acquisition rate (Hz) (optional)
%    smooth_win_s     Smoothing window (s) (optional)
%    npl_frac         Fraction of neuropil to remove (optional)
%
% Output:
%    dFF              Matrix of dFF (num ROIs x num timepoints)
%    F_raw            Matrix of raw fluorescence (optional)
%    F_npl            Matrix of neuropil fluorescence (optional)


function [dFF,F_raw,F_npl] = get_dFF(Ain,Ain_npl,Y,acq_rate,smooth_win_s,npl_frac)

% Default smoothing for gcamp6f timescale
if nargin < 5 || isempty(smooth_win_s)
    smooth_win_s = .2; % 200 ms
end

% Default - subtract 70% of neuropil
if nargin < 6 || isempty(npl_frac)
    npl_frac = 0.7; 
end

if npl_frac <0 || npl_frac > 1
    error('Fraction of neuropil removed must be between 0 and 1');
end

% If smooth_win_s is set to 0 (i.e., no smoothing)
% Then replace with [] (otherwise smooth function will throw an error)
if smooth_win_s == 0
    smooth_win_s = [];
end

% If acq_rate is not given, do not smooth
if nargin < 4 || acq_rate == []
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
dFF = zeros(N_ROIs,T);

for k = 1:N_ROIs
    
        % Average F over all pixels within masks 
        F_raw(k,:) = mean(Y(Ain(:,k)==1,:),1);
        F_npl(k,:) = mean(Y(Ain_npl(:,k)==1,:),1);
        
        % Subtract neuropil
        F = F_raw(k,:) - npl_frac * F_npl(k,:);

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