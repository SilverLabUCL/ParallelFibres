%
% This function calculates correlations for axon grouping
% Any correlation less than a threshold (rho_min) is set to 0
%
% Input:
%    dFF                Raw fluorescence matrix (num pixels x num timepoints)
%    rho_min            Minimum correlation coefficient to consider
% 
% Output:
%    C            Correlation

function C = get_corrs(dFF,rho_min)
    
    C = corrcoef(dFF');
    C = C-diag(diag(C)); % remove ones on the diagonal
    C = max(C,rho_min);  % remove any corerlations less than rho_thresh
    C(C==rho_min) = 0; % set min corrs to 0 
    C = triu(C,1); %  Keep only upper triangular part
    