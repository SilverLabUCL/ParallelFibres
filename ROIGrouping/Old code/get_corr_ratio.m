%
% This function calculates C_ratio matrix to be used for ROI grouping
%
% Input:
%    dFF                Raw fluorescence matrix (num pixels x num timepoints)
%    rho_min            Minimum correlation coefficient to consider
%    K                  Number of segments to split data for calculation of
%                       the correlation coefficient std (or, if K = 0 or 1,
%                       just returns the raw correlation coefficient)
% 
% Output:
%    C_ratio            Matrix of the following ratio:
%                       (avg CC) / (std of CC) for each pair of ROIs
%                       (where CC = correlation coefficient)

function C_ratio = get_corr_ratio(dFF,rho_min,K)
    
    if nargin < 2 || isempty(K)
        K = 4;
    end
    
    [N_ROIs,T] = size(dFF);
    
    C = corrcoef(dFF');
    C = C-diag(diag(C)); % remove ones on the diagonal
    C = max(C,rho_min);  % remove any corerlations less than rho_thresh
    %C = C-tril(C); % Remove redundancy
    C(C==rho_min) = 0; % set min corrs to 0 
    
    if K > 1

        C_std = zeros(N_ROIs^2,K);

        for t = 0:K-1
            temp = dFF(:, round(t/K *T)+1 : round((t+1)/K *T)  );
            C_std_temp = corrcoef(temp');
            C_std(:,t+1) = C_std_temp(:);
        end

        C_std = std(C_std,[],2); C_std = reshape(C_std,N_ROIs,N_ROIs);

        C_ratio = C./C_std;
        
    else
        
        C_ratio = C;
        
    end
    
    %  Keep only upper triangular part
    C_ratio = triu(C_ratio,1); 
    