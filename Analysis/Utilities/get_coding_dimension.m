
% This function calculates the vector separating QW and A periods/manifolds
% aka the coding dimension (following Li et al 2016)
% Used for coloring plots since A/QW segmentation doesn't label all timepts
%
% Input:
%    dFF      Population activity
%    A        Start & end times of active periods
%    QW       Start & end times of QW periods
%
% Output:
%    w        Normalized vector corresponding tp the separation between A 
%             and QW manifolds (N_ROIs x 1)

function w = get_coding_dimension(dFF,A,QW)

    % All indices for A 
    A_ix = [];
    for k = 1:length(A)
        A_ix = [A_ix, A(k,1):A(k,2)];
    end

    % All indices for QW
    QW_ix = [];
    for k = 1:length(QW)
        QW_ix = [QW_ix, QW(k,1):QW(k,2)];
    end

    % Mean of each manifold in activity space
    mu_A = mean(dFF(:,A_ix),2);
    mu_QW = mean(dFF(:,QW_ix),2);

    % Vector separating means of QW / A
    w = mu_A - mu_QW;
    w = w / norm(w); % normalize
