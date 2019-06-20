
% This function calculates the ratio between the variance of the deviations
% from the linear fit between to ROIs, to the same variance only for the
% baseline activity.
%
% First, 'baseline' activity for the two ROIs is determined by fitting a
% mixture of two gaussians to the time series of both ROIs (k=2). The
% gaussian with lower mean is considered the baseline. The variance ratio
% is the variance of all data points projected onto the vector that is the
% best linear fit between the ROIs, divided by the variance from just the
% baseline activity.
%
% Input:
%    xx                ROI 1 time series
%    yy                ROI 2 time series
%    manual            Set to 1 to plot, 0 otherwise
%
% Output:
%    var_ratio         Ratio of variance of projected data, to variance of
%                       projected baseline data

function var_ratio = get_var_ratio(xx,yy,manual)

if nargin < 3 || isempty(manual)
    manual = 0;
end

% Turn into column if it isn't already a column
if isrow(xx)
    xx = xx';
end
if isrow(yy)
    yy = yy';
end

% Data matrix
X = [xx,yy];

% Fit 2 Gaussians to capture the baseline data + events
options = statset('MaxIter',1000); % Increase number of EM iterations
gmfit = fitgmdist(X,2,'CovarianceType','full','SharedCovariance',false,'Options',options);

% Identify baseline distribution
test = (gmfit.mu(2,:) > gmfit.mu(1,:)) ;
if test(1) == 1 && test(2) == 1
    S_bl = gmfit.Sigma(:,:,1); % baseline cov matrix
    m_bl = gmfit.mu(1,:); % baseline mean
elseif test(1) == 0 && test(2) == 0
    S_bl = gmfit.Sigma(:,:,2); % baseline cov matrix
    m_bl = gmfit.mu(2,:); % baseline mean
else % If no real baseline for both ROIs
    error
end

% Find scaling vector u via regression coefficients
x = [xx,ones(size(xx))];
b = (x'*x)\(x'*yy);
u = [1;b(1)]; u = u/norm(u);

% Orthogonal vector to u
% This is the vector we will project onto
v = [-u(2);u(1)]; v = v/norm(v);

% Ratio of variances
num = var(X*v); % variance of 
den = v'*S_bl*v; % variance of baseline gaussian projected onto v
var_ratio = num / den;

% Plot everything 
if manual
    figure(2), hold off, plot(xx,yy,'.'), hold on
    
    % To plot error ellipse
    % Get eigenvectors & eigenvalues of the fitted baseline distribution
    [V_bl, D_bl] = eig(S_bl);

    % Sort eigenvalues in descending order
    [D_bl,ind] = sort(diag(D_bl),'descend');
    V_bl = V_bl(:,ind);

    % Calculate the angle between the x-axis and the principal eigenvector
    angle_bl = atan2(V_bl(2,1), V_bl(2,2));

    % Get the 95% confidence interval error ellipse
    chisquare_val = sqrt(chi2inv(.95, 2)) ;
    theta_grid = linspace(0,2*pi);
    a_bl = chisquare_val*sqrt(D_bl(1));
    b_bl = chisquare_val*sqrt(D_bl(2));

    r_ellipse_bl = [a_bl*cos( theta_grid );b_bl*sin( theta_grid )]' *  [ cos(angle_bl) sin(angle_bl); -sin(angle_bl) cos(angle_bl) ];
    plot(r_ellipse_bl(:,1) + m_bl(1),r_ellipse_bl(:,2) + m_bl(2),'k','LineWidth',2)
    
    % Plot u and v vetors
    plot(m_bl(1)+[0,u(1)],m_bl(2)+[0,u(2)],'k','LineWidth',3)
    plot(m_bl(1)+[0,v(1)],m_bl(2)+[0,v(2)],'r','LineWidth',3)
    
    title(var_ratio), axis tight equal
    set(gca,'FontSize',18)
end


%%
% % Centre and rotate data 
% Y = (X-m_bl) * R';
% 
% % Normalize lengths
% Y(:,1) = Y(:,1) / a_bl;
% Y(:,2) = Y(:,2) / b_bl;
% 
% % Points that are outside the confidence interval
% ix = find( Y(:,1).^2 + Y(:,2).^2 > 1);
% 
% % Find scaling vector via regression coefficients
% x = [xx(ix),ones(numel(ix),1)];
% b = (x'*x)\(x'*yy(ix));

%%
% T  = size(xx,1);
% X = [xx,yy];
% 
% % Fit 2 Gaussians to capture the baseline data + 
% options = statset('MaxIter',1000); % Increase number of EM iterations
% gmfit = fitgmdist(X,2,'CovarianceType','full','SharedCovariance',false,'Options',options);
% 
% % Get eigenvalues
% test = (gmfit.mu(2,:) > gmfit.mu(1,:)) ;
% if test(1) == 1 && test(2) == 1
%     S_bl = gmfit.Sigma(:,:,1); % baseline Gaussian
%     S_event = gmfit.Sigma(:,:,2); % large event Gaussian
%     m_bl = gmfit.mu(1,:);
%     m_event = gmfit.mu(2,:);
% elseif test(1) == 0 & test(2) == 0
%     S_event = gmfit.Sigma(:,:,1); % baseline Gaussian
%     S_bl = gmfit.Sigma(:,:,2); % large event Gaussian
%     m_event = gmfit.mu(1,:);
%     m_bl = gmfit.mu(2,:);
% else
%     error
% end
% 
% % Get eigenvectors & eigenvalues of the fitted distriutions
% [V_bl, D_bl] = eig(S_bl);
% [V_event, D_event] = eig(S_event);
% 
% % Sort eigenvalues in descending order
% [D_bl,ind] = sort(diag(D_bl),'descend');
% V_bl = V_bl(:,ind);
% 
% [D_event,ind] = sort(diag(D_event),'descend');
% V_event = V_event(:,ind);
% 
% % Calculate the angle between the x-axis and the largest eigenvector
% angle_bl = atan2(V_bl(2,1), V_bl(2,2));
% angle_event = atan2(V_bl(2,1), V_bl(2,2));
% 
% % This angle is between -pi and pi.
% % Let's shift it such that the angle is between 0 and 2pi
% if(angle_bl < 0)
%     angle_bl = angle_bl + 2*pi;
% end
% if(angle_event < 0)
%     angle_event = angle_event + 2*pi;
% end
% 
% % Get the 95% confidence interval error ellipse
% chisquare_val = 2;
% theta_grid = linspace(0,2*pi);
% a_bl = chisquare_val*sqrt(D_bl(1));
% b_bl = chisquare_val*sqrt(D_bl(2));
% a_event = chisquare_val*sqrt(D_event(1));
% b_event = chisquare_val*sqrt(D_event(2));
% 
% %let's rotate the ellipse to some angle phi
% r_ellipse_bl = [a_bl*cos( theta_grid );b_bl*sin( theta_grid )]' *  [ cos(angle_bl) sin(angle_bl); -sin(angle_bl) cos(angle_bl) ];
% r_ellipse_event = [a_event*cos( theta_grid );b_event*sin( theta_grid )]' *  [ cos(angle_event) sin(angle_event); -sin(angle_event) cos(angle_event) ];
% 
% % Draw the error ellipse
% figure, hold on, plot(xx,yy,'.')
% plot(r_ellipse_bl(:,1) + m_bl(1),r_ellipse_bl(:,2) + m_bl(2),'k','LineWidth',2)
% plot(r_ellipse_event(:,1) + m_event(1),r_ellipse_event(:,2) + m_event(2),'k','LineWidth',2)
% title(D_event(2)/mean(D_bl))
% axis tight equal