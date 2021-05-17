
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
%    dFF1                ROI 1 time series
%    dFF2                ROI 2 time series
%    manual              Set to 1 to plot, 0 otherwise
%
% Output:
%    var_ratio         Ratio of variance of projected data, to variance of
%                       projected baseline data

function var_ratio = get_var_ratio(dFF1,dFF2,manual)

if nargin < 3 || isempty(manual)
    manual = 0;
end

%% Choose higher variance ROI as regressor - improves fit
if var(dFF1) > var(dFF2)
    xx = dFF1; yy = dFF2;
else
    xx = dFF2; yy = dFF1;
end

% Turn into column if it isn't already a column
if isrow(xx)
    xx = xx';
end
if isrow(yy)
    yy = yy';
end

% Data matrix
X = [xx,yy]; % Tx2 matrix

% Fit 2 Gaussians to capture the baseline data + events
options = statset('MaxIter',1000); % Increase number of EM iterations
gmfit = fitgmdist(X,2,'CovarianceType','diagonal','SharedCovariance',false,'Options',options);

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

S_bl = diag(S_bl);

% Find scaling vector u via regression coefficients
x = [xx,ones(size(xx))];
b = (x'*x)\(x'*yy);
u = [1;b(1)]; u = u/norm(u);
u0 = u;
v0 = [-u0(2);u0(1)];

% Find minimum projection onto a vector
% This is the vector we will project onto
f = @(vv)var(X*[vv(1);vv(2)]/norm(vv));
v = fminsearch(f,v0);
v = v/norm(v);

% Ratio of variances
% All data
var_proj_data = var(X*v); % variance of projection of X onto v
% Baseline distribution
mean_proj_baseline = (m_bl*v)'; 
var_proj_baseline = v'*S_bl*v; % variance of baseline gaussian projected onto v
var_ratio = var_proj_data / var_proj_baseline;

% Plot everything 
if manual
    figure(2), hold off, plot(xx,yy,'-','Color',[.3,.3,.3]), hold on
    
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
    plot(m_bl(1)+[0,u(1)]/3,m_bl(2)+[0,u(2)]/3,'k','LineWidth',3)
    plot(m_bl(1)+[0,v(1)]/3,m_bl(2)+[0,v(2)]/3,'r','LineWidth',3)
    
    title(var_ratio), axis equal
    set(gca,'FontSize',15)
    xlabel('\Delta F/F ROI 1')
    ylabel('\Delta F/F ROI 2')
    
    figure(3), hold off, histogram(X*v,'Normalization','pdf','FaceColor',[.3,.3,.3]), hold on
    plot(-.5:.01:.5,normpdf(-.5:.01:.5,mean_proj_baseline,sqrt(var_proj_baseline)),'k','LineWidth',2)
    set(gca,'FontSize',15)
    xlim([-.5,.5])
    ylabel('PDF')
    xlabel('Projected value')
end

