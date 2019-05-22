
% This function calculates a weighted error ratio between the errors
% during high-amplitude events to errors during baseline.
%
% Reasoning: If two ROIs have medium correlation value, how to distinguish
% between (1) two ROIs on the same axon whose correlations are corrupted by
% noise, and (2) two ROIs on different axons who are correlated? In the
% former case, would expect noise to similarly corrupt low and high
% amplitude events, but in the latter case would expect high amplitude
% events to sometimes be missed. Proposed solution is to compare accuracy
% of high amplitude events and normalize to baseline noise.
% 
% Input:
%    xx                Vector of regressor values
%    yy                Vector of data to be fit
%
% Output:
%    we                Weighted error ratio

function we = get_we(xx,yy)

% Turn into column if it isn't already a column
if isrow(xx)
    xx = xx';
end
if isrow(yy)
    yy = yy';
end

T  = size(xx,1);

%% Use linear regression on xx to predict yy
x = [xx,ones(T,1)];
y = yy;

% w is weighting factor, piecewise linear from 20th to 90th percentiles
min_ = GetSn(y)*2;%prctile(y,20);
max_ = GetSn(y)*5;%prctile(y,80);
w = y; w(w > max_) = max_; w(w<min_) = min_;
w = (w-min_)/(max_-min_);

% Normalize w
w = w/sum(w);

% Linear regression coefficients
b = (x'*x)\(x'*y);
%W = diag(w);
%b = (x'*W*x)\(x'*W*y);

% Estimate of y
yhat = x*b;

% Error
e = (y-yhat).^2 ;

% For ratio of weighted error for large events
n = sum(w .* e);

% Denominator is afiguverage error for baseline
min_ = prctile(y,10);
d = sum((y <= min_) .* e) / sum((y <= min_));

% Ratio numerator / denominator
we = (n / d); 

 