% This function uses linear regression on the (thresholded) 
% spatial filters to infer the fibre direction
%
% Input:
%    A   (d1*d2) x 1 spatial vector of ROIs
%    d   vector of patch size (d=[d1,d2])
%
% Output:
%    vector_axon    vector giving inferred direction of the fibre
%    x_intercept    x-intercept, so [x_intercept,0] is a point on the fibre
%    endpoints      2x2 matrix, of start (row 1) and end (row 2) points for
%                    plotting during regrouping

function [vector_axon,x_intercept,endpoints] = get_axon_vector(A,d)

d1 = d(1);
d2 = d(2);

A = reshape(A,d1,d2);
[y,x] = find(A>0);
p = polyfit(y,x,1);
rmse = sqrt(sum((p(1)*y+p(2) - x).^2));

% also try transposed version - this wors better when rois are horizontal
p2 = polyfit(x,y,1);
rmse2 = sqrt(sum((p2(1)*x+p2(2) - y).^2));

if rmse <= rmse2

    y = [0,d1];
    x = p(1)*y+p(2);
    vector_axon = [diff(x),diff(y)];
    vector_axon = vector_axon/norm(vector_axon);
    
else
    % if rmse is lower for transposed version use that
    x = [0,d2];
    y = p2(1)*x+p2(2);
    vector_axon = [diff(x),diff(y)];
    vector_axon = vector_axon/norm(vector_axon);
    
end

% Convention, make 1st element of vector positive
% since direction doesn't matter
if vector_axon(1) < 0
    vector_axon = -vector_axon;
end

x_intercept = p(2);

% Endpoints for plotting during axon grouping
endpoints = [x',y'];
