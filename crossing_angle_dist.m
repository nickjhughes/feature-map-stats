
function [angles, counts, crossing_angles] = crossing_angle_dist(od, op, mask, nangles, bins)
%CROSSING_ANGLE_DIST Calculate contour crossing angle distribution.
%
% [angles, counts, crossing_angles] = crossing_angle_dist(od, op, mask, nangles, bins)
%
% Calculates a distribution of intersection angles between the OP and OD maps.
% Uses od_op_crossing to find the crossing angles. A binary mask can be used
% to restrict the analysis. Both the distribution, raw counts in each bin,
% and all the crossing angles are returned.
%
% See also:
% od_op_crossing

% Input validation and defaults
if nargin < 2
    error('At least two inputs are required.');
end
if nargin < 3
    mask = true(size(op));
end
if nargin < 4
    nangles = 4;
end
if nargin < 5
    bins = [0:10:80, 91];
end

% Mask
op(~mask) = nan;
od(~mask) = nan;

% Calculate crossing angles of OP and OD maps
crossing_angles = od_op_crossing(od, op, nangles);

% Calculate distribution
a = histc(crossing_angles*180/pi, bins);
n = a(1:end-1);
counts = n;
angles = n/sum(n);
