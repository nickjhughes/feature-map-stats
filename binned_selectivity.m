
function [mean_sel, std_sel, all_sel] = binned_selectivity(op, od, mask, nbins)
%BINNED_SELECTIVITY Calculate average orientation selectivity in OD bins.
%
% [mean_sel, std_sel, all_sel] = binned_selectivity(op, od, mask, nbins)
%
% Takes the same input as pinwod. Calculates the average OP selectivity in each
% of the OD bins defined in pinwod. Returns the mean and standard deviation of
% each bin, as well as all the selectivity values.
%
% See also:
% pinwod

% Input validation and defaults
if nargin < 2
    error('At least two inputs are required.');
end
if nargin < 3
    mask = true(size(op));
end
if nargin < 4
    nbins = 5;
end

% Convert OP to selectivity
op = abs(op);

% Mask
op(~mask) = nan;
od(~mask) = nan;

% Split the OD map into positive and negative sections
od_pos = od;
od_pos(od <= 0) = nan;
od_neg = od;
od_neg(od > 0) = nan;

% Divide the two OD map parts into nbins quantiles
% Positive
if nbins == 2
    q = quantile(od_pos(:), 2);
    q = q(1);
else
    q = quantile(od_pos(:), nbins-1);
end
bin_od_pos = nan(size(od_pos));
bin_od_pos(od_pos < q(1)) = nbins;
for k = 2:nbins-1
    bin_od_pos(od_pos < q(k) & od_pos > q(k-1)) = nbins - k + 1;
end
bin_od_pos(od_pos > q(nbins-1)) = 1;
% Negative
if nbins == 2
    q = quantile(od_neg(:), 2);
    q = q(1);
else
    q = quantile(od_neg(:), nbins-1);
end
bin_od_neg = nan(size(od_neg));
bin_od_neg(od_neg < q(1)) = nbins;
for k = 2:nbins-1
    bin_od_neg(od_neg < q(k) & od_neg > q(k-1)) = nbins - k + 1;
end
bin_od_neg(od_neg > q(nbins-1)) = 1;

% Calculate average selectivity in OD quantiles
mean_sel = nan(nbins,1);
std_sel = nan(nbins,1);
all_sel = cell(nbins,1);
for j = 1:nbins
    % bin_od_pos => 1 (centre), n (border)
    % bin_od_neg => n (centre), 1 (border)
    mean_sel(j) = mean([op(bin_od_pos == j); op(bin_od_neg == (nbins-j+1))]);
    std_sel(j) = std([op(bin_od_pos == j); op(bin_od_neg == (nbins-j+1))])/sqrt(length([op(bin_od_pos == j); op(bin_od_neg == (nbins-j+1))]));
    all_sel{j} = [op(bin_od_pos == j); op(bin_od_neg == (nbins-j+1))];
end
