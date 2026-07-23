function validate_bins(bins, m)
% VALIDATE_BINS Validate a histogram bin specification.
%
% validate_bins(bins, m) throws an error if BINS is not a valid histogram
% specification for M variables. It accepts:
%   * a positive-integer scalar (same number of bins in every dimension),
%   * a positive-integer vector of length M (bins per dimension),
%   * a cell array of length M holding the bin edges of each dimension.
%
% See also histogram_pdf, renyi_entropy_binning, tsallis_entropy_binning

% Copyright (c) 2026, Kristína Mezeiová
% Licensed under the BSD 3-Clause License.

if iscell(bins)
    if numel(bins) ~= m
        error('The number of edge vectors must equal the number of variables.');
    end
    for i = 1:m
        if numel(bins{i}) < 2
            error('Each edge vector must contain at least two elements.');
        end
    end
    return
end

if ~isnumeric(bins) || ~isreal(bins) || any(~isfinite(bins(:)))
    error('Parameter "bins" must be numeric, finite and real, or a cell array of edges.');
end

if isscalar(bins)
    checkPositiveInteger(bins);
    return
end

if isvector(bins)
    if numel(bins) ~= m
        error('The number of bins must be scalar or have one element per variable.');
    end
    for i = 1:numel(bins)
        checkPositiveInteger(bins(i));
    end
    return
end

error('Parameter "bins" must be a scalar, a vector, or a cell array of edges.');
end

function checkPositiveInteger(v)
if ~(v > 0 && v == floor(v))
    error('Parameter "bins" must contain positive integers.');
end
end
