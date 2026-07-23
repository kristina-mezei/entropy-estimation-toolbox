function [p,counts] = histogram_pdf(X, bins)
% HISTOGRAM_PDF Estimate a probability distribution using histogram binning.
% 
% Syntax
% p = histogram_pdf(X,bins)
% [p,counts] = histogram_pdf(X,bins)
% 
% Inputs
%   X - N-by-m matrix, where rows are observations and columns are
%          variables.
%
%   bins - Histogram specification:
%          * scalar: same number of bins in every dimension.
%          * vector of length m: number of bins for each dimension.
%          * cell array of length m: bin edges for each dimension.
% 
% Outputs
%   p      - probability estimated by dividing the number of elements in each bin by the total number of elements in the input data.
%   counts - frequency of observations of the data in m-dimensional bins
% 
% Examples
%   % Same number of bins in all dimensions
%   X = randn(1000, 2);
%   [p, counts] = histogram_pdf(X, 50);
%
%   % Different number of bins in each dimension
%   [p, counts] = histogram_pdf(X, [40 60 30]);
%
%   % User-defined bin edges
%   edges = {
%       linspace(-3,3,51), ...
%       linspace(0,10,31), ...
%       [-5 -2 -1 0 1 2 5]
%   };
%   [p, counts] = histogram_pdf(data, edges);
% 
% References
% 
% See also
% histogram_pdf
% renyi_entropy_binning
% renyi_entropy_knn
% tsallis_entropy_binning
% tsallis_entropy_knn


% Copyright (c) 2026, Kristína Mezeiová
% Licensed under the BSD 3-Clause License.
% 
%% Get dimensions of the signal
[N, m] = size(X);

%% Construct bin edges

if iscell(bins)

    if numel(bins) ~= m
        error('The number of edge vectors must equal the number of variables.');
    end

    edges = bins;
    sz = zeros(1,m);

    for i = 1:m
        if numel(edges{i}) < 2
            error('Each edge vector must contain at least two elements.');
        end

        sz(i) = numel(edges{i}) - 1;
    end

else

    if isscalar(bins)
        bins = repmat(bins,1,m);
    elseif numel(bins) ~= m
        error('The number of bins must be scalar or have one element per variable.');
    end

    edges = cell(1,m);

    for i = 1:m

        xmin = min(X(:,i));
        xmax = max(X(:,i));

        % Prevent zero-width interval
        if xmin == xmax
            xmax = xmax + eps(xmax);
        end

        edges{i} = linspace(xmin,xmax,bins(i)+1);

    end

    sz = bins;

end
%% Assign samples to histogram bins

subs = zeros(N,m);

for i = 1:m
    subs(:,i) = discretize(X(:,i),edges{i});
end

%% Remove samples outside the specified edges

valid = all(subs>0,2);
subs = subs(valid,:);
%% N-dimensional histogram

%
if(m == 1)
    counts = accumarray(subs,1,[sz, 1]);
elseif(m>1)
    counts = accumarray(subs,1,sz);
end


%% Probability distribution
p = counts/N;
end


