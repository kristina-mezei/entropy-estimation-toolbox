function Hq = renyi_entropy_binning(X, q, bins)

% RENYI_ENTROPY_BINNING estimates Rényi's entropy using histogram binning.
% Reference:
% A. Rényi: On measures of entropy and information, in Proceedings of the Fourth Berkeley Symposium on Mathematical
% Statistics and Probability, Volume 1: Contributions to the Theory of Statistics,
% Vol. 4 (University of California Press, 1961) pp. 547–562.


% Entropy values are returned in natural units (nats).
%
%
%  Inputs:
%   X - N-by-m matrix, where rows are observations and columns are
%          variables.
%
%    q    : Rényi entropy order(s). Can be a scalar or a vector.
%
%
%   bins - Histogram specification:
%          * scalar: same number of bins in every dimension.
%          * vector of length m: number of bins for each dimension.
%          * cell array of length m: bin edges for each dimension.
%
%
%
%  Outputs:
%    Hq            : Rényi's entropy values. Rows correspond to columns of
%                   data, columns correspond to q values.
%



% Input checking

% arguments
%     X  (:,:) double {mustBeNonempty}
%     q  (1,:) double {mustBeNonempty}
% end

%% Get dimensions of the signal
[N, m] = size(X);

if isscalar(bins)
    if ~(isnumeric(bins) && isscalar(bins) && isfinite(bins) && bins > 0 && bins == floor(bins))
        error('Parameter "bins" must be a positive integer scalar.');
    end
end
if isvector(bins) && any(size(bins) - [1 1])
    bins = bins(:);
    if numel(bins) ~= m
        error('The number of bins must be scalar or have one element per variable.');
    end
    for i = 1:m
        if ~(isnumeric(bins(i)) && isscalar(bins(i)) && isfinite(bins(i)) && bins(i) > 0 && bins(i) == floor(bins(i)))
            error('Parameter "bins" must be a positive integer scalar.');
        end
    end

end

if iscell(bins)

    if numel(bins) ~= m
        error('The number of edge vectors must equal the number of variables.');
    end


    for i = 1:m
        if numel(bins{i}) < 2
            error('Each bins vector must contain at least two elements.');
        end
    end
end






% Initialization of variables
Hq = NaN(length(q),1);



%  Calculation of probability density function
p = histogram_pdf(X, bins);
p_nonzero = p(p>0);

% Skip the calculation for a wrong measurement (all values are NaN)
if isempty(p_nonzero)
    error('Wrong measurement, all values are NaN. Probability function is zero.')
end

% Calculation of Renyi entropy:
for j = 1: length(q)

    % Limit case when q = 1 is the Shannon entropy:
    if abs(q(j)-1) <  1e-10
        Hq(j) = -sum(p_nonzero .* log(p_nonzero));

        % Min-entropy (the most probable state)
    elseif isinf(q(j))

        Hq(j) = -log(max(p_nonzero));

        % For q = 0 Hq counts the nonzero bins
    elseif q(j)==0
        Hq(j) = log(numel(p_nonzero));

    else
        % Rényi's entropy
        Hq(j) = log(sum(p_nonzero.^q(j)))/(1 - q(j));
    end
end