function Sq = tsallis_entropy_binning(X, q, bins)

% TSALLIS_ENTROPY_BINNING estimates Tsallis' entropy using histogram
% binning.
%
% Reference:
% C. Tsallis: Possible generalization of Boltzmann-Gibbs statistics. Jour-
% nal of statistical physics 52, 1988, pp. 479–487.
%


% Entropy values are returned in natural units (nats).
%
%
%  Inputs:
%   X - N-by-m matrix, where rows are observations and columns are
%          variables.
%
%    q    : entropy order(s). Can be a scalar or a vector.
%
%
%   bins - Histogram specification:
%          * scalar: same number of bins in every dimension.
%          * vector of length m: number of bins for each dimension.
%          * cell array of length m: bin edges for each dimension.
%
%
%  Outputs:
%    Sq    : Tsallis entropy values. 

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
Sq = NaN(1, length(q));



%  Calculation of probability density function
p = histogram_pdf(X, bins);
p_nonzero = p(p>0);

% Skip the calculation for a wrong measurement (all values are NaN)
if isempty(p_nonzero)
    error('Wrong measurement, all values are NaN. Probability function is zero.')
end


% Calculation:
for j = 1: length(q)

    % Limit case when q = 1 is the Shannon entropy:
    if abs(q(j)-1) <  1e-10
        Sq(j) = -sum(p_nonzero .* log(p_nonzero));

        % Min-entropy (the most probable state)
    elseif isinf(q(j))

        Sq(j) = 0;

        % For q = 0, Tsallis entropy equals the number of occupied bins minus one.
    elseif q(j)==0
        Sq(j) = numel(p_nonzero)-1;

    else
        % Tsallis' entropy
        Sq(j) = (1-sum(p_nonzero.^q(j)))/(q(j)-1);
    end
end