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

arguments
    X  (:,:) double {mustBeNonempty}
    q  (1,:) double {mustBeNonempty, mustBeReal}
    bins
end

%% Get dimensions of the signal
[~, m] = size(X);

% Validate the bin specification
validate_bins(bins, m);

% Initialization of variables (row vector: one entropy value per q)
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