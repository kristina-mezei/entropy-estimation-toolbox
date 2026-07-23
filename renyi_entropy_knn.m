function [H, N] = renyi_entropy_knn(X, k, q, metric)
% RENYI_ENTROPY_KNN estimates Rényi's entropy using knn-algorithm (LPS)
% Syntax
%
% H = renyi_entropy_knn(X,k,q)
% H = renyi_entropy_knn(X,k,q,metric)
% [H,N] = renyi_entropy_knn(__)
% 
% Inputs
%  X - Data containing N_all observations in m-dimensional space.
%       vector | matrix
%  k - order of the nearest neighbours
%       positive integer
%  q - Rényi parameters
%       vector
%  metric - distance metrics to be used for calculation, can be chebychev or euclidean,
%           default value is chebychev
% 
% Outputs
%   H - Rényi's entropy calculated for each value of q
%   N - number of points actualy used for calculation
% 
% Examples
% X = randn(1000,1);
% H = renyi_entropy_knn(X, 15, [0.5:0.5:5]);
% 
% [H, N] = renyi_entropy_knn(X, 15, [0.5:0.5:5], "euclidean");
% 
% References
% Leonenko N., Pronzato L., Savani V. (2008).
% A Class of Rényi Information Estimators for Multidimensional Densities.
% The Annals of Statistics, 36(5), 2153–2182.
% DOI: 10.1214/07-AOS539
% 
% Kozachenko L. F., Leonenko N. (1987).
% Sample Estimate of the Entropy of a Random Vector
% Probl. Peredachi Inf., 23(2), 9–16. 
% 
% 
% See also
% histogram_pdf
% renyi_entropy_binning
% renyi_entropy_knn
% tsallis_entropy_binning
% tsallis_entropy_knn
%
% for q != 1 function calculates Leonenko, Prozanto, Savani (LPS) algorithm for Rényi entropy
% for q = 1 function calculates Kozachenko - Leonenko estimate for Shannon entropy
%




% Input checking

arguments
    X  (:,:) double {mustBeNonempty}
    k (1,1) double {mustBeNonempty, mustBePositive, mustBeInteger}
    q  (1,:) double {mustBeNonempty, mustBeReal}
    metric string {mustBeMember(metric, ["chebychev","euclidean"])} = 'chebychev'
end


% N_all is the number of observations, m is the dimension of the space
[N_all,m] = size(X);

% nq is number of Rényi parameters
nq = length(q);


% Preallocation of memory
H = NaN(1, nq);


% check of parameters
if any(k >= N_all)
    error('k must be smaller than the number of points!');
end


if any(k + 1 - q <= 0)
    warning('Parameter q is too large for given k (k+1-q must be > 0).');
end





% Selection of appropiate values of q, so that q < k+1
valid_q_ind = q < (k + 1);
valid_q = q(valid_q_ind);

if(isempty(valid_q))
    error("Parameter q must be less than k+1!")
end



% Selection of indices for which SPL or KL algorithms are calculated
% q≈1 is evaluated using the Kozachenko–Leonenko estimator
% to avoid numerical instability of the Rényi formula.
tol = 1e-8;

indKL = abs(valid_q-1) < tol;
indLPS  = abs(valid_q-1) >= tol;




% Distances to knn:
[~, dist] = knnsearch(X, X, 'K', k+1, "Distance", lower(metric));
% rho_k - distances to the kth nearest neighbours:
rho_k = dist(:,k+1);   % k-th neightbour, 1st is the point itself

% check for duplicate points
nDuplicates = sum(rho_k == 0);

if nDuplicates > 0
    warning('%d observations with zero kNN distance were removed.', nDuplicates);
    rho_k = rho_k(rho_k > 0);
end

N = numel(rho_k);


if N < 2 || N <= k 
    error('Too few non-duplicate observations remain.');
end

log_rho = log(rho_k);


% log_norm
log_norm = log(N-1);

% Volume of unit ball in m-dimensional space
% Vm differs for various metrics
% log of the volume of the unite ball is calculated

switch lower(metric)
    case "chebychev"
        logVm = m*log(2);
    case "euclidean"
        logVm = (1/2)*m*log(pi) - gammaln(1 + m/2);
end


% LPS algorithm
if any(indLPS)

    % Iq calculation:
    % zeta is calculated with the use of logarithm
    alpha = reshape(1-valid_q(indLPS),1,[]);
    % log_zeta = (1-valid_q(indSPL)).*log_norm + gammaln(k) - gammaln(k+1-valid_q(indSPL)) + (1-valid_q(indSPL)).*logVm + m*(1-valid_q(indSPL)).*log_rho;
    log_zeta = alpha.*log_norm + gammaln(k) - gammaln(k+alpha) + alpha.*logVm + m*alpha.*log_rho;

    % The computation of the Rényi functional was implemented using a numerically stable log-sum-exp transformation to avoid overflow in large-sample regimes.
    a = max(log_zeta);              % reference shift
    log_mean_exp = a + log(mean(exp(log_zeta - a), 1, 'omitnan'));

    % Rényi entropy
    H(indLPS) = log_mean_exp ./ (1 - valid_q(indLPS));

end

% KL algorithm
if any(indKL)
    H(indKL) = psi(N) - psi(k) + logVm + m*mean(log_rho, 1, 'omitnan');
    % Hy= - psi(k) + psi(N) + M*mean(log(2*dd), 'omitnan');
end