function tests = entropyToolboxTest
% ENTROPYTOOLBOXTEST Unit tests for the entropy-estimation-toolbox.
%
% Run all tests from the repository root (or from anywhere) with:
%   results = runtests('tests/entropyToolboxTest.m')
% or
%   results = runtests('entropyToolboxTest')   % if tests/ is on the path
%
% The tests cover:
%   * histogram_pdf normalisation and counts,
%   * analytic uniform-distribution values for the binning estimators,
%   * the special cases q = 1 (Shannon), q = 0 and q = Inf,
%   * the differential Shannon entropy of a Gaussian via the kNN estimators,
%   * a regression test for the q-alignment fix in the kNN estimators,
%   * output shapes.

tests = functiontests(localfunctions);
end

%% Fixtures ---------------------------------------------------------------

function setupOnce(testCase)
% Make the toolbox functions available regardless of the current folder.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
end

%% Helpers ----------------------------------------------------------------

function [X, B, r] = uniformData()
% Data that fills B equal-width bins with exactly r observations each, so
% the estimated distribution is exactly uniform (p = 1/B in every bin).
B = 8;
r = 100;
centres = (0.5:1:(B-0.5))';          % 0.5, 1.5, ..., B-0.5
X = repmat(centres, r, 1);
end

function edges = unitEdges(B)
edges = {0:B};                         % B unit-width bins [0,1),...,[B-1,B)
end

%% histogram_pdf ----------------------------------------------------------

function testHistogramNormalisation(testCase)
X = randn(5000, 1);
edges = {linspace(-6, 6, 41)};         % wide enough to contain all samples
[p, counts] = histogram_pdf(X, edges);
verifyEqual(testCase, sum(p(:)), 1, 'AbsTol', 1e-12);
verifyEqual(testCase, sum(counts(:)), numel(X));
end

function testHistogramUniformCounts(testCase)
[X, B, r] = uniformData();
[p, counts] = histogram_pdf(X, unitEdges(B));
verifyEqual(testCase, counts(:), repmat(r, B, 1));
verifyEqual(testCase, p(:), repmat(1/B, B, 1), 'AbsTol', 1e-12);
end

function testHistogram2D(testCase)
X = [ (0.5:1:3.5)', (0.5:1:3.5)' ];    % 4 points on the diagonal
[~, counts] = histogram_pdf(X, {0:4, 0:4});
verifyEqual(testCase, size(counts), [4 4]);
verifyEqual(testCase, trace(counts), 4);   % all mass on the diagonal
verifyEqual(testCase, sum(counts(:)), 4);
end

%% renyi_entropy_binning --------------------------------------------------

function testRenyiBinningUniform(testCase)
% For a uniform distribution over B states every Rényi order gives log(B).
[X, B, ~] = uniformData();
q = [0 0.5 1 2 5 Inf];
Hq = renyi_entropy_binning(X, q, unitEdges(B));
verifyEqual(testCase, Hq, repmat(log(B), 1, numel(q)), 'AbsTol', 1e-10);
end

function testRenyiBinningShannonLimit(testCase)
% q = 1 must reproduce the plug-in Shannon entropy.
X = randn(4000, 1);
edges = {linspace(-6, 6, 33)};
p = histogram_pdf(X, edges);
pnz = p(p > 0);
shannon = -sum(pnz .* log(pnz));
Hq = renyi_entropy_binning(X, 1, edges);
verifyEqual(testCase, Hq, shannon, 'AbsTol', 1e-10);
end

function testRenyiBinningShape(testCase)
[X, B, ~] = uniformData();
q = [0.5 1 2 3];
Hq = renyi_entropy_binning(X, q, unitEdges(B));
verifyEqual(testCase, size(Hq), [1 numel(q)]);
end

%% tsallis_entropy_binning ------------------------------------------------

function testTsallisBinningUniform(testCase)
% Uniform over B states: S_q = (1 - B^(1-q)) / (q - 1); S_1 = log(B);
% S_0 = B - 1; S_Inf = 0 (special case implemented in the function).
[X, B, ~] = uniformData();
q = [0.5 2 5];
Sq = tsallis_entropy_binning(X, q, unitEdges(B));
expected = (1 - B.^(1 - q)) ./ (q - 1);
verifyEqual(testCase, Sq, expected, 'AbsTol', 1e-10);

verifyEqual(testCase, tsallis_entropy_binning(X, 1, unitEdges(B)), log(B), 'AbsTol', 1e-10);
verifyEqual(testCase, tsallis_entropy_binning(X, 0, unitEdges(B)), B - 1, 'AbsTol', 1e-10);
verifyEqual(testCase, tsallis_entropy_binning(X, Inf, unitEdges(B)), 0, 'AbsTol', 1e-12);
end

function testTsallisBinningShape(testCase)
[X, B, ~] = uniformData();
q = [0.5 1 2 3];
Sq = tsallis_entropy_binning(X, q, unitEdges(B));
verifyEqual(testCase, size(Sq), [1 numel(q)]);
end

%% renyi_entropy_knn ------------------------------------------------------

function testRenyiKnnGaussianShannon(testCase)
% Differential Shannon entropy of a standard normal is 0.5*log(2*pi*e).
rng(42);
X = randn(8000, 1);
H = renyi_entropy_knn(X, 10, 1);
verifyEqual(testCase, H, 0.5*log(2*pi*exp(1)), 'AbsTol', 0.05);
end

function testRenyiKnnShape(testCase)
rng(7);
X = randn(1000, 1);
q = 0.5:0.5:3;
H = renyi_entropy_knn(X, 10, q);
verifyEqual(testCase, size(H), [1 numel(q)]);
verifyTrue(testCase, all(isfinite(H)));
end

function testRenyiKnnQAlignment(testCase)
% Regression test for the q-alignment fix: an invalid q (>= k+1) placed in
% the MIDDLE of the vector must yield NaN at exactly that position, while
% the surrounding valid orders remain finite and correctly placed.
rng(1);
X = randn(1500, 1);
k = 5;                                  % valid requires q < k+1 = 6
q = [0.5 10 2];                         % q(2) = 10 is invalid
H = suppressExpectedWarning(@() renyi_entropy_knn(X, k, q));

verifyTrue(testCase, isfinite(H(1)), 'H(1) should be finite (q=0.5).');
verifyTrue(testCase, isnan(H(2)),    'H(2) must be NaN for the invalid q=10.');
verifyTrue(testCase, isfinite(H(3)), 'H(3) should be finite (q=2).');

% The valid entries must equal the result of calling with the valid q only.
Href = renyi_entropy_knn(X, k, [0.5 2]);
verifyEqual(testCase, [H(1) H(3)], Href, 'AbsTol', 1e-12);
end

function testRenyiKnnQAlignmentWithShannon(testCase)
% Same as above but the surviving order includes the Shannon case q = 1,
% which is handled by a different (KL) branch.
rng(2);
X = randn(1500, 1);
k = 5;
q = [0.5 10 1];
H = suppressExpectedWarning(@() renyi_entropy_knn(X, k, q));
verifyTrue(testCase, isfinite(H(1)));
verifyTrue(testCase, isnan(H(2)));
verifyTrue(testCase, isfinite(H(3)));
end

%% tsallis_entropy_knn ----------------------------------------------------

function testTsallisKnnShape(testCase)
rng(9);
X = randn(1000, 1);
q = 0.5:0.5:3;
S = tsallis_entropy_knn(X, 10, q);
verifyEqual(testCase, size(S), [1 numel(q)]);
verifyTrue(testCase, all(isfinite(S)));
end

function testTsallisKnnQAlignment(testCase)
rng(3);
X = randn(1500, 1);
k = 5;
q = [0.5 10 2];
S = suppressExpectedWarning(@() tsallis_entropy_knn(X, k, q));
verifyTrue(testCase, isfinite(S(1)));
verifyTrue(testCase, isnan(S(2)));
verifyTrue(testCase, isfinite(S(3)));

Sref = tsallis_entropy_knn(X, k, [0.5 2]);
verifyEqual(testCase, [S(1) S(3)], Sref, 'AbsTol', 1e-12);
end

%% Error handling ---------------------------------------------------------

function testKnnAllInvalidQErrors(testCase)
X = randn(200, 1);
verifyError(testCase, @() suppressExpectedWarning(@() renyi_entropy_knn(X, 3, [10 20])), ...
    ?MException);
end

function testBinsValidation(testCase)
X = randn(100, 2);
verifyError(testCase, @() histogram_pdf(X, [4 5 6]), ?MException);  % wrong length
verifyError(testCase, @() histogram_pdf(X, -3), ?MException);       % not positive
verifyError(testCase, @() histogram_pdf(X, 2.5), ?MException);      % not integer
end

%% Local utilities --------------------------------------------------------

function out = suppressExpectedWarning(fcn)
% Run FCN while silencing the expected "q too large" warning so that it does
% not clutter the test output.
w = warning('off', 'all');
cleanup = onCleanup(@() warning(w));
out = fcn();
end
