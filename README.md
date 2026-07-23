# entropy-estimation-toolbox

A small MATLAB toolbox for estimating **Rényi** and **Tsallis** entropies of
one- and multi-dimensional data. Two families of estimators are provided:

- **Histogram binning** – the probability distribution is estimated from an
  N-dimensional histogram and the entropy is computed from the discrete
  probabilities.
- **k-nearest-neighbours (kNN)** – a non-parametric estimator based on the
  distances to the k-th nearest neighbour. For `q ≠ 1` the
  Leonenko–Pronzato–Savani (LPS) estimator is used, and for `q = 1` the
  Kozachenko–Leonenko (KL) estimator of the Shannon entropy is used as the
  limit case.

All entropy values are returned in **natural units (nats)**.

## Contents

| File | Description |
| --- | --- |
| `histogram_pdf.m` | Estimates an N-dimensional probability distribution using histogram binning. |
| `renyi_entropy_binning.m` | Rényi entropy from a histogram-based probability distribution. |
| `renyi_entropy_knn.m` | Rényi entropy using the kNN (LPS / KL) estimator. |
| `tsallis_entropy_binning.m` | Tsallis entropy from a histogram-based probability distribution. |
| `tsallis_entropy_knn.m` | Tsallis entropy using the kNN (LPS / KL) estimator. |

## Requirements

- **MATLAB** (R2019b or newer recommended – the kNN functions use the
  `arguments` validation block).
- **Statistics and Machine Learning Toolbox** – required by the kNN
  estimators for `knnsearch`.

The histogram functions only use base MATLAB (`discretize`, `accumarray`).

## Installation

Clone the repository and add the folder to your MATLAB path:

```matlab
addpath('path/to/entropy-estimation-toolbox');
```

## Data conventions

- `X` is an `N × m` matrix: **rows are observations, columns are variables**
  (dimensions). A column vector is treated as a 1-D signal.
- The estimators compute the **joint** entropy over all `m` columns, not one
  value per column.
- `q` is the entropy order and can be a scalar or a vector; the entropy is
  evaluated for every value of `q`.

## Function reference

### `histogram_pdf`

```matlab
p            = histogram_pdf(X, bins)
[p, counts]  = histogram_pdf(X, bins)
```

Estimates the probability distribution of `X` by binning.

**Inputs**

- `X` – `N × m` data matrix.
- `bins` – histogram specification:
  - scalar: same number of bins in every dimension,
  - vector of length `m`: number of bins per dimension,
  - cell array of length `m`: explicit bin edges per dimension.

**Outputs**

- `p` – estimated probability of each bin.
- `counts` – observation frequency in each `m`-dimensional bin.

### `renyi_entropy_binning`

```matlab
Hq = renyi_entropy_binning(X, q, bins)
```

Rényi entropy of order(s) `q` from a histogram estimate of the distribution.
Special cases handled explicitly:

- `q → 1` : Shannon entropy `-Σ p·log p`,
- `q = 0` : `log` of the number of occupied bins (Hartley / max-entropy),
- `q = Inf` : min-entropy `-log(max p)`.

### `tsallis_entropy_binning`

```matlab
Sq = tsallis_entropy_binning(X, q, bins)
```

Tsallis entropy of order(s) `q` from a histogram estimate of the distribution.
Special cases:

- `q → 1` : Shannon entropy,
- `q = 0` : number of occupied bins minus one,
- `q = Inf` : `0`.

### `renyi_entropy_knn`

```matlab
H       = renyi_entropy_knn(X, k, q)
H       = renyi_entropy_knn(X, k, q, metric)
[H, N]  = renyi_entropy_knn(__)
```

Rényi entropy estimated from k-nearest-neighbour distances.

**Inputs**

- `X` – `N × m` data.
- `k` – order of the nearest neighbour (positive integer, `k < N`).
- `q` – vector of Rényi orders. Each `q` must satisfy `q < k + 1`.
- `metric` – `"chebychev"` (default) or `"euclidean"`.

**Outputs**

- `H` – Rényi entropy for each value of `q`.
- `N` – number of points actually used (after removing duplicates with zero
  neighbour distance).

### `tsallis_entropy_knn`

```matlab
S       = tsallis_entropy_knn(X, k, q)
S       = tsallis_entropy_knn(X, k, q, metric)
[S, N]  = tsallis_entropy_knn(__)
```

Same interface and constraints as `renyi_entropy_knn`, returning the Tsallis
entropy instead.

## Examples

```matlab
% ----- Histogram-based estimators -----
X = randn(10000, 2);                       % 2-D Gaussian data
q = [0.5 1 2 Inf];

Hq = renyi_entropy_binning(X, q, 50);      % 50 bins per dimension
Sq = tsallis_entropy_binning(X, q, 50);

% ----- kNN-based estimators -----
x = randn(2000, 1);
q = 0.5:0.5:5;

H  = renyi_entropy_knn(x, 15, q);                  % Chebychev metric
[S, N] = tsallis_entropy_knn(x, 15, q, "euclidean");
```

For a 1-D standard normal distribution the (differential) Shannon entropy is
`0.5*log(2*pi*e) ≈ 1.4189` nats, which is a convenient sanity check for the
`q = 1` case of the kNN estimators.

## Notes and limitations

- Histogram estimators are sensitive to the number of bins and to the amount
  of data; too many bins on too little data lead to many empty bins and biased
  estimates.
- kNN estimators require `q < k + 1`; values that violate this constraint are
  dropped and a warning is raised.
- Duplicate observations (zero k-th neighbour distance) are removed before the
  kNN estimation; `N` reports how many points remained.

## References

- A. Rényi, *On measures of entropy and information*, Proc. 4th Berkeley
  Symp. on Math. Statist. and Prob., Vol. 1 (1961), pp. 547–562.
- C. Tsallis, *Possible generalization of Boltzmann–Gibbs statistics*,
  Journal of Statistical Physics 52 (1988), pp. 479–487.
- N. Leonenko, L. Pronzato, V. Savani, *A Class of Rényi Information
  Estimators for Multidimensional Densities*, The Annals of Statistics 36(5)
  (2008), pp. 2153–2182. DOI: 10.1214/07-AOS539.
- L. F. Kozachenko, N. Leonenko, *Sample Estimate of the Entropy of a Random
  Vector*, Probl. Peredachi Inf. 23(2) (1987), pp. 9–16.

## License

BSD 3-Clause License. See [`LICENSE`](LICENSE) for details.

Copyright (c) 2026, Kristína Mezeiová.
