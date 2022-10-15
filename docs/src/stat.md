# Statistical tools

## Running mean and variance

There are plenty of implementations of mean and variance of numbers stored in a vector.  However, sometimes  it is more convenient to maintain running estimates of mean and/or variance, adding new data as it pops up, without storing the whole data set in an array.  One way of achieving this is with recursion formulae proposed by D. H. D. West.

### West's recursion formulae

These are recursion relations that allow to hold running estimates of mean and variance without storing the data points.  Additionally, they provide a good numerical estimate of the variance (performing as well as the two-pass algorithm without need for a second pass).

For unweighted data, the formulae are (``x_N`` is the ``N``-th data point, and ``\mu_N`` and ``\sigma^2_N`` are the estimates of mean and variance with ``N`` points, respectively):
```math
\begin{align*}
  \mu_N & = \mu_{N-1} + \frac{1}{N} \left(x_N -\mu_{N-1} \right), \\
  S_N & = S_{N-1} + \frac{N-1}{N} \left(x_N -\mu_{N-1}\right)^2, & \sigma^2_N &= \frac{S_N}{N-1}.
\end{align*}
```

The weighted versions (using ``w_N`` for the weight of ``x_N``) are:
```math
\begin{align*}
   \mu_N & = \mu_{N-1} + \frac{w_N}{\sum_{i=1}^N w_i} \left(x_N -\mu_{N-1} \right), \\
   S_N &= S_{N-1} + \frac{w_N \sum_{i=1}^{N-1} w_i}{\sum_{i=1}^N w_i} \left(x_N -\mu_{N-1}\right)^2, 
   &  \sigma^2_N &= \frac{S_N}{\frac{N-1}{N}\sum_{i=1}^N w_i}.
\end{align*}
```

#### Reference

 - D. H. D. West, Updating mean and variance estimates: an improved method. _Communications of the ACM_ __22,__ 532 (1979).


```@meta
CurrentModule = BioStatPhys
```

### API

#### Unweighted data
```@docs
MeanVar
push!(MV::MeanVar,x)
mean(MV::MeanVar)
var(MV::MeanVar)
```

#### Weighted data
```@docs
WMeanVar
```
- see [`MeanVar`](@ref)
```@docs
push!(MV::WMeanVar,x,w)
mean(MV::WMeanVar)
var(MV::WMeanVar)
```


## Histograms

A type for computation of histograms, with track of outliers.  Provides access to bin counts or probabilities.

### API

```@docs
Histogram
push!(::Histogram,::AbstractFloat)
outliers
area
counts
prob
median
```


## Time correlation functions

The _time (auto-)correlation function_ of a signal ``a(t)`` (which is a noisy quantity that can be found to take different values if it is measured several times, i.e. it is a stochastic process) is defined as
```math
  C(t_0,t) = \langle a^*(t_0) a(t_0+t) \rangle,
```
where ``a^*`` is the complex conjugate and brackets stands for a  history-average (an average over many realizations of the stochastic process, i.e. over many repetitions of the experiment, resetting the initial conditions).

The _connected time correlation_, or _auto-covariance_, is
```math
   C_c(t_0,t) = \left\langle \delta a^*(t_0) \delta a(t_0+t) \right\rangle = 
   C(t_0,t) - \left\langle a^*(t_0) \right\rangle \left\langle a(t_0+t) \right\rangle,
```
where ``\delta a(t)=a(t)-\langle a(t) \rangle``.

Assume the process ``a(t)`` is sampled uniformly in time with interval ``\Delta t``.  Calling ``N`` the number of time samples and ``M``the number of experiments (i.e. different measurements of ``a(t)``  after resetting the initial conditions), we represent these data as ``M`` sequences ``A_i^{(n)}`` with ``i=1,\ldots,N``,  ``n=1,\ldots,M``.  Then the following statistical estimators of the correlation functions can be computed.

In the general case,
```math
\begin{align*}
  \langle a(t_i) \rangle & \approx \overline{A_i} =\frac{1}{M} \sum_{n=1}^M A_i^{(n)}, \\
  C_c(t_i,t_k) & \approx \hat C^{(c)}_{i,k} =\frac{1}{M}\sum_n^M \delta A_i^{(n)}
  \delta A_{i+k}^{(n)}  , \qquad \delta A_i^{(n)} = A_i^{(n)} - \overline{A_i}.
\end{align*}
```

If the process ``a(t)`` is stationary, then the correlation functions time-translation invariant (i.e. independent of ``t_0``), and an estimate can be obtained with a single sample, averaging over the time origin:
```math
\begin{align*}
  \overline{A} &= \frac{1}{N}\sum_{i=1}^N A_i, \\
  \hat C^{(c)}_k &= \frac{1}{N-k} \sum_{j=1}^{N-k} \delta A_j \delta A_{j+k},
  \qquad
  \delta A_j = A_j - \overline{A}.
\end{align*}
```
Note that the estimate obtained this way is noisier the larger the value of ``k``.

If the process is not stationary, then the first estimator, which needs several samples, must be used.  Although we implement this here, note that in the connected case the covariance function `cov` from the `Statistics` package can be used instead and is more convenient.  We also implement the estimator for the stationary case.

For details, and caveats when using these estimators, see the review article

 - T. S. Grigera, Correlation functions as a tool to study collective behaviour phenomena in biological systems. _J. Phys. Complex._ __2,__ 045016 (2021).  [[DOI](https://doi.org/10.1088/2632-072X/ac2b06)]


### API

The stationary (or TTI) estimate for a single real or complex sequence is implemented with an algorithm that uses fast-Fourier-transforms, giving ``O(N\log N)`` performance.

```@docs
time_correlation
```
