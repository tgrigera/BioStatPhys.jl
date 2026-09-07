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


## Mean and variance in growing windows

`GeoAve` is a type to average time series in geometrically growing
windows.  The `push!` function takes pairs of numbers ``x,y`` and
averages together all points whose ``x``-coordinate (assume it is a
time) falls within a window.  The first window starts at ``t_0`` and
is of length ``b``, successive windows grow geometrically by factor
``w``.  In other words, window boundaries are located at
```math
 t_n = t_0 + b \frac{w^n - 1}{w-1} 
```
When returning the averages (`get_mean()`), the  time assigned to the
average value is the center of the window, i.e. the  abscissas `s_n` are
```math
      s_n = t_n + \frac{1}{2} b w^n 
```

Points may be given in any order, and averages are available at any
time.  Points can be added after querying for the accumulated mean and
variance.  The initial time `t0` is handled separately, and values
corresponding to `t0` are averaged without windowing.

The case ``w=1`` is supported and handled as a special case, yielding
windows of fixed width equal to ``b``. ``w<1`` is not recommended.

For example:
```julia
G = GeoAve(t0=0,wfactor=1.2,base=0.5)
for i=1:100 push!(G,i*rand(),rand()) end
t,m,v = get_mean(G)
```

### API

```@docs
GeoAve
get_mean(::GeoAve)
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


## Mann-Whitney test

The Mann-Whitney test is designed to check for an "effect" of the
variation of some experimental condition on the experimental result.
This test assumes a real variable is measured under "control" and
"experimental" conditions several times in two groups of otherwise identical
experiments.  The aim is to establish to what degree the variable
measured under experimental conditions is likely to exceed the value
measured under control conditions.

The version of the test implemented here is Bayesian, basically meaning that one asks how the probability distribution of the quantity of interest changes after we know the experimental result.  The quantity of interest of here is ``\Theta``, namely the probability that a random experimental observation ``X_E`` will be larger than a random control observation ``X_C``,
```math
    \Theta = P(X_E > X_C).
```
If ``\Theta>1/2`` then the experimental conditions are likely to
yield higher values of ``X`` than the control conditions.

One starts from a _prior_ probability distribution, which encodes what one knows or expects about the quantity before the experiment; here the user can give a prior for ``\Theta``, or by default a flat prior (i.e. a uniform distribution over the interval ``[0,1]``) will be taken, expressing the fact that we know nothing before the experiment is done.  The experimental data then allows computation of the _posterior_ probability ``P(\Theta | U_E)``, i.e. the distribution of ``\Theta``, given the value ``U_E`` (the number of times the experiment exceeds control) determined in the experiment.  The method also needs the _likelihood_ ``P(U_E|\Theta)``, i.e. the probability of measuring ``U_E`` for a
given value of ``\Theta``.  This is computed either through a Monte
Carlo or with an analytical approximation, as discussed by Chechile (2020).

The output of the test is the posterior ``P(\Theta | U_E)``, as well as the median and the confidence or credibility interval (the interval encompassing 95% of the probability around the median).

Example usage:

```julia
using Distributions
using GLMakie
using BioStatPhys
using Random

rng = Random.Xoshiro(260833)

XC = rand(rng,12)   # Results of the "control" experiment
XE = rand(rng,7)    # Results of the "experimental" experiment

r = BioStatPhys.Bayesian_U_test(XC,XE,rng=rng,method=:MC)

fig = Figure()
ax = Axis(fig[1,1], xlabel="Θ",ylabel="P(Θ|U)")
lines!(ax,r.Θbins,r.PΘ_U./(r.Θbins[2]-r.Θbins[1]))
display(fig)
```


#### Reference

- Richard A. Chechile, _Bayesian Statistics for Experimental
   Scientists,_ MIT Press, Cambdrige, MA (2020), Chapter 7.

### API

```@docs
Bayesian_U_test
```
