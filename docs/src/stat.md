# Statistical tools

## Running mean and variance

There are plenty of implementations of mean and variance of numbers stored in a vector.  However, sometimes  it is more convenient to maintain running estimates of mean and/or variance, adding new data as it pops up, without storing the whole data set in an array.  One way of achieving this is with recursion formulae proposed by D. H. D. West.

### West's recursion formulae

These are recursion relations that allow to hold running estimates of mean and variance without storing the data points.  Additionally, they provide a good numerical estimate of the variance (performing as well as the two-pass algorithm without need for a second pass).

For unweighted data, the formulae are:
```math
\overline{X}_N = \overline{X}_{N-1} + \frac{1}{N} \left(x_N -\overline{X}_{N-1} \right),
```
for the mean ``\overline{X}_N``, where ``x_N`` is the ``N``-th data point, and
```math
T_N = T_{N-1} + \frac{N-1}{N} \left(x_N -\overline{X}_{N-1}\right)^2, \qquad S^2_N = \frac{T_N}{N-1},
```
for the variance ``S^2_N``.

The weighted versions (using ``w_N`` for the weight of ``x_N``) are:
```math
 \overline{X}_N = \overline{X_{N-1}} + \frac{w_N}{\sum_{i=1}^N w_i} \left(x_N -\overline{X}_{N-1} \right),
```
```math
 T_N = T_{N-1} + \frac{w_N \sum_{i=1}^{N-1} w_i}{\sum_{i=1}^N w_i} \left(x_N -\overline{X}_{N-1}\right)^2, \qquad S^2_N = \frac{T_N}{\frac{N-1}{N}\sum_{i=1}^N w_i},
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

