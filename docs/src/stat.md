# General statistical tools

## Running mean and variance

There are plenty of implementations of mean and variance of numbers stored in a vector.  However, sometimes  it is more convenient to maintain running estimates of mean and/or variance, adding new data as it pops up, without storing the whole dataset in an array.  One way of achieving this is with West's recursion formula, which we implement here.

```@meta
CurrentModule = BioStatPhys
```

### Unweighted data
```@docs
MeanVar
push!(MV::MeanVar,x)
mean(MV::MeanVar)
var(MV::MeanVar)
```

### Weighted data
```@docs
WMeanVar
```
- see [`MeanVar`](@ref)
```@docs
push!(MV::WMeanVar,x,w)
mean(MV::WMeanVar)
var(MV::WMeanVar)
```
