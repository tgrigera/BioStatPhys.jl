# Miscellaneous tools

## Real-indexed vector

This is a simple implementation of a vector that can be accessed with real (floating point) sub-indices.  It does simply a linear binning of a specified (semi-open) interval ``[min,max)`` on the real line, and provides a convenient `Vector`-like syntax to access elements.  The stored values can be of any type, but the interval and number of bins are fixed at construction.

```@repl 1
    push!(LOAD_PATH,homedir()*"/software/BioStatPhys.jl"); # hide
    using BioStatPhys; # hide
	A = BinnedVector{String}(5,min=0.,max=1.)
	A[0.1]="hello"
	A[0.2]="bye bye"
	A
```
	
To initialise elements with a default value, give function with a signature like `zeros`: 
```@example 1
B = BinnedVector{Rational{Int64}}(10,min=2.,max=8.,init=zeros)
B[2.1]=3//5
B[5.3]=4//1
B
```

Alternatively, it is possible to construct specifying bin width and the interval.  In this case, one of the values will need to be rounded to get an integer number of bins:
```@repl 1
	C = BinnedVector{Int}(Δ=0.15,min=5.,max=10.,round_Δ=RoundUp); (C.Δ, interval(C))
	D = BinnedVector{Int}(Δ=0.15,min=5.,max=10.,round_max=RoundUp); (D.Δ, interval(D))
	E = BinnedVector{Int}(Δ=0.15,min=5.,max=10.,round_min=RoundUp); (E.Δ, interval(E))
```
`RoundDown` is also recognised.

Indexing with integers refers to bin numbers.
```@repl 1
	A[1:2]
	B[1:5]
```

Outliers are mapped to bin 0 if below range, or -1 if above range.
```@repl 1
B[0.]=10//3 ; B[40.]=99//2;
B[0]
B[-1]
```

As an example of use, the [`Histogram`](@ref) type is implemented using `BinnedVector`.


### Special treatment of zero

In some cases, like when computing space correlation functions, it is convenient to treat the case `A[0.]` separately from the bin that starts at 0.  This is what the `ZBinnedvector` type does.

```@repl 1
Z = ZBinnedVector{Int}(Δ=0.2,max=0.5,round_max=RoundUp,init=zeros)
Z[0.] += 10
Z[0.1] += 10
Z[0.15] += 10
Z
```

Note that there is no `min` argument in the construction of `ZBinnedvector`, since it is always zero.  In other respects, `ZBinnedvector` behaves like `BinnedVector`, except that bin centres are not equally spaced: the first (pesudo-bin) is at 0, the next bin centre is `delta(Z)/2`, and the rest are `delta(Z)` units apart.


### API

```@docs
BinnedVector
ZBinnedVector
interval
delta
nbins
bin
binc
Base.range(::BinnedVector)
```


## Distance binning

The `distance_binning` function takes a series of points in (2- or 3-``d``) space and creates a `BinnedVector` that classifies the pairs according to their Euclidean distance.

```@docs
DistanceBinning
distance_binning
```


## Regions

Types and functions to describe periodic and non-periodic regions of space.

```@docs
Region
Rectangle
PeriodicRectangle
dimension
volume
dborder
distancesq
```
