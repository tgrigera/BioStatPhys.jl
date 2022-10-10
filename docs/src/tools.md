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


### API

```@docs
BinnedVector
interval
delta
nbins
bin
binc
```
