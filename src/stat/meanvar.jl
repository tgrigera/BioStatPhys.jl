###############################################################################
#
# MeanVar
#
"""
    MeanVar

Type to allow running computation of mean and variance (i.e. without
storing the data points).

After the object is created with `MeanVar()`, data are fed one at a
time with `push!`, causing the internal state to be updated but
without storing the new data point (`MeanVar` requires a fixed amount
of storage: two real numbers and an integer).  Mean and variance
estimates from a `MeanVar` object are obtained calling `mean` and
`var`.

`MeanVar` uses West's recursion formula, see

 - D. H. D. West, _Communications of the ACM_ __22__, 532 (1979)

# Example

    mv=MeanVar()
    push!(mv,a)
    push!(mv,b)
    ...
    println("Mean = \$(mean(mv))")
    println("Variance = \$(var(mv))")

# See also

See `WMeanVar` for the case of weighted data.
"""
mutable struct MeanVar
    mean::Float64
    pvar::Float64
    N::Int
end

MeanVar()=MeanVar(0.,0.,0)

"""
    function push!(MV::MeanVar,x)

Add data point `x` to `MeanVar` object `MV`
"""
function Base.push!(MV::MeanVar,x)
    MV.N+=1
    Q::Float64=x-MV.mean
    R=Q/MV.N
    MV.mean+=R
    MV.pvar+=Q*R*(MV.N-1)
    return MV
end

"Compute mean estimate (sample mean) from a `MeanVar` object"
@inline mean(MV::MeanVar) = MV.mean
"Compute variance estimate (population variance) from a `MeanVar` object"
@inline var(MV::MeanVar) = MV.pvar/(MV.N-1)

function Base.show(io::IO,mv::MeanVar)
    println(io,"MeanVar object with $(mv.N) datapoints")
end

###############################################################################
#
# WMeanVar
#
"""
    WMeanVar

Type to allow running computation of mean and variance of weighted
data (without storing the data points).

The object is created with `WMeanVar()`, data are fed one at a time
with `push!`, and mean and variance are obtained calling `mean` and
`var`.

Interface is the same as for the unweighted version `MeanVar', which
see for an example.

`WMeanVar` uses West's recursion formula, see

 - D. H. D. West, _Communications of the ACM_ __22__, 532 (1979)
"""
mutable struct WMeanVar
    mean::Float64
    pvar::Float64
    sumW::Float64
    N::Int
end

WMeanVar()=WMeanVar(0.,0.,0.,0)

function Base.push!(MV::WMeanVar,x,w)
    nsW=MV.sumW+w
    Q::Float64=x-MV.mean
    R=Q*w/nsW
    MV.mean+=R
    MV.pvar+=Q*R*MV.sumW
    MV.sumW=nsW
    MV.N+=1
    return MV
end

@inline mean(MV::WMeanVar) = MV.mean
@inline var(MV::WMeanVar) = MV.pvar*MV.N/(MV.sumW*(MV.N-1))

function Base.show(io::IO,wmv::WMeanVar)
    println(io,"WMeanVar object with $(wmv.N) datapoints")
end
