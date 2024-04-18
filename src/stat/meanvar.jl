# MeanVar.jl -- Running mean and variance with West's recursion
#
# Copyright (C) 2022 by Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# For details see the file LICENSE in the root directory, or check
# <https://www.gnu.org/licenses/>.

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
    function push!(MV::MeanVar,X::AbstractVector)

Add data point `x` to `MeanVar` object `MV`.  If `X` is a vector,
iteratively `push!` all its elements.
"""
function Base.push!(MV::MeanVar,x)
    MV.N+=1
    Q::Float64=x-MV.mean
    R=Q/MV.N
    MV.mean+=R
    MV.pvar+=Q*R*(MV.N-1)
    return MV
end

Base.push!(MV::MeanVar,X::AbstractVector) = for x ∈ X push!(MV,x) end

"Compute mean estimate (sample mean) from a `MeanVar` object"
@inline mean(MV::MeanVar) = MV.mean
"Compute variance estimate (population variance) from a `MeanVar` object"
@inline var(MV::MeanVar) = MV.pvar/(MV.N-1)

"""
    empty!(mv)

Clear all data from `MeanVar` or `WMeanVar` object `mv`
"""
@inline function Base.empty!(mv::MeanVar)
    mv.mean = 0.
    mv.pvar = 0.
    mv.N = 0
end


"""
    isempty(mv::MeanVar) -> Bool

Determine whether the `MeanVar` object is empty (has not been fed any
data via `push!`
"""
@inline Base.isempty(mv::MeanVar) = mv.N==0

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

@inline function Base.empty!(mv::WMeanVar)
    mv.mean = 0.
    mv.pvar = 0.
    mv.sumW = 0.
    mv.N = 0
end

"""
    function push!(WMV::WMeanVar,x)

Add data point `x` with weight `W` to `WMeanVar` object `WMV`
"""
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

"Compute mean estimate (sample mean) from a `WMeanVar` object"
@inline mean(MV::WMeanVar) = MV.mean
"Compute variance estimate (population variance) from a `WMeanVar` object"
@inline var(MV::WMeanVar) = MV.pvar*MV.N/(MV.sumW*(MV.N-1))

function Base.show(io::IO,wmv::WMeanVar)
    println(io,"WMeanVar object with $(wmv.N) datapoints")
end
