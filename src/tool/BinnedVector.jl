# BinnedVector.jl
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

"""
    mutable struct BinnedVector{T}

Simple real-indexed vector.  Maps a real interval to an integer range
`1:nbins` then used to index a vector of arbitrary type.  Can be used
e.g.  to build histograms.  The real interval and number of bins
`nbins` are fixed at the outset.

# Example

Create with

    A = BinnedVector{Int}(nbins,min=1.,max=10.,init=zeros)

for fixed number of bins, or

    A = BinnedVector{Int}(;Δ=0.1,min=1.,max=10.init=zeros,round_Δ=RoundUp)
    A = BinnedVector{Int}(;Δ=0.1,min=1.,max=10.init=zeros,round_min=RoundUp)
    A = BinnedVector{Int}(;Δ=0.1,min=1.,max=10.init=zeros,round_max=RoundUp)

to get fixed `(min, max)`, `(Δ,max)` or `(min,Δ)` respectively
(rounding mode `RoundDown` is also recognised).

`init` is optional, defaults to leave elements undefined.

Access or write as a vector:

    A[5.2] += 2

Numbers above and below range map to two special bins.

If indexed with integers, these are interpreted as bin numbers.
`A[0]` and `A[-1]` are the outlier bins (below and above, respectively).
"""
struct BinnedVector{T} <: AbstractArray{T,1}
    min::Float64
    max::Float64
    nbins::Int
    Δ::Float64
    data::Vector{T}
end

BinnedVector{T}(nbins::Integer;min::AbstractFloat,max::AbstractFloat,init=nothing) where {T} =
    BinnedVector{T}(min,max,nbins,(max-min)/nbins,
                    isnothing(init) ? Vector{T}(undef,nbins+2) : init(T,nbins+2) )

function BinnedVector{T}(;Δ::AbstractFloat,min::AbstractFloat,max::AbstractFloat,
                         round_Δ::Union{RoundingMode,Nothing}=nothing,
                         round_max::Union{RoundingMode,Nothing}=nothing,
                         round_min::Union{RoundingMode,Nothing}=nothing,
                         init=nothing) where {T}

    if !isnothing(round_Δ)
        nbins = round_Δ==RoundDown ? Int(ceil((max-min)/Δ)) : Int(floor((max-min)/Δ))
    elseif !isnothing(round_max)
        nbins = round_max==RoundUp ? Int(ceil((max-min)/Δ)) : Int(floor((max-min)/Δ))
        max = min + Δ*nbins
    elseif !isnothing(round_min)
        nbins = round_min==RoundDown ? Int(ceil((max-min)/Δ)) : Int(floor((max-min)/Δ))
        min = max - Δ*nbins
    else
        throw("BinnedVector: one of round_Δ, round_min, round_max must be different from nothing")
    end
    return BinnedVector{T}(min,max,nbins,(max-min)/nbins,
                    isnothing(init) ? Vector{T}(undef,nbins+2) : init(T,nbins+2) )
end    

Base.size(A::BinnedVector{T}) where {T} = tuple(A.nbins)
Base.size(A::BinnedVector{T},dim) where {T} = dim==1 ? A.nbins : 1

"""
    interval(A::BinnedVector{T})

Return tuple `(min,max)` giving the extrema of the real interval
mapped to the array bins.
"""
interval(A::BinnedVector{T}) where {T} = tuple(A.min,A.max)

"""
    nbins(A::BinnedVector{T})

Return number of bins of `A`
"""
nbins(A::BinnedVector{T}) where {T} = A.nbins

"""
    delta(A::BinnedVector{T})

Return the with of the bins of `A`
"""
delta(A::BinnedVector{T}) where {T} = A.Δ

"""
    bin(A::BinnedVector{T},x::Float64) where {T}

Map real value `x` to bin number.  Return 0 if below range,
or -1 if above range.
"""
function bin(A::BinnedVector{T},x::Float64)::Int where {T}
    b::Int = floor(Int,(x-A.min)/A.Δ)+1
    if b<1 return 0
    elseif b>A.nbins return -1
    else return b
    end
end

"""
    binc(A::BinnedVector{T},i::Int) where {T} 

Return center of bin `i`, which must be in range `1:size(A,1)`. Does
not perform range check.
"""
binc(A::BinnedVector{T},bin::Int) where {T} = A.min + (bin-0.5)*A.Δ

"""
    Base.range(A::BinnedVector{T}) where{T}

Return range spanning bin centers of `A`.  With `r=collect(range(A))`
one obtains a `Vector` such that `r[i]=binc(A,i)`
"""
Base.range(A::BinnedVector{T}) where{T} = range(start=A.min+0.5*A.Δ,step=A.Δ,length=A.nbins)

Base.getindex(A::BinnedVector{T},x::Float64) where {T} =
    A.data[bin(A,x)+2]

Base.getindex(A::BinnedVector{T},i::Int) where {T} =
    A.data[i+2]

Base.setindex!(A::BinnedVector{T},value::T,x::Float64) where {T} =
    A.data[bin(A,x)+2]=value

Base.setindex!(A::BinnedVector{T},value::T,i::Int) where {T} =
    A.data[i+2]=value
