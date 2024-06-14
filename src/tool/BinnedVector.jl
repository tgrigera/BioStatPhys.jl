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
    BBinnedVector{T,SpecialZero} <: AbstractArray{T,1}

This is a "base" type that fully implements `BinnedVector` and
`ZBinnedVector`, which are simply aliases for `BBinnedVector{T,false}`
and `BinnedVector{T,true}'.

See `BinnedVector` and `ZBinnedVector` for documentation.
"""
mutable struct BBinnedVector{T,SpecialZero} <: AbstractArray{T,1}
    min::Float64
    max::Float64
    nbins::Int
    Δ::Float64
    data::Vector{T}
end

BBinnedVector{T,SZ}(nbins::Integer;min::AbstractFloat=0.,max::AbstractFloat,
                init=nothing) where {T,SZ} =
                    BBinnedVector{T,SZ}(SZ ? 0. : min,max,nbins,(max-min)/nbins,
                                    isnothing(init) ? Vector{T}(undef,nbins+2) : init(T,nbins+2) )

"""
    BinnedVector{T}

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
const BinnedVector{T} = BBinnedVector{T,false}

function BinnedVector{T}(;Δ::Real,min::Real,max::Real,
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

"""
    ZBinnedVector{T}

This type works like `BinnedVector{T}` except that the range always
starts with `0.`, and that the `0.` is treated specially.  The first
bin excludes the zero (i.e., unlike the rest of the bins is an open
interval on both sides), and the values indexed with 0. are treated
separately.  When using integer indices, index 1 corresponds to 0.,
index 2 to the first bin, and so on.

This is useful when dealing for example with space correlation
functions, where a zero distance implies a correlation of a particle
with itself, which is often convenient to treat separately from the
inter-particle correlations.
"""
const ZBinnedVector{T} = BBinnedVector{T,true}

function ZBinnedVector{T}(;Δ::Real,max::Real,
                         round_Δ::Union{RoundingMode,Nothing}=nothing,
                         round_max::Union{RoundingMode,Nothing}=nothing,
                         init=nothing) where {T}

    if !isnothing(round_Δ)
        nbins = round_Δ==RoundDown ? Int(ceil((max)/Δ)) : Int(floor((max)/Δ))
    elseif !isnothing(round_max)
        nbins = round_max==RoundUp ? Int(ceil((max)/Δ)) : Int(floor((max)/Δ))
        max = Δ*nbins
    else
        throw("ZBinnedVector: one of round_Δ, round_max must be different from nothing")
    end
    nbins += 1
    return ZBinnedVector{T}(0.,max,nbins,max/(nbins-1),
                            isnothing(init) ? Vector{T}(undef,nbins+2) : init(T,nbins+2) )
end    

Base.size(A::BBinnedVector{T,SZ}) where {T,SZ} = tuple(A.nbins)
Base.size(A::BBinnedVector{T,SZ},dim) where {T,SZ} = dim==1 ? A.nbins : 1

"""
    interval(A::BinnedVector{T})
    interval(A::ZBinnedVector{T})

Return tuple `(min,max)` giving the extrema of the real interval
mapped to the array bins.
"""
interval(A::BBinnedVector{T,SZ}) where {T,SZ} = tuple(A.min,A.max)

"""
    nbins(A::BinnedVector{T})
    nbins(A::ZBinnedVector{T})

Return number of bins of `A`
"""
nbins(A::BBinnedVector{T,SZ}) where {T,SZ} = A.nbins

"""
    delta(A::BinnedVector{T})
    delta(A::ZBinnedVector{T})

Return the with of the bins of `A`
"""
delta(A::BBinnedVector{T,SZ}) where {T,SZ} = A.Δ

"""
    bin(A::BinnedVector{T},x::Float64) where {T}
    bin(A::ZBinnedVector{T},x::Float64) where {T}

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

function bin(A::ZBinnedVector{T},x::Float64)::Int where {T}
    if x==0. return 1 end
    b::Int = floor(Int,(x-A.min)/A.Δ)+2
    if b<2 return 0
    elseif b>A.nbins return -1
    else return b
    end
end

"""
    binc(A::BinnedVector{T},i::Int) where {T} 
    binc(A::ZBinnedVector{T},i::Int) where {T} 

Return center of bin `i`, which must be in range `1:size(A,1)`. Does
not perform range check.
"""
binc(A::BinnedVector{T},bin::Int) where {T} = A.min + (bin-0.5)*A.Δ

binc(A::ZBinnedVector{T},bin::Int) where {T} = bin==1 ? 0. : (bin-1-0.5)*A.Δ

Base.getindex(A::BBinnedVector{T,SZ},x::Float64) where {T,SZ} =
    A.data[bin(A,x)+2]

Base.getindex(A::BBinnedVector{T,SZ},i::Int) where {T,SZ} =
    A.data[i+2]

Base.setindex!(A::BBinnedVector{T,SZ},value,x::Float64) where {T,SZ} =
    A.data[bin(A,x)+2]=value

Base.setindex!(A::BBinnedVector{T,SZ},value,i::Int) where {T,SZ} =
    A.data[i+2]=value

"""
    Base.range(A::BinnedVector{T}) where{T}
    Base.range(A::ZBinnedVector{T}) where{T}

Return iterator spanning bin centers of `A`.  With `r=collect(range(A))`
one obtains a `Vector` such that `r[i]=binc(A,i)`.

For `BinnedVector` this is a standard range, for `ZBinnedVector` it is
a range-compatible custom iterator since the distance from 0. to the
first bin center is half the distance between bin centers.
"""
Base.range(A::BinnedVector{T}) where{T} = range(start=A.min+0.5*A.Δ,step=A.Δ,length=A.nbins)

struct ZBViter range end

Base.range(A::ZBinnedVector{T}) where{T} =
    ZBViter(range(start=A.min+0.5*A.Δ,step=A.Δ,length=A.nbins-1))

import Base.eltype, Base.show

Base.eltype(::ZBViter) = Float64

Base.show(io::IO, r::ZBViter) = print(io, "Iterator for ZBinnedVector, [0.,",r.range,"]")

import Base.iterate

Base.iterate(::ZBViter) = (0.,nothing)

Base.iterate(it::ZBViter,state) = state === nothing ? iterate(it.range) : iterate(it.range,state)

import Base.length

Base.length(it::ZBViter) = 1+length(it.range)
