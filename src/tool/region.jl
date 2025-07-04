# region.jl -- Geometrical regions with different boundary conditions
#
# Copyright (C) 2023 by Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
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

using StaticArrays

"""
    abstract type Region

Base type for D-dimensional regions.
"""
abstract type Region{D} end

"""
    dimension(r<:Region)

Return the dimension of region `r`, in the sense of dimension of a
space or a manifold.
"""
dimension(::Region{D}) where D = D

###############################################################################
#
# Non-periodic regions
#

abstract type NonPeriodicRegion{D} <: Region{D} end

distance(::NonPeriodicRegion,x::AbstractVector{<:Number},y::AbstractVector{<:Number}) =
    LinearAlgebra.norm(x.-y)

"""
    HyperCube{D} <: NonPeriodicRegion

Describes a nonperiodic D-dimensional region delimited by
perpendicular planes.  Size _and_ origin are arbitrary.  The size
along each dimension can be different.

To create a hypercube with side 10 and origin at `x0` do e.g.

    HC = HyperCube{4}(10.,10.,10.,10.,x0=(0.,0.,2.,2.))

"""
struct HyperCube{D} <: NonPeriodicRegion{D}
    x0::SVector{D,Float64}
    L::SVector{D,Float64}
end

HyperCube{D}(L...;x0) where D = HyperCube{D}(SVector(x0),SVector(L...))

"Return the volume of the given region"
function volume(r::HyperCube{D}) where D
    vol = 1.
    for i ∈ 1:D vol *= r.L[i]  end
    return vol
end

"""
    Rectangle <: NonPeriodicRegion

Describes a rectangular non-periodic 2-d region with an arbitrary origin and
size.  Create with

    rec = Rectangle(Lx,Ly,x0=(x0,y0))
"""
Rectangle = HyperCube{2}

function HyperCube{D}(pos::ConfigurationT) where D
    x0 = MVector{D,Float64}(undef)
    L = MVector{D,Float64}(undef)
    for i ∈ 1:D
        x0[i],L[i] = extrema(map(x->x[i],pos))
        L[i] -= x0[i]
    end
    return HyperCube{D}(L...,x0=x0)
end


"""
    linear_size(r<:Region)

Return minimum and maximum linear size of region
"""
linear_size(r::Region) = extrema(r.L)


"""
    dborder(r<:NonPeriodicRegion,p::AbstractVector)

Return the distance from the point `p` (assumed included in region
`r`) to the nearest border.
"""
function dborder(r::HyperCube{D},p::AbstractVector{<:Number}) where D
    dis = MVector{2*D,Float64}(undef)
    dis[1:D] .= p .- r.x0
    dis[D+1:end] .= r.x0 .+ r.L .- p
    return minimum(dis)
end

dborder(r::Rectangle,p::AbstractVector{<:Number}) =
    minimum( SVector( p[1]-r.x0[1], r.x0[1]+r.L[1]-p[1],
        p[2]-r.x0[2], r.x0[2]+r.L[2]-p[2] ) )

"""
    struct Cube <: NonPeriodicRegion

Describes a non-periodic 3-d region with cubic symmetry (a
rectangular prism).
"""
Cube = HyperCube{3}

dborder(r::Cube,p::AbstractVector{<:Number}) =
    minimum( SVector( p[1]-r.x0[1], r.x0[1]+r.L[1]-p[1],
        p[2]-r.x0[2], r.x0[2]+r.L[2]-p[2],
        p[3]-r.x0[3], r.x0[3]+r.L[3]-p[3] ) )

###############################################################################
#
# Periodic regions
#

abstract type PeriodicRegion{D} <: Region{D} end

"""
    PeriodicHyperCube{D} <: PeriodicRegion

Describes a periodic D-dimensional region.  Actually the size along each
dimension can be different.
"""
struct PeriodicHyperCube{D} <: PeriodicRegion{D}
    L::SVector{D,Float64}
end

PeriodicHyperCube{D}(L...) where D = PeriodicHyperCube{D}(SVector(L...))

dimension(::PeriodicHyperCube{D}) where D = D

function volume(r::PeriodicHyperCube{D}) where D
    vol = 1.
    for i ∈ 1:D vol *= r.L[i]  end
    return vol
end

"""
    PeriodicRectangle <: PeriodicRegion

Describes periodic rectangular 2-d region with arbitrary size.  Alias
for `PeriodicHyperCube{2}`.
"""
PeriodicRectangle = PeriodicHyperCube{2}

"""
    PeriodicCube(Lx,Ly,Lz)

Return a periodic cube (actually rectangular cuboid).  Alias for `PeriodicHyperCube{3}`
"""
PeriodicCube = PeriodicHyperCube{3}

function ddiff(a,b,box_length)
  temp = a-b
  return temp - box_length*round(temp/box_length)
end

"""
    distancesq(r<:PeriodicRegion,x::AbstractVector{<:Number},y::AbstractVector{<:Number})

Return the periodic squared distance between points `x` and `y` in the periodic region `r`
"""
function distancesq(r::PeriodicHyperCube{D},x::AbstractVector{<:Number},
                    y::AbstractVector{<:Number}) where D
    dsq = 0.
    for i ∈ 1:D
        dx = ddiff(x[i],y[i],r.L[i])
        dsq += dx*dx
    end
    return dsq
end

distance(r::PeriodicHyperCube{D},x::AbstractVector{<:Number},y::AbstractVector{<:Number}) where D =
    sqrt(distancesq(r,x,y))


fold!(x::AbstractArray,box::AbstractArray) =
    for i ∈ eachindex(x) x[i] -= box[i]*floor(x[i]/box[i]) end

"""
    fold!(r::PeriodicRegion,conf::ConfigurationT)

Translate periodically the positions in configuration `conf` so that
all of them fall within the bounds of region `r` (i.e. in the
'original' periodic volume as defined in `r`).
"""    
function fold!(r::PeriodicRegion,conf::ConfigurationT)
    for n ∈ eachindex(conf)
        fold!(conf[n],r.L)
    end
end
