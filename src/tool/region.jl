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

"""
    abstract type Region

Base type for regions.
"""
abstract type Region end

###############################################################################
#
# Non-periodic regions
#

abstract type NonPeriodicRegion <: Region end

distance(::NonPeriodicRegion,x::AbstractVector{<:Number},y::AbstractVector{<:Number}) =
    LinearAlgebra.norm(x.-y)

"""
    Rectangle <: NonPeriodicRegion

Describes a rectangular non-periodic 2-d region with an arbitrary origin and
size.
"""
struct Rectangle <: NonPeriodicRegion
    x0::Float64   # Origin
    y0::Float64
    Lx::Float64   # Length
    Ly::Float64
end

function Rectangle(pos::ConfigurationT)
    xmin,xmax = extrema(map(x->x[1],pos))
    ymin,ymax = extrema(map(x->x[2],pos))
    return Rectangle(xmin,ymin,xmax-xmin,ymax-ymin)
end

"""
    dimension(r<:Region)

Return the dimension of region `r`, in the sense of dimension of a
space or a manifold.
"""
dimension(::Rectangle) = 2

"""
    dborder(r<:NonPeriodicRegion,p::AbstractVector)

Return the distance from the point `p` (assumed included in region
`r`) to the nearest border.
"""
dborder(r::Rectangle,p::AbstractVector{<:Number}) =
    minimum( [ p[1]-r.x0, r.x0+r.Lx-p[1], p[2]-r.y0, r.y0+r.Ly-p[2]] )

"Return volume of region"
volume(r::Rectangle) = r.Lx * r.Ly

"""
    struct CubicBox <: NonPeriodicRegion

Describes a non-periodic 3-d region with cubic symmetry (a
rectangular prism).
"""
struct CubicBox <: NonPeriodicRegion
    x0::Float64   # Origin
    y0::Float64
    z0::Float64
    Lx::Float64   # Length
    Ly::Float64
    Lz::Float64
end

function CubicBox(pos::ConfigurationT)
    xmin,xmax = extrema(map(x->x[1],pos))
    ymin,ymax = extrema(map(x->x[2],pos))
    zmin,zmax = extrema(map(x->x[3],pos))
    return CubicBox(xmin,ymin,zmin,xmax-xmin,ymax-ymin,zmax-zmin)
end

dimension(::CubicBox) = 3

dborder(r::CubicBox,p::AbstractVector{<:Number}) =
    minimum( [p[1]-r.x0, r.x0+r.Lx-p[1], p[2]-r.y0, r.y0+r.Ly-p[2], p[3]-r.z0, r.z0+r.Lz-p[3]] )

volume(r::CubicBox) = r.Lx * r.Ly * r.Lz


###############################################################################
#
# Periodic regions
#

abstract type PeriodicRegion <: Region end

"""
    PeriodicRectangle <: PeriodicRegion

Describes periodic rectangular 2-d region with arbitrary size.
"""
struct PeriodicRectangle <: PeriodicRegion
    Lx::Float64
    Ly::Float64
end

dimension(::PeriodicRectangle) = 2

volume(r::PeriodicRectangle) = r.Lx * r.Ly

function ddiff(a,b,box_length)
  temp = a-b
  return temp - box_length*round(temp/box_length)
end

"""
    distancesq(r<:PeriodicRegion,x::AbstractVector{<:Number},y::AbstractVector{<:Number})

Return the periodic squared distance between points `x` and `y` in the periodic region `r`
"""
function distancesq(r::PeriodicRectangle,x::AbstractVector{<:Number},y::AbstractVector{<:Number})
    dx = ddiff(x[1],y[1],r.Lx)
    dy = ddiff(x[2],y[2],r.Ly)
    return dx*dx + dy*dy
end
