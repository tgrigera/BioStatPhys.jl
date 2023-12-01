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

abstract type Region end

struct Rectangle <: Region
    x0::Float64   # Origin
    y0::Float64
    Lx::Float64   # Length
    Ly::Float64
end

function Rectangle(pos::Matrix{<:Number})
    xmin,xmax = extrema(pos[:,1])
    ymin,ymax = extrema(pos[:,2])
    return Rectangle(xmin,ymin,xmax-xmin,ymax-ymin)
end

dborder(r::Rectangle,x::Number,y::Number) = minimum( [x-r.x0, r.x0+r.Lx-x, y-r.y0, r.y0+r.Ly-y] )

volume(r::Rectangle) = r.Lx * r.Ly

abstract type PeriodicRegion <: Region end

struct PeriodicRectangle <: PeriodicRegion
    Lx::Float64
    Ly::Float64
end

volume(r::PeriodicRectangle) = r.Lx * r.Ly

function ddiff(a::Float64,b::Float64,box_length::Float64)
  temp = a-b
  return temp - box_length*round(temp/box_length)
end

function distancesq(r::PeriodicRectangle,x::Number,y::Number)
    dx = ddiff(x[1],y[1],r.Lx)
    dy = ddiff(x[2],y[2],r.Ly)
    return dx*dx + dy*dy
end

