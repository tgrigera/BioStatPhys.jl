# binnings.jl -- Some special cases of BinnedVector
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
    distance_binning(pos,Δr;rmin=0.,rmax=nothing)

Create a BinnedVector of tuples `(i,j)` with `i!=j`, where each bin
contains all pairs of positions in `pos` whose distance is within the
bin extrema.

`pos` is expected to be a `Matrix` where each row `pos[i,:]` is a 2-d
or 3-d position vector.
"""
function distance_binning(pos,Δr;rmin=0.,rmax=nothing)
    if isnothing(rmax)
        xmin,xmax=extrema(pos[:,1])
        ymin,ymax=extrema(pos[:,2])
        zmin,zmax=0,0
        if size(pos,2)==3
            zmin,zmax=extrema(pos[:,3])
        end
        vmin=[xmin,ymin,zmin]
        vmax=[xmax,ymax,zmax]
        rmax=LinearAlgebra.norm(vmax-vmin)
    end
    TII=Tuple{Int,Int}
    binning=BinnedVector{Vector{TII}}(Δ=Δr,min=rmin,max=rmax,round_max=RoundUp,
                               init=(x,n)->[TII[] for _=1:n])
    for i=1:size(pos,1)-1,j=i+1:size(pos,1)
        push!(binning[LinearAlgebra.norm( pos[i,:]-pos[j,:] )],(i,j))
    end
    return binning
end

"Alias for a `BinnedVector` of `Vector{Tuple{Int,Int}}`"
DistanceBinning=BinnedVector{Vector{Tuple{Int,Int}}}
