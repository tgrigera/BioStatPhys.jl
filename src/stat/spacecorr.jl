# spacecorr.jl -- computation of space correlations
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
    space_correlation(binning::DistanceBinning, X; connected=true, normalized=false, Xmean=nothing)

Compute space correlation of scalar signal `X` in an isotropic,
space-translation-invariant system using the given `DistanceBinning`.
Since the system is assumed to be isotropic and homogeneous,
correlations should depend only on the distances among the pairs, so
return a tuple `(r,C)` where `r` is a vector of distances and `C` the
correlation function at the corresponding distance.

Positions are assumed to be fixed in all configurations, with only the
signal changing value.

 - `X`: a matrix where each column holds the signal at different
   positions in space and at a given time.  All columns are used to
   average the correlation estimate over configurations.

 - `binning`: instead of positions, this function expects a pre-built
   [`DistanceBinning`](@ref) object holding lists of pairs separated
   by distances in the given range, such as the binnings created by
   [`distance_binning`](@ref).

 - `connected`: If true, compute the connected correlation,
   i.e. subtracting the mean of the signal.  If `Xmean` is `nothing`,
   space-averaging is used, i.e. at each configuration (column) the
   mean of `X` is computed and subtracted.  For phase averaging, give
   a precomputed mean in `Xmean`.

 - `normalized`: If `true`, normalize by ``C(r=0)``.

"""
function space_correlation(binning::DistanceBinning,X;
                           connected=true,normalized=false,Xmean=nothing)
    Cr=zeros(Float64,size(binning,1))
    Cr0=0
    nt=size(X,2)
    for time=1:nt
	sg::Vector{Float64}=X[:,time]
        if connected
            if isnothing(Xmean)
	        sg .-= Statistics.mean(sg)
            else
                sg .-= Xmean
            end
        end
	Cr0+= Statistics.mean(sg.*sg)
	for ib in eachindex(binning)
	    for (i,j) in binning[ib]
		Cr[ib]+=sg[i]*sg[j]
	    end
	end
    end
    Cr0/=nt
    if normalized
	Cr./=Cr0*nt*length.(binning)
	Cr0=1
    else
	Cr./=nt*length.(binning)
    end
    return vcat(0,collect(range(binning))),vcat(Cr0,Cr)
end

#    Ccr_space(positions,signal,Δr;rmax=2000,normalize=false)


"""
    correlation_length_r0(r,C)

Compute the correlation length proxy ``r_0`` from ``C_c(r)`` given as
vectors `r` and `C` (as obtained e.g. from `space_correlation`(@ref).

Note that ``r_0`` is **not** a correlation length, just a proxy that
scales with system size ``L`` as ``\\log L`` or ``L`` if the actual
correlation length is much shorter or much larger than ``L``
respectively.
"""
function correlation_length_r0(r,C)
    ir=findfirst(x->x<0,C)
    if !isnan(C[ir-1]) return r[ir-1] - C[ir-1]*(r[ir]-r[ir-1])/(C[ir]-C[ir-1])
    else return r[ir]
    end
end
