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


#
# This must be extended to provide a push-interface 
#

include("../tool/region.jl")

mutable struct Density_correlation{R<:Region}
    region::R
    npart::Int
    nconf::Int
    npr::BinnedVector{Int}     # number of pairs at distance r
    Cnpr::BinnedVector{Int}    # cumulative number of pairs
    centers::BinnedVector{Int}     # valid centers
end

# function density_correlation_rdf(region::PeriodicRectangle,pos,Δr;rmax=nothing)
#     rdf = Density_correlation_rdf(
#         region, 0,
#         BinnedVector{Int}(Δ=Δr,min=0.,max=rmax,round_max=RoundUp,init=zeros),
#         BinnedVector{Int}(Δ=Δr,min=0.,max=rmax,round_max=RoundUp,init=zeros),
#         BinnedVector{Int}(Δ=Δr,min=0.,max=rmax,round_max=RoundUp,init=zeros),
#     )
#     return density_correlation_rdf(rdf,pos)
# end

# function density_correlation_rdf(rdf::Density_correlation_rdf{R},pos)  where R <: PeriodicRegion
#     rdf.nconf += 1
#     for i ∈ 1:size(pos,1)-1, j ∈ i+1:size(pos,1)
#         r = sqrt(distancesq(rdf.region,pos[i,:],pos[j,:]))
#         b = bin(rdf.npr,r)
#         rdf.npr[b] += 2
#         if b<0 b=nbins(rdf.Cnpr) end
#         rdf.Cnpr[1:b] .+= 2
#     end
#     return rdf
# end

# function ufa
#     gr = zeros(Float64,size(binning,1))
#     gr .= length.(binning)
#     gr ./= N
#     Δ = delta(binning)
#     for (i,r) ∈ enumerate(range(binning))
#         vol = π * ( (r+Δ)^2 - (r-Δ)^2 )
#         gr[i] /= vol * Nc
#     end
#     return gr
# end

function density_correlation(region::Rectangle,pos,Δr;rmax=nothing)
    if isnothing(rmax)
        rmax=LinearAlgebra.norm([region.Lx,region.Ly])
    end
    dcorr = Density_correlation(
        region, size(pos,1),0,
        BinnedVector{Int}(Δ=Δr,min=0.,max=rmax,round_max=RoundUp,init=zeros),
        BinnedVector{Int}(Δ=Δr,min=0.,max=rmax,round_max=RoundUp,init=zeros),
        BinnedVector{Int}(Δ=Δr,min=0.,max=rmax,round_max=RoundUp,init=zeros),
    )
    return density_correlation(dcorr,pos)
end

function density_correlation(dcorr::Density_correlation{Rectangle},pos)
    @assert dcorr.npart==size(pos,1)
    dcorr.nconf += 1
    for i ∈ 1:size(pos,1)
        dbd = dborder(dcorr.region,pos[i,1], pos[i,2])
        lsb=bin(dcorr.npr,dbd)-1
        dcorr.centers[1:lsb] .+= 1
        for j ∈ 1:size(pos,1)
            if i==j continue end
            r = LinearAlgebra.norm( pos[i,:]-pos[j,:] )
            b = bin(dcorr.npr,r)
            if 1 <= b <= lsb
                dcorr.Cnpr[1:b] .+= 1
                dcorr.npr[b] += 1
            else
                dcorr.Cnpr .+= 1
            end
        end
    end
    return dcorr
end

function rdf(dcorr::Density_correlation{R}) where R <: Region
    rmax = interval(dcorr.npr)[2]
    Δ = delta(dcorr.npr)
    gr = BinnedVector{Float64}(Δ=Δ,min=0.,max=rmax,round_max=RoundUp,init=zeros)
    Cr = BinnedVector{Float64}(Δ=Δ,min=0.,max=rmax,round_max=RoundUp,init=zeros)
    ρ = dcorr.npart / volume(dcorr.region)
    for (i,r) ∈ enumerate(range(dcorr.npr))
        vol = π * ( (r+Δ/2)^2 - (r-Δ/2)^2 )
        gr[i] = dcorr.npr[i] / (dcorr.npart * ρ * vol * dcorr.centers[i] * dcorr.nconf)
        Cr[i] = Cr[i-1] + vol*gr[i]
        # Cr[i] = rdf.Cnpr[i] / (rdf.centers[i]* rdf.nconf)
    end
    return gr,Cr
end

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
