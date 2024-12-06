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

import SpecialFunctions

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

mutable struct KSpace_correlation_vectorqty
    k::AbstractRange{Float64}
    const Cv::Vector{Float64}
    const Cs::Vector{Float64}
    const Cmod::Vector{Float64}
end

mutable struct Space_correlation_vectorqty{R<:Region}
    connected::Bool
    const region::R
    nconf::Int
    npart::Int  # Total number of particles processed (not just particles per config)
    const npairs::ZBinnedVector{Int}         # Number of pairs at given distance
    const npairs_nz::ZBinnedVector{Int}      # Number of pairs with non-zero modulus
    const Cv::ZBinnedVector{Float64}
    const Cs::ZBinnedVector{Float64}
    const Cmod::ZBinnedVector{Float64}
    bin::ZDistanceBinning
    kspacedata::Union{KSpace_correlation_vectorqty,Nothing}
end

"""
    space_correlation(region::Region,Δr;rmax=nothing,connected=false)
    space_correlation(region::Region,pos::ConfigurationT,Δr;rmax=nothing,connected=false)

Build and return a `Space_correlation_vectorqty` object to compute
space correlations of some vectorial observable.  If `pos` is given, it
will be assumed that positions will be the same for all configurations,
and a `DistanceBinning` object will be built to save time later.  If posiitons
change with every configuration, then omit this argument.

To add data to an object `corr`, call

    space_correlation!(corr,pos,vec)

or

    space_correlation!(corr,vec)

if particles are static.  Correlation functions are obtained calling
`correlations(corr)`.

"""
function space_correlation(region::Region,Δr;rmax=nothing,connected=false,kspace=false)
    if isnothing(rmax)
        rmax,_ = linear_size(region)
        if isa(region,PeriodicRegion) rmax /= 2. end
    end
    scorr = Space_correlation_vectorqty{Region}(
        connected, region, 0, 0,
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Float64}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Float64}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Float64}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZDistanceBinning(Δ=1.,max=1.,round_max=RoundUp),
        nothing
    )
    if kspace 
        Δk = 2π/rmax
        kmax = 2π/delta(scorr.npairs)
        Nk = ceil(Int,kmax/Δk)
        kcorr = KSpace_correlation_vectorqty(
            range(start=0,step=Δk,length=Nk),
            zeros(Float64,Nk),
            zeros(Float64,Nk),
            zeros(Float64,Nk)
        )
        scorr.kspacedata = kcorr

    end
    return scorr
end

function space_correlation!(corr::Space_correlation_vectorqty{R},pos::ConfigurationT,
    vec::ConfigurationT) where R<:Region

    sv = vec ./ LinearAlgebra.norm.(vec)
    sv = replace(sv,[NaN,NaN]=>missing)
    if corr.connected
        vmean = Statistics.mean(vec)
        skp=collect(skipmissing(sv))
        if length(skp)>0 smean = Statistics.mean(skp) end
        modmean = Statistics.mean(x->LinearAlgebra.norm(x),vec)
    else
        vmean, smean, modmean = [0., 0.], [0.,0.], 0.
    end
    
    corr.nconf += 1
    N = size(pos,1)
    Nnz = count(!ismissing,sv)
    if !isnothing(corr.kspacedata)
        Ckv = zeros(Float64,length(corr.kspacedata.Cv))
        Cks = zeros(Float64,length(corr.kspacedata.Cs))
        Ckmod = zeros(Float64,length(corr.kspacedata.Cmod))
    end

    for i ∈ eachindex(pos)
        corr.npart += 1
        ri = pos[i]
        vi = vec[i] - vmean
        if !ismissing(sv[i]) si = sv[i] - smean end
        modi = LinearAlgebra.norm(vec[i]) - modmean

        corr.npairs[0.] += 1
        corr.Cv[0.] += vi ⋅ vi
        if !ismissing(sv[i]) 
            corr.npairs_nz[0.] += 1
            corr.Cs[0.] += si ⋅ si
        end
        corr.Cmod[0.] += modi * modi

        if !isnothing(corr.kspacedata)
            Ckv .+= (vi ⋅ vi)
            Ckmod .+= (modi*modi)
            if !ismissing(sv[i])
                Cks .+= (si ⋅ si)
            end
        end

        for j ∈ i+1:size(pos,1)
            dr = distance(corr.region,pos[j],ri)
            corr.npairs[dr] += 2
            vj = vec[j] - vmean
            corr.Cv[dr] += 2 * vi ⋅ vj
            if !ismissing(sv[i]) && !ismissing(sv[j])
                corr.npairs_nz[dr] += 2
                sj = sv[j] - smean
                corr.Cs[dr] += 2 * si ⋅ sj
            end
            modj = LinearAlgebra.norm(vec[j]) - modmean
            corr.Cmod[dr] += 2 * modi * modj

            if !isnothing(corr.kspacedata)
                kr = dr .* corr.kspacedata.k
                xx = x->isoexp(x,corr.region)
                Ckv .+= 2 .* (vi ⋅ vj .* xx.(kr) )
                Ckmod .+= 2 .* (modi*modj .* xx.(kr) )
                if !ismissing(sv[i]) && !ismissing(sv[j])
                    Cks .+= 2 .* (si ⋅ sj .* xx.(kr) )
                end
            end
        end
    end

    if !isnothing(corr.kspacedata)
        corr.kspacedata.Cv .+= Ckv ./ N
        if Nnz>0 corr.kspacedata.Cs .+= Cks ./ Nnz end
        corr.kspacedata.Cmod .+= Ckmod ./ N
    end

end

isoexp(_,::Region{D}) where D = throw(ErrorException("kspace correlation not implemented for dimension $D"))

isoexp(kr,::Region{1}) = exp(im*kr)

isoexp(kr,::Region{2}) = SpecialFunctions.besselj0(kr)

isoexp(kr,::Region{3}) = sinc(kr/π)

function space_correlation(region::Region,pos::ConfigurationT,Δr;rmax=nothing,connected=false)
    if isnothing(rmax)
        rmax,_ = linear_size(region)
        if isa(region,PeriodicRegion) rmax /= 2. end
    end
    scorr = space_correlation(region,Δr,rmax=rmax,connected=connected)
    scorr.bin = distance_binning(region,pos,Δr,rmax=rmax)
    return scorr
end

function space_correlation!(corr::Space_correlation_vectorqty{R},
    vec::ConfigurationT) where R<:Region

    sv = vec ./ LinearAlgebra.norm.(vec)
    sv = replace(sv,[NaN,NaN]=>missing)
    if corr.connected
        vmean = Statistics.mean(vec)
        skp=collect(skipmissing(sv))
        if length(skp)>0 smean = Statistics.mean(skp) end
        modmean = Statistics.mean(x->LinearAlgebra.norm(x),vec)
    else
        vmean, smean, modmean = [0., 0.], [0.,0.], 0.
    end
    
    corr.nconf += 1
    for ib in eachindex(corr.bin)
        for (i,j) in corr.bin[ib]
            corr.npairs[ib] += 1

            vi = vec[i] - vmean
            vj = vec[j] - vmean
            corr.Cv[ib] += vi ⋅ vj

            if !ismissing(sv[i]) si = sv[i] - smean
                if !ismissing(sv[j])
                    corr.npairs_nz[ib] += 1
                    sj = sv[j] - smean
                    corr.Cs[ib] += si ⋅ sj
                end
            end
            
            modi = LinearAlgebra.norm(vec[i]) - modmean
            modj = LinearAlgebra.norm(vec[j]) - modmean
            corr.Cmod[ib] += modi * modj
        end
    end

end

function correlations(C::Space_correlation_vectorqty{R}; normalize=false) where R<:Region
    r = collect(range(C.Cv))
    Cv = C.Cv ./ C.npairs
    Cs = C.Cs ./ C.npairs_nz
    Cmod = C.Cmod ./ C.npairs

    if normalize
        Cv[2:end] ./= Cv[1]
        Cv[1] = 1.
        Cs[2:end] ./= Cs[1]
        Cs[1] = 1.
        Cmod[2:end] ./= Cmod[1]
        Cmod[1] = 1.
    end

    if isnothing(C.kspacedata) return r,Cv,Cs,Cmod end
    k = collect(C.kspacedata.k)
    Ckv = C.kspacedata.Cv ./ C.nconf
    Cks = C.kspacedata.Cs ./ C.nconf
    Ckmod = C.kspacedata.Cmod ./ C.nconf
    return (r,Cv,Cs,Cmod),(k,Ckv,Cks,Ckmod)
end

###############################################################################
#
# Density correlations
#

"""
    Density_correlation{R<:Region}

Datatype holding information about about pair distribution from which
the density correlations are computed.  Not to be manipulated
directly, but objects of this type are returned by the
`density_correlation` methods, and are taken as argument to compute
the actual correlation functions, e.g. `rdf`.
"""
mutable struct Density_correlation{R<:Region}
    region::R
    ρ::MeanVar
    nconf::Int
    npr::ZBinnedVector{Int}         # number of pairs at distance r
    centers::ZBinnedVector{Int}     # valid centers
end

"""
    density_correlation(region,pos;Δr,rmax)

Return a `Density_correlation` object for the given `region` and
configuration `pos`, with resolution `Δr` and maximum range `rmax`.
The returned object can be passed to the functions that compute the
final correlations, such as `rdf`, or more configurations can be added
to it by calling `density_correlation(::DensityCorrelation,pos)`.

`pos` must be a vector of vectors, i.e. each element of `pos` must be
a vector of the appropriate dimensionality.  Performance advantage may
be obtained using `StaticArrays` to represent individual positions.
More precisely, the type of `pos` is

    AbstractVector{T} where T<:AbstractVector{W} where W<:Number

"""
function density_correlation(region::PeriodicRegion;Δr,rmax=nothing)
    if isnothing(rmax)
        rmax = 0.5 * volume(region)^(1/dimension(region))
    end
    dcorr = Density_correlation(
        region,MeanVar(),0,
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros)
    )
    return dcorr
end

density_correlation(region::PeriodicRegion,pos::ConfigurationT;Δr,rmax=nothing) =
    density_correlation!(density_correlation(region,Δr=Δr,rmax=rmax),pos)

"""
    density_correlation!(dcorr::Density_correlation{R},pos) where R <: Region

Use configuration `pos` to add statistics to the `DensityCorrelation`
object `dcorr`.  The configuration need not have the same number of particles
as those used previously, but check with functions `rdf` and the like how these
fluctuations are treated.
"""
function density_correlation!(dcorr::Density_correlation{R},pos::ConfigurationT) where R <: PeriodicRegion
    dcorr.nconf += 1
    npart = size(pos,1)
    dcorr.npr[0.] += size(pos,1)
    for i ∈ 1:npart-1, j ∈ i+1:npart
        r = sqrt(distancesq(dcorr.region,pos[i],pos[j]))
        dcorr.npr[r] += 2
    end
    push!(dcorr.ρ,npart/volume(dcorr.region))
    return dcorr
end

"""
    rdf(dcorr::Density_correlation{R};two_particle_density=false)

If `two_particle_density`is `false` (default), compute the radial
distribution function ``g(r)`` and the correlation integral ``C(r)``.
Both are returned as a tuple of `ZBinnedVector`.

If `two_particle_density` is `true`, return the density-density
correlation function ``G(r)`` instead of ``g(r)``.  For the isotropic
case assumed here, ``G(r) = \\rho^2 g(r) + \\rho \\delta(r)/4\\pi r^2``
(in 3-d).  Since ``G(r)`` is singular at ``r=0``, here for ``r=0`` we
return ``\\rho``, i.e. the integral of ``G(r)`` in a very small volume
around the origin.

If the `dcorr` object was fed with configurations with fluctuating
number of particles, then the radial distribution function will be
computed using the definition for the grand canonical ensamble (see
J.-P. Hansen and I. R. McDonald, _Theory of Simple Liquids_, Academic
Press (2005)).  This may or may not be what you want.  See the online
documentation of the package for more details.
"""
function rdf(dcorr::Density_correlation{R};two_particle_density=false) where R <: PeriodicRegion
    rmax = interval(dcorr.npr)[2]
    Δ = delta(dcorr.npr)
    gr = ZBinnedVector{Float64}(Δ=Δ,max=rmax,round_max=RoundUp,init=zeros)
    Cr = ZBinnedVector{Float64}(Δ=Δ,max=rmax,round_max=RoundUp,init=zeros)
    ρ = mean(dcorr.ρ)
    fac = two_particle_density ? volume(dcorr.region) * dcorr.nconf :
          ρ^2 * volume(dcorr.region) * dcorr.nconf
    Cr[1] = dcorr.npr[1] / dcorr.nconf # Because C(r) does not exclude self-correlation
    for (i,r) ∈ enumerate(range(dcorr.npr))
        vols = shell_volume(r,Δ,Val(dimension(dcorr.region)))
        if i==1 continue end
        gr[i] = dcorr.npr[i] / (fac * vols)
        Cr[i] = Cr[i-1] + dcorr.npr[i] / dcorr.nconf
    end
    gr[1] = two_particle_density ? ρ : 0.
    Cr[1] = 0.
    Cr ./= (ρ*volume(dcorr.region))^2
    return gr,Cr
end

function density_correlation(region::NonPeriodicRegion,pos::ConfigurationT;Δr,rmax=nothing)
    if isnothing(rmax)
        rmax = maximum(region.L-region.x0)
    end
    dcorr = Density_correlation(
        region,MeanVar(),0,
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros)
    )
    return density_correlation!(dcorr,pos)
end

function density_correlation!(dcorr::Density_correlation{<:NonPeriodicRegion},pos::ConfigurationT)
    npart = size(pos,1)
    push!(dcorr.ρ,npart/volume(dcorr.region))
    dcorr.nconf += 1
    for i ∈ 1:size(pos,1)
        dbd = dborder(dcorr.region,pos[i])
        lsb=bin(dcorr.npr,dbd)-1
        dcorr.centers[1:lsb] .+= 1
        for j ∈ 1:size(pos,1)
            r = LinearAlgebra.norm( pos[i]-pos[j] )
            b = bin(dcorr.npr,r)
            if 1 <= b <= lsb
                dcorr.npr[b] += 1
            end
        end
    end
    return dcorr
end

shell_volume(_,Δ,::Val{1}) = Δ
shell_volume(r,Δ,::Val{2}) = π * ( (r+Δ/2)^2 - (r-Δ/2)^2 )
shell_volume(r,Δ,::Val{3}) = (4π/3) * ( (r+Δ/2)^3 - (r-Δ/2)^3 )

function rdf(dcorr::Density_correlation{R};two_particle_density=false) where R <: NonPeriodicRegion
    rmax = interval(dcorr.npr)[2]
    Δ = delta(dcorr.npr)
    gr = ZBinnedVector{Float64}(Δ=Δ,max=rmax,round_max=RoundUp,init=zeros)
    Cr = ZBinnedVector{Float64}(Δ=Δ,max=rmax,round_max=RoundUp,init=zeros)
    ρ = mean(dcorr.ρ)
    Cr[1] = dcorr.npr[1] /  dcorr.centers[1]  # Because C(r) does not exclude self-correlation
    for (i,r) ∈ enumerate(range(dcorr.npr))
        vol = shell_volume(r,Δ,Val(dimension(dcorr.region)))
        if i==1 continue end
        gr[i] = dcorr.npr[i] / (ρ * vol * dcorr.centers[i])
        Cr[i] = Cr[i-1] + dcorr.npr[i] / dcorr.centers[i]
    end
    if two_particle_density
        gr .*= ρ^2
        gr[1] = ρ
    end
    Cr[1] = 0.
    Cr ./= (ρ*volume(dcorr.region))^2
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
