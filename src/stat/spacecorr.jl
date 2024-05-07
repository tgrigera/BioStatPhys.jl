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

mutable struct Space_correlation_vectorqty{R<:Region}
    connected::Bool
    const region::R
    nconf::Int
    npart::Int  # Total number of particles processed (not just particles per config)
    const npairs::ZBinnedVector{Int}         # Number of pairs at given distance
    const npairs_nz::ZBinnedVector{Int}       # Number of pairs with non-zero modulus
    const Cv::ZBinnedVector{Float64}
    const Cs::ZBinnedVector{Float64}
    const Cmod::ZBinnedVector{Float64}
    bin::ZDistanceBinning
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
function space_correlation(region::Region,Δr;rmax=nothing,connected=false)
    if isnothing(rmax)
        rmax,_ = linear_size(region)
        if isa(region,PeriodicRegion) rmax /= 2. end
    end
    scorr = Space_correlation_vectorqty(
        connected, region, 0, 0,
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Float64}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Float64}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Float64}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZDistanceBinning(Δ=1.,max=1.,round_max=RoundUp)
    )
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
    for i ∈ eachindex(pos)
        corr.npart += 1
        ri = pos[i]
        vi = vec[i] - vmean
        if !ismissing(sv[i]) si = sv[i] - smean end
        modi = LinearAlgebra.norm(vec[i]) - modmean

        for j ∈ i:size(pos,1)
            dr = distance(corr.region,pos[j],ri)
            corr.npairs[dr] += 1
            vj = vec[j] - vmean
            corr.Cv[dr] += vi ⋅ vj
            if !ismissing(sv[i]) && !ismissing(sv[j])
                corr.npairs_nz[dr] += 1
                sj = sv[j] - smean
                corr.Cs[dr] += si ⋅ sj
            end
            modj = LinearAlgebra.norm(vec[j]) - modmean
            corr.Cmod[dr] += modi * modj
        end
    end

end

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
    return r,Cv,Cs,Cmod
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
    npart::Int
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
function density_correlation(region::PeriodicRegion,pos::ConfigurationT;Δr,rmax=nothing)
    if isnothing(rmax)
        rmax = 0.5 * volume(region)^(1/dimension(region))
    end
    dcorr = Density_correlation(
        region, size(pos,1),0,
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros)
    )
    return density_correlation!(dcorr,pos)
end

"""
    density_correlation!(dcorr::Density_correlation{R},pos) where R <: Region

Use configuration `pos` to add statistics to the `DensityCorrelation`
object `dcorr`.
"""
function density_correlation!(dcorr::Density_correlation{R},pos::ConfigurationT) where R <: PeriodicRegion
    @assert dcorr.npart==size(pos,1)
    dcorr.nconf += 1
    for i ∈ 1:size(pos,1)-1, j ∈ i+1:size(pos,1)
        r = sqrt(distancesq(dcorr.region,pos[i],pos[j]))
        dcorr.npr[r] += 2
    end
    dcorr.npr[0.] += dcorr.npart
    return dcorr
end

"""
    rdf(dcorr::Density_correlation{R})

Compute the radial distribution function ``g(r)`` and the correlation
integral ``C(r)``.  Both are returned as a tuple of `ZBinnedVector`.
"""
function rdf(dcorr::Density_correlation{R}) where R <: PeriodicRegion
    rmax = interval(dcorr.npr)[2]
    Δ = delta(dcorr.npr)
    gr = ZBinnedVector{Float64}(Δ=Δ,max=rmax,round_max=RoundUp,init=zeros)
    Cr = ZBinnedVector{Float64}(Δ=Δ,max=rmax,round_max=RoundUp,init=zeros)
    ρ = dcorr.npart / volume(dcorr.region)
    for (i,r) ∈ enumerate(range(dcorr.npr))
        vol = shell_volume(r,Δ,Val(dimension(dcorr.region)))
        if i==1 continue end
        gr[i] = dcorr.npr[i] / (ρ * vol * dcorr.npart * dcorr.nconf)
        Cr[i] = Cr[i-1] + dcorr.npr[i] / (ρ * dcorr.npart * dcorr.nconf)
   end
    Cr[2] += dcorr.npr[1] / (ρ * dcorr.centers[1] * dcorr.nconf)  # Because C(r) does not exclude self-correlation
    Cr ./= volume(dcorr.region)
    return gr,Cr
end

""""
    densdens(dcorr::Density_correlation{R<:Region})

Compute the density-density correlation function ``G(r)`` from a
`Density_correlation` object.
"""
function densdens(dcorr::Density_correlation{R}) where R <: PeriodicRegion
    rmax = interval(dcorr.npr)[2]
    Δ = delta(dcorr.npr)
    G = BinnedVector{Float64}(Δ=Δ,min=0.,max=rmax,round_max=RoundUp,init=zeros)
    for (i,r) ∈ enumerate(range(dcorr.npr))
        if i==1 continue end
        svol = shell_volume(r,Δ,Val(dimension(dcorr.region)))
        if i==2 
            G[1] = (dcorr.npr[1] + dcorr.npr[2]) / (svol * dcorr.nconf)
        else
            G[i-1] = dcorr.npr[i] / (svol * dcorr.nconf)
        end
    end
    G ./= volume(dcorr.region)
    return G
end    

function density_correlation(region::NonPeriodicRegion,pos::ConfigurationT;Δr,rmax=nothing)
    if isnothing(rmax)
        rmax = region.Lx-region.x0
    end
    dcorr = Density_correlation(
        region, size(pos,1),0,
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros),
        ZBinnedVector{Int}(Δ=Δr,max=rmax,round_max=RoundUp,init=zeros)
    )
    return density_correlation!(dcorr,pos)
end

function density_correlation!(dcorr::Density_correlation{<:NonPeriodicRegion},pos::ConfigurationT)
    @assert dcorr.npart==size(pos,1)
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

shell_volume(r,Δ,::Val{2}) = π * ( (r+Δ/2)^2 - (r-Δ/2)^2 )
shell_volume(r,Δ,::Val{3}) = (4π/3) * ( (r+Δ/2)^3 - (r-Δ/2)^3 )

function rdf(dcorr::Density_correlation{R}) where R <: NonPeriodicRegion
    rmax = interval(dcorr.npr)[2]
    Δ = delta(dcorr.npr)
    gr = ZBinnedVector{Float64}(Δ=Δ,max=rmax,round_max=RoundUp,init=zeros)
    Cr = ZBinnedVector{Float64}(Δ=Δ,max=rmax,round_max=RoundUp,init=zeros)
    ρ = dcorr.npart / volume(dcorr.region)
    for (i,r) ∈ enumerate(range(dcorr.npr))
        vol = shell_volume(r,Δ,Val(dimension(dcorr.region)))
        if i==1 continue end
        gr[i] = dcorr.npr[i] / (ρ * vol * dcorr.centers[i])
        Cr[i] = Cr[i-1] + dcorr.npr[i] / (ρ * dcorr.centers[i])
    end
    Cr[2] += dcorr.npr[1] / (ρ * dcorr.centers[1])  # Because C(r) does not exclude self-correlation
    Cr ./= volume(dcorr.region)
    return gr,Cr
end

# Maybe this has a problem with r=0?
function densdens(dcorr::Density_correlation{R}) where R <: NonPeriodicRegion
    rmax = interval(dcorr.npr)[2]
    Δ = delta(dcorr.npr)
    G = BinnedVector{Float64}(Δ=Δ,min=0.,max=rmax,round_max=RoundUp,init=zeros)
    ρ = dcorr.npart / volume(dcorr.region)
    for (i,r) ∈ enumerate(range(dcorr.npr))
        if i==1 continue end
        svol = shell_volume(r,Δ,Val(dimension(dcorr.region)))
        if i==2 
            G[1] = ρ * (dcorr.npr[1] + dcorr.npr[2]) /
                (svol * (dcorr.centers[1] + dcorr.centers[2]) * dcorr.nconf)
        else
            G[i-1] = ρ *  dcorr.npr[i] / (svol * dcorr.centers[i] * dcorr.nconf)
        end
    end
    return G
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
