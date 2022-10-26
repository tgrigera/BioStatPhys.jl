# timecorr.jl -- computation of time correlations
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

import AbstractFFTs
import FFTW

"""
    time_correlation(A; connected=true, normalized=false, nt=nothing, i0=nothing, Amean=nothing)

Compute the time (auto-)correlation function for signal `A`.  Returns
a vector `C[1:nt]`

 - `A`: the time signal, assumed to be sampled at evenly-spaced
   intervals.  If `i0==nothing`, it must be a `Vector` (real or
   complex), which is a single realisation or measurment of the random
   signal, otherwise it must be a matrix, where columns represent times
   and rows are different realizations of the random signal.

 - If `i0==nothing` (default), then `A` is assumed stationary (or
   time-translation invariant, TTI), and an estimate employing a
   single sequence is computed.  Otherwise, it must be an integer in
   the range `1<=i0<=size(A,2)` and is interpreted as an index for the
   desired reference time.

 - `connected`: if true, compute the connected (i.e. subtracting the
   time average) correlation.

 - `normalized`: if true, return the result normalized by `C[1]`.  Not
   recommended if non-TTI.

 - `nt`: the maxium time difference to compute in the TTI case,
   otherwise ignored.  Default `size(A,1)÷2`.

 - `Amean`: if `connected` is requested, then the signal mean can be
   given if known, otherwise it will be computed.

In the TTI case, an FFT-based implementation is used.

In the non-TTI, connected, case, it is probably better to use the
covariance function of the `Statistics` package, as `cov(A,dims=1)`.
"""
function time_correlation(X;connected=true,normalized=false, i0=nothing,nt=-1,Amean=nothing)
    if connected ave = isnothing(Amean) ? Statistics.mean(X) : Amean
    else ave = 0.
    end
    N=size(X,1)
    if nt<0 nt=N÷2 end
    if nt>N nt=N end
    if isnothing(i0) time_correlation_tti_fft(X,connected=connected,normalized=normalized,
                                              nt=nt,Xmean=ave)
    else time_correlation_tw_direct(X,i0=i0,connected=connected,normalized=normalized,
                                    Xmean=ave)
    end
end

"Compute the time correlation of real data for the TTI case with the FFT algorithm."
function time_correlation_tti_fft(X::Vector{<:Real};connected,normalized,Xmean,nt)
    N=size(X,1)

    # Subtract mean and pad data with 0s
    pdata= vcat(X,zero(X))
    if connected pdata[1:N].-=Xmean end

    # Compute IFFT(FFT(data)*conj(FFT(data)))
    fdata=AbstractFFTs.rfft(pdata)
    fdata.*=conj.(fdata)
    pdata=AbstractFFTs.irfft(fdata,size(pdata,1))
    # Normalize
    pdata[1]/=N
    if normalized
        for i=1:nt-1 pdata[i+1]/=pdata[1]*(N-i) end
        pdata[1]=1
    else
        for i=1:nt-1 pdata[i+1]/=N-i end
    end
    return pdata[1:nt]
end

"Compute the time correlation of complex data for the TTI case with the FFT algorithm."
function time_correlation_tti_fft(X::Vector{<:Complex};connected,normalized,Xmean,nt)
    N=size(X,1)

    # Subtract mean and pad data with 0s
    pdata= vcat(X,zero(X))
    if connected pdata[1:N].-=Xmean end

    fdata=AbstractFFTs.fft(pdata)
    fdata.*=conj.(fdata)
    pdata=AbstractFFTs.ifft(fdata)
    pdata[1]/=N
    if normalized
        for i=1:nt-1 pdata[i+1]/=pdata[1]*(N-i) end
        pdata[1]=1
    else
        for i=1:nt-1 pdata[i+1]/=N-i end
    end
    return pdata[1:nt]
end


"""
    time_correlation_tti_direct(X;connected,normalized,Xmean,nt)

Compute time correlation in the TTI case (i.e. averaging over time
origin) by the direct ``O(N^2)`` method.  Provided mostly as a check of
the algoritm that uses FFT.  Not exported.
"""
function time_correlation_tti_direct(X;connected,normalized,Xmean,nt)
    if !connected Xmean=0. end
    n=size(X,1)
    C=zeros(eltype(X),nt)

    for i=1:n
        xi = X[i]-Xmean
        for j=i:min(n,nt+i-1)
            C[j-i+1] += conj(xi) * (X[j]-Xmean)
        end
    end
    
    f = normalized ? C[1]/n : 1.
    for i=1:nt
        C[i]/=f*(n-i+1)
    end

    return C
end

"""
    time_correlation_tw_direct(X::Matrix{<:Number};i0,connected,
    normalized,Xmean)

Compute time correlation for fixed time origin (tw), specified by `i0`.

`X` must be a numerical matrix, with each row representing one sample
of the signal (i.e. columns are times, rows are different experiments).

`Xmean` is the mean of X at each time, e.g. `Xmean=Statistics.mean(X,dims=1)`.
If `nothing`, it will be computed.
"""
function time_correlation_tw_direct(X;i0,connected,normalized,Xmean)
    n=size(X,2)
    nt=n-i0+1
    C=zeros(eltype(X),nt)

    if connected
        if isnothing(Xmean) Xmean=Statistics.mean(X,dims=1) end
        for i=0:nt-1
            C[i+1] = Statistics.mean(conj.(X[:,i0].-Xmean[i0]) .* (X[:,i0+i].-Xmean[i0+i]))
        end
    else
        for i=0:nt-1
            C[i+1] = Statistics.mean(conj.(X[:,i0]) .* (X[:,i0+i]))
        end
    end
    
    if normalized C./=C[1] end

    return C
end

import Roots

"""
    correlation_time_spectral(C,Δt)

Compute the spectral correlation time ``\\tau`` from the connected time
correlation `C` (e.g. as computed by [`time_correlation`](@ref)).  The
scalar `Δt` is the time step for the sampling of `C`.

The exact definition ``\\tau`` is

```math
 \\int_0^\\infty \\!\\!dt \\, \\frac{C_c(t)}{C_c(0)} \\frac{\\sin t/\\tau}{t} = \\frac{\\pi}{4}
```
"""
function correlation_time_spectral(C,Δt)
    t0=Δt/1000
    t1=Δt*size(C,1)-1
    in0=intHH(t0,C,Δt)
    in1=intHH(t1,C,Δt)
    if in0*in1>0
        throw("Problem in correlation_time_spectral: failed to bracket: int($t0)=$in0,  int($t1)=$in1")
    end
    return Roots.find_zero(x->intHH(x,C,Δt),(t0,t1),Roots.A42())
end

"Compute ``\\int_0^T C(t) \\frac{ \\sin(t/tau) / t} \\, dt - \\pi/4`` "
function intHH(tau,C,Δt)
    I0::Float64=C[1]*Δt/tau
    for i=2:length(C)
        I0+= C[i]*sin(i*Δt/tau)/i
    end
    return I0/C[1]-π/4
end
