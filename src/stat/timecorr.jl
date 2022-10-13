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

# function time_correlation(X;connected=true,normalized=false,Xmean=nothing,
#                           i0::Nothing,nt=nothing)
#     time_correlation_tti_fft(X,connected=connected,normalized=normalized,Xmean=Xmean,nt=nt)
# end

function time_correlation_tti_fft(X::Vector{<:Real};connected=true,
                                  normalized=false,Xmean=nothing,nt=-1)
    N=size(X,1)
    if nt<0 nt=N÷2 end
    if nt>N nt=N end
    ave = connected ? (isnothing(Xmean) ? Statistics.mean(X) : Xmean) : 0.

    # Subtract mean and pad data with 0s
    pdata= vcat(X,zero(X))
    if connected pdata[1:N].-=ave end

    fdata=AbstractFFTs.rfft(pdata)
    fdata.*=conj.(fdata)
    pdata=AbstractFFTs.irfft(fdata,size(pdata,1))
    if normalized
        pdata[1]/=N
        for i=1:nt-1 pdata[i+1]/=pdata[1]*(N-i) end
        pdata[1]=1
    else
        for i=0:nt-1 pdata[i+1]/=N-i end
    end
    return pdata[1:nt]
end

function time_correlation_tti_fft(X::Vector{<:Complex};connected=true
                                  ,normalized=false,Xmean=nothing,nt=-1)
    N=size(X,1)
    if nt<0 nt=N÷2 end
    if nt>N nt=N end
    ave = connected ? (isnothing(Xmean) ? Statistics.mean(X) : Xmean) : 0.

    # Subtract mean and pad data with 0s
    pdata= vcat(X,zero(X))
    if connected pdata[1:N].-=ave end

    fdata=AbstractFFTs.fft(pdata)
    fdata.*=conj.(fdata)
    pdata=AbstractFFTs.ifft(fdata)
    if normalized
        pdata[1]/=N
        for i=1:nt-1 pdata[i+1]/=pdata[1]*(N-i) end
        pdata[1]=1
    else
        for i=0:nt-1 pdata[i+1]/=N-i end
    end
    return pdata[1:nt]
end


"""
    time_correlation_tti_direct(X;connected=true,normalized=false,Xmean=nothing,nt=-1)

Compute time correlation in the TTI case (i.e. averaging over time
origin) by the direct ``O(N^2)`` method.  Provided mostly as a check of
the algoritm that uses FFT.  Not exported.
"""
function time_correlation_tti_direct(X;connected=true,normalized=false,Xmean=nothing,nt=-1)
    n=size(X,1)
    if nt<0 nt=n÷2 end
    if nt>n nt=n end
    C=zeros(eltype(X),nt)
    ave = connected ? (isnothing(Xmean) ? Statistics.mean(X) : Xmean) : 0.

    for i=1:n
        xi = X[i]-ave
        for j=i:min(n,nt+i-1)
            C[j-i+1] += conj(xi) * (X[j]-ave)
        end
    end
    
    f = normalized ? C[1]/n : 1.
    for i=1:nt
        C[i]/=f*(n-i+1)
    end

    return C
end

"""
    time_correlation_tw_direct(X;i0=1,connected=true,normalized=false,Xmean=nothing)

Compute time correlation for fixed time origin (tw), specified by `i0`.
"""
function time_correlation_tw_direct(X;i0=1,connected=true,normalized=false,Xmean=nothing)
    n=size(X,1)
    nt=n-i0
    C=zeros(Float64,nt)
    ave = connected ? (isnothing(Xmean) ? Statistics.mean(X) : Xmean) : 0.

    for i=0:nt-1
        C[i+1] += conj((X[i0]-ave)) * (X[i0+i]-ave)
    end
    
    if normalized C./=C[1] end

    return C
end
