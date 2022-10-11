# Histogram.jl -- Histogram type, based on BinnedVector
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

mutable struct Histogram
    ndata::BigInt
    counts::BinnedVector{BigInt}
end

Histogram(nbins;min::AbstractFloat,max::AbstractFloat) =
    Histogram(0,BinnedVector{BigInt}(nbins,min=min,max=max,init=zeros))

function push!(his::Histogram,datum::AbstractFloat)
    his.ndata+=1
    his.counts[datum]+=1
end

outliers(his::Histogram)=his.counts[0]+his.counts[-1]

area(his::Histogram)=(his.ndata-outliers(his))/his.ndata;

binc(his::Histogram,bin)=binc(his.counts,bin)

prob(his::Histogram,bin)= his.counts[bin]/(his.ndata*delta(his.counts))

function median(his::Histogram)
    cc=his.counts[0]
    m=0
    while cc<his.ndata/2
        m+=1
        cc+=his.counts[m]
    end
    return binc(his,m)
end

function prob(his::Histogram)
    min,max=interval(his.counts)
    x=collect(min+delta(his.counts)/2:delta(his.counts):max)
    f = 1/(his.ndata*delta(his.counts))
    return x,f*his.counts
end

function counts(his::Histogram)
    x=collect(his.counts.min+his.counts.delta/2:his.counts.delta:his.counts.max)
    return x,his.counts
end
