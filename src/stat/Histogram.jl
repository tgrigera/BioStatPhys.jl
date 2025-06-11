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

"""
    struct Histogram

Type to build histograms, based on `BinnedVector`.  Create a
`Histogram` object with

    his=Histogram(nbis,min=interval_min,max=interval_max)

Add data with `push!` and access with `area`, `outliers`, `prob` and
`counts`.
"""
mutable struct Histogram
    ndata::BigInt
    counts::BinnedVector{BigInt}
end

Histogram(nbins;min::Number,max::Number) =
    Histogram(0,BinnedVector{BigInt}(nbins,min=min,max=max,init=zeros))

"""
    push!(his::Histogram,datum::AbstractFloat)

Add new data point to histogram `his`
"""
function push!(his::Histogram,datum::AbstractFloat)
    his.ndata+=1
    his.counts[datum]+=1
end

push!(his::Histogram,datum::Integer) = push!(his,Float64(datum))


"""
    outliers(his::Histogram)

Return number of points outside histogram interval.
"""
outliers(his::Histogram)=his.counts[0]+his.counts[-1]

"""
    area(his::Histogram)

Compute total area under histogram within the interval defined at
creation.  This interprets histogram as a probability distribution, so
that `area` is bounded by 1.  It can be less than 1 due to outlier
points.
"""
area(his::Histogram)=(his.ndata-outliers(his))/his.ndata;

"""
    binc(his::Histogram,bin)

Return position of center of bin `bin`.  No range check performed.
"""
binc(his::Histogram,bin)=binc(his.counts,bin)

"""
    counts(his::Histogram)

Return tuple `(x,counts)` where `x` is a vector with the position of
bin centers and `counts` is a vector of bin counts.
"""
function counts(his::Histogram)
    return collect(range(his.counts)),his.counts
end

"""
    prob(his::Histogram,bin)

Return probability density for `bin`, i.e. counts for `bin`
divided by total data points and by the bin width.
"""
prob(his::Histogram,bin)= his.counts[bin]/(his.ndata*delta(his.counts))

"""
    prob(his::Histogram)

Return tuple `(x,prob)` where `x` is a vector with the position of bin
centers and `prob` is a vector of probability densities.  `prob` is to
be interpreted as needing integration over an appropriate interval,
e.g. the sum of `prob` elements multiplied by the bin with is equal to
`area(his)`.
"""
function prob(his::Histogram)
    f = 1/(his.ndata*delta(his.counts))
    return collect(range(his.counts)),f*his.counts
end

"""
    median(his::Histogram)

Return an approximation to the median as the center of the lowest bin
`b` such that the sum of the low outliers plus the counts of the bins
up to `b` exceeds half the data points.
"""
function median(his::Histogram)
    cc=his.counts[0]
    m=0
    while cc<his.ndata/2
        m+=1
        cc+=his.counts[m]
    end
    return binc(his,m)
end
