# BioStatPhys.jl -- module definition
#
# Copyright (C) 2022 Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
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

module BioStatPhys

import Statistics,LinearAlgebra

import Base.push!,Base.show
export MeanVar,mean,var,WMeanVar
include("./stat/meanvar.jl")

import Base.range
export BinnedVector,bin,binc,interval,delta,nbins
include("./tool/BinnedVector.jl")

export distance_binning,DistanceBinning
include("./tool/binnings.jl")

export Histogram,outliers,area,binc,prob,median,counts
include("./stat/Histogram.jl")

export time_correlation
include("./stat/timecorr.jl")

export space_correlation, correlation_length_r0
include("./stat/spacecorr.jl")

end
