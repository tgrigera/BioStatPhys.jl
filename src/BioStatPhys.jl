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

import Statistics
using  LinearAlgebra

ConfigurationT = AbstractVector{T} where T<:AbstractVector{W} where W<:Number

import Base.push!,Base.show
export MeanVar,mean,var,WMeanVar
include("./stat/meanvar.jl")

export BinnedVector,BBinnedVector,ZBinnedVector,bin,binc,interval,delta,nbins
include("./tool/BinnedVector.jl")

include("./tool/region.jl")
export Region, PeriodicRegion, NonPeriodicRegion, HyperCube, Rectangle, Cube
export PeriodicHyperCube, PeriodicRectangle, PeriodicCube
export dimension, volume, linear_size, dborder, distance, distancesq, fold!

include("./tool/binnings.jl")
export distance_binning,DistanceBinning,ZDistanceBinning

export Histogram,outliers,area,binc,prob,median,counts
include("./stat/Histogram.jl")

export time_correlation, correlation_time_spectral
include("./stat/timecorr.jl")

export GeoAve, get_mean
include("./stat/GeoAve.jl")

export space_correlation, space_correlation!, correlations, correlation_length_r0
export density_correlation, density_correlation!, rdf
include("./stat/spacecorr.jl")

end
