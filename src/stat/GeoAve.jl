# GeoAve.jl -- Running mean and variance with West's recursion
#
# Copyright (C) 2024 by Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
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
    GeoAve(;t0=0,wfactor=1.5,base=1.)

Return a `GeoAve` object.  This is used to takes pairs of numbers (x,y)
(through `push(G,x,y)`) and average together all points whose x-coordinate (assume
it is a time) falls within a window.  The first window starts
at `t_0` and is of length `base`, successive windows grow
geometrically by a factor `wfactor`.

The initial time `t0` is handled separately, and values corresponding
to `t0` are averaged without windowing.

After filling the object with points, use `get_mean` to obtain the
averages.  `wfactor=1` is supported and handled as a special case.
`wfactor<1` is not recommended.

"""
mutable struct GeoAve
    base    ::Float64
    t0      ::Float64
    wfactor ::Float64
    logwf   ::Float64
    read_fb ::Float64
    mean    ::Vector{MeanVar}

    GeoAve(;t0=0,wfactor=1.5,base=1.) =
        new(base,t0,wfactor,log(wfactor),(wfactor-1)/base,MeanVar[])
end

function Base.push!(gav::GeoAve,time::Number,e::Number)
    if time<gav.t0 return end
    if time==gav.t0 n=1
    else
      n = gav.wfactor == 1 ? floor(Int, (time-gav.t0)/gav.base) :
        floor(Int, log(gav.read_fb*(time-gav.t0)+1)/gav.logwf )
      n += 2
    end
    if n > length(gav.mean)
        ol =  length(gav.mean)
        resize!(gav.mean,n+10)
        for i =ol+1:length(gav.mean) gav.mean[i] = MeanVar() end
    end
    push!(gav.mean[n],e)
    return nothing
end


"""
    get_mean(gav::GeoAve)

Return a tuple `(t,X,v)` of arrays, where `X` holds the average
of all data over a time window, `v` is the respective variance,
and `t` gives the centre of the window (except for `t=gav.t0`, which
is not windowed).
"""
function get_mean(gav::GeoAve)
    time = Float64[]
    ave  = Float64[]
    vr   = Float64[]

    if !isempty(gav.mean[1])
        push!(time,gav.t0)
        push!(ave,mean(gav.mean[1]))
        push!(vr,var(gav.mean[1]))
    end
  
    t = gav.t0 + 0.5*gav.base
    Δt = 0.5*gav.base*(1+gav.wfactor)
    for i = 2:length(gav.mean)
        if !isempty(gav.mean[i])
            push!(time,t)
            push!(ave,mean(gav.mean[i]))
            push!(vr,var(gav.mean[i]))
        end
        t += Δt
        Δt *= gav.wfactor
    end
    return time, ave, vr
end
