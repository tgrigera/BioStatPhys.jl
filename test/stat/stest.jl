# stest.jl -- tests for MeanVar and time correlations
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For details, see the file LICENSE in the root directory, or
# check <https://www.gnu.org/licenses/>.

function test_MeanVar(file::String)
    data=readdlm(file,comments=true)
    mv=MeanVar()
    map(x->push!(mv,x),data)
    return mean(mv),var(mv)
end

function test_WMeanVar(file::String;equalweights=false)
    data=readdlm(file,comments=true)
    mv=WMeanVar()
    if equalweights
        map(x->push!(mv,x,1),data)
    else
        for i=1:size(data,1)
            push!(mv,data[i,1],data[i,2])
        end
    end
    return mean(mv),var(mv)
end

# Load dictionary with the results for three NIST datasets 

test_MeanVar_dict=Dict{String,NamedTuple}(
    "./stat/NumAcc1.dat"=>(mean=10000002., sd=1.),
    "./stat/NumAcc2.dat"=>(mean=1.2, sd=0.1),
    "./stat/NumAcc4.dat"=>(mean=10000000.2, sd=0.1) )

test_WMeanVar_dict=Dict{String,NamedTuple}(
    "./stat/NumAcc1.dat"=>(mean=10000002., sd=1., ew=true),
    "./stat/NumAcc2.dat"=>(mean=1.2, sd=0.1, ew=true),
    "./stat/NumAcc4.dat"=>(mean=10000000.2, sd=0.1, ew=true),
    "./stat/wtest1.dat"=>(mean=80.30555555555556, sd=15.422269466672065, ew=false)
)

function test_tcorr(X)
    C1=BioStatPhys.time_correlation_tti_direct(X,connected=true)
    C2=BioStatPhys.time_correlation_tti_fft(X,connected=true)
    @test C1 ≈ C2

    Xc=X .+ im.*rand(size(X))
    C1=BioStatPhys.time_correlation_tti_direct(Xc,connected=true)
    C2=BioStatPhys.time_correlation_tti_fft(Xc,connected=true)
    @test C1 ≈ C2
end
