# runtests.jl -- run all tests for BioStatPhys.jl
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

using BioStatPhys
using DelimitedFiles
using Test

include("./stat/stest.jl")

@testset "BioStatPhys/stat" begin

    for (file,res) in test_MeanVar_dict
        m,v = test_MeanVar(file)
        @test m≈res.mean && sqrt(v)≈res.sd
    end

    for (file,res) in test_WMeanVar_dict
        m,v = test_WMeanVar(file,equalweights=res.ew)
        @test m≈res.mean && sqrt(v)≈res.sd
    end

end
