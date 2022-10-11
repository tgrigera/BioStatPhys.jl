# histotest.jl -- tests for Histogram
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

import SpecialFunctions

gauss_prob(x)=0.5*SpecialFunctions.erf(x/sqrt(2))  # Gaussian probability from 0 to x

function test_histogram()
    his=Histogram(21,min=-1.,max=1.)
    rng=Random.MersenneTwister(12560)
    gauss=Distributions.Normal(0.,1.)
    for _=1:1000000
        push!(his,rand(rng,gauss))
    end
    @test isapprox(area(his),2*gauss_prob(1.),rtol=1e-3)
    pr=prob(his,11)  # centered at 0
    δ=delta(his.counts)
    @test isapprox(δ*pr, 2*gauss_prob(δ/2), atol=1e-3 )
    _,pr=prob(his)
    @test δ*sum(pr) ≈ area(his)
end
