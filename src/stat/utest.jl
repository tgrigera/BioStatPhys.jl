# utest.jl -- Mann-Whitney test
#
# Copyright (C) 2026 by Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
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

using Random
using Distributions

"""
    compute_UE(control,experiment)

Given two vectors with the results of measurements ``X_C`` under
control conditions and measurements ``X_E`` under experimental
conditions, compute the Mann-Whitney statistic ``U_E``, namely the
number of times that ``X_E`` exceeds ``X_C`` (considering all possible
pairs of control and experimental measurements).

"""
function compute_UE(control,experiment)
    UE = 0
    for XE ∈ experiment, XC ∈ control 
	if XE>XC UE += 1 end
    end
    return UE
end

"""
    Bayesian_U_test_small(UE,nE,nC;bins=200,rng=Random.default_rng())

Perform a Bayesian-style Mann-Whitney test for small sample size.
This test assumes a real variable is measured under "control" and
"experimental" conditions several times in identical experiments in
each condition.  The aim is to establish to what degree the variable
measured under experimental conditions is likely to exceed the value
measured under control conditions.

 - `UE`: measured value of the ``U_E`` statistic, i.e. the number of
    times that experiment exceeds control (considering all pairs).
    Can be computed with `compute_UE`.

 - `nE`: number of "experimerimental" observations ``n_E``.

 - `nC`: number of "control" observations ``n_C``.

The quantity of interest of the test is the probability that a random
experimental observation ``X_E`` will be larger than a random control
observation ``X_C``,
```math
    \\Theta = P(X_E > X_C).
```
If ``\\Theta>1/2`` then the experimental conditions are likely to
yield higher values of ``X`` than the control conditions.


In the Bayesian spirit, the methods seeks to compute the posterior
probability ``P(\\Theta | U_E)``, i.e. the probability distribution of
``\\Theta``, given the observed value of ``U_E``.  The prior
``P(\\Theta)`` (i.e. the knowledge about ``\\Theta`` before performing
the experiment) is assumed uniform in ``[0,1]``(i.e. no knowledge of
``\\Theta``).  The method also needs the _likelihood_
``P(U_E|\\Theta)``, i.e. the probability of measuring ``U_E`` for a
given value of ``\\Theta``.  This is computed in this function through
Monte Carlo, dividing the unit interval in `bins` bins and using the
provided random number generator `rng`.  For this reason this function
is rather slow if `nC` and `nE` are large.  If both numbers are larger
than about 15, then another method is recommended (see
`Bayesian_U_test_large`.

Returns a tuple with

 - `PΘ_U`: the posterior ``P(\\Theta | U_E)``,

 - `Θav`: average of ``\\Theta`` under the computed posterior,

 - `FΘ_U`: cumulative posterior probability for ``\\Theta``,

 - `pH1`: probability that ``\\Theta>1/2``, i.e. ``1-F_\\Theta(1/2)``.
   This can be regarded as the probability of ``H_1``, i.e. the
   hypothesis that the effect of the experimental conditions is to
   increase ``X``.

Method and code after

 - Richard A. Chechile, _Bayesian Statistics for Experimental Scientists,_ MIT Press, Cambdrige, MA (2020).  Chapter 7.

"""
function Bayesian_U_test_small(UE,nE,nC;bins=200,rng=Random.default_rng())

    Θbins = range(start=.0025,stop=.9975,length=bins)
    PΘ_U = zeros(Float64,bins)
    XE=zeros(Float64,nE)
    XC=zeros(Float64,nC)
    uexp = Exponential(1)

    # Compute the likelihood P(U|Θ).  The prior is taken flat (uniform in [0,1]).
    # It is stored in the same vector that will hold the posterior
    for j ∈ 1:bins
        Θ = Θbins[j]
        kΘ = Θ / (1 − Θ)   # this is scale, not rate, i.e. 1/rate
        rexp = Exponential(kΘ)
        for _ in 1:10000   # Compute the likelihood of UE for this value of θ using MonteCarlo
            Uz = 0
            for L in 1:nE
                XE[L] = rand(rng, rexp)
            end
            for L in 1:nC
                XC[L] = rand(rng, uexp)
            end
            Uz = compute_UE(XC,XE)
            if Uz == UE
                PΘ_U[j] = PΘ_U[j] + 1
            end
        end
    end
    PΘ_U ./= sum(PΘ_U) # Normalize.  Since the prior is uniform, this is now the posterior P(Θ|U_E)
    Θav = sum(PΘ_U.*Θbins)  # posterior mean
    FΘ_U = cumsum(PΘ_U)  # cumulative distribution
    pH1 = 1−FΘ_U[bins÷2]  # probability for E>C

    return (;PΘ_U, Θav, FΘ_U, pH1)
end

function Bayesian_U_test_small(XC::Vector{<:Real},XE::Vector{<:Real};
    bins=200,rng=Random.default_rng())
    UE = compute_UE(XC,XE)
    return Bayesian_U_test_small(UE,length(XE),length(XC),bins=bins,rng=rng)
end

function Bayesian_U_test_large(treatment,control)
    UE = 0.0 # Treatment beats Control
    UC = 0.0 # Control beats Treatment

    for t in treatment, c in control
        if t > c
            UE += 1.0
        elseif t < c
            UC += 1.0
        else # Split ties equally
            UE += 0.5
            UC += 0.5
        end
    end

    # 3. Formulate the Bayesian update using Beta(α, β)
    # Prior: Beta(1,1) represents a perfectly uniform/flat initial belief
    prior_α, prior_β = 1.0, 1.0
    posterior = Beta(prior_α + E, prior_β + C)

    # 4. Extract metrics
    θ_median = median(posterior)
    # 95% Credible Interval (Highest Density or Equal-Tailed)
    ci_lower = quantile(posterior, 0.025)
    ci_upper = quantile(posterior, 0.975)

    # Directional Test: Probability that treatment genuinely improves the metric (θ > 0.5)
    prob_improves = 1.0 - cdf(posterior, 0.5)

    return (; UE, UC, θ_median, ci_lower, ci_upper, prob_improves)
    #println("Estimated Effect Size (θ): ", round(θ_median, digits=3))
    #println("95% Credible Interval: [", round(ci_lower, digits=3), ", ", round(ci_upper, digits=3), "]")
    #println("Probability that Treatment improves quantity: ", round(prob_improves * 100, digits=1), "%")
end
