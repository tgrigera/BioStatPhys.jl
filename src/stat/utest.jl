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
using Base: return_types

"""
    compute_UC_UE(control,experiment)

Given two vectors with the results of measurements ``X_C`` under
control conditions and measurements ``X_E`` under experimental
conditions, compute the Mann-Whitney statistic ``U_E``, namely the
number of times that ``X_E`` exceeds ``X_C`` (considering all possible
pairs of control and experimental measurements).

"""
function compute_UC_UE(control,experiment)
    UE = 0
    UC = 0
    for XE ∈ experiment, XC ∈ control 
	if XE>XC UE += 1
        elseif XC>XE UC +=1
        end
    end
    return (;UC, UE)
end

"""
    Bayesian_U_test_MC(UE,nC,nE;bins=200,rng=Random.default_rng())

Perform a Bayesian-style Mann-Whitney test, computing the likelihood
``P(U_E|\\Theta)`` by a Monte Carlo method, dividing the unit interval
in `bins` bins and using the provided random number generator `rng`.
For this reason this function is rather slow if `nC` and `nE` are
large.  If both numbers are larger than about 15, then an analytical
approximation for likelihood is recommended (see
`Bayesian_U_test_beta`).

Arguments:

 - `UE`: measured value of the ``U_E`` statistic, i.e. the number of
    times that experiment exceeds control (considering all pairs).
    Can be computed with `compute_UC_UE`.

 - `nE`: number of "experimerimental" observations ``n_E``.

 - `nC`: number of "control" observations ``n_C``.

 - `bins`: number of bins to use to generate the likelihood with Monte
    Carlo.

 - `rng`: A random number generator to use for the Monte Carlo computation.

Returns a named tuple with

 - `PΘ_U`: the posterior ``P(\\Theta | U_E)`` (a vector of reals),

 - `Θbins`: a range giving the ``\\Theta`` values of the elements in the `PΘ_U` vector.

The method is by R. A. Chechile,

 - Richard A. Chechile, _Bayesian Statistics for Experimental Scientists,_ MIT Press, Cambdrige, MA (2020), Chapter 7,

which is also inspiration for the code. 

"""
function Bayesian_U_test_MC(UE,nC,nE;bins=200,rng=Random.default_rng())

    Δ = 1/bins
    Θbins = range(start=Δ/2,stop=1-Δ/2,length=bins)
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
        for _ in 1:100000   # Compute the likelihood of UE for this value of θ using MonteCarlo
            Uz = 0
            for L in 1:nE
                XE[L] = rand(rng, rexp)
            end
            for L in 1:nC
                XC[L] = rand(rng, uexp)
            end
            _, Uz = compute_UC_UE(XC,XE)
            if Uz == UE
                PΘ_U[j] = PΘ_U[j] + 1
            end
        end
    end
    PΘ_U ./= sum(PΘ_U) # Normalize.  Since the prior is uniform, this is now the posterior P(Θ|U_E)

    return (;Θbins, PΘ_U) 
end

"""
    Bayesian_U_test_beta(UE,nC,nE;prior_α=1.,prior_β=1.)

Perform a Bayesian-style Mann-Whitney test, using an analytical
approximation for the likelihood ``P(U_E|\\Theta)``.

Arguments:

 - `UE`: measured value of the ``U_E`` statistic, i.e. the number of
    times that experiment exceeds control (considering all pairs).
    Can be computed with `compute_UC_UE`.

 - `nE`: number of "experimerimental" observations ``n_E``.

 - `nC`: number of "control" observations ``n_C``.

 - `prior_α`, `prior_β` (optional): Values of the ``\alpha`` and
   ``\\beta`` parameters of the ``\beta`` distribution used as prior.
   Default is to take a flat (uniform) prior over ``[0,1]``,
   i.e. ``\alpha=\beta=1``.

Returns a named tuple with

 - `PΘ_U`: the posterior ``P(\\Theta | U_E)`` (a `Beta` object from
   `Distributions`).

Reference for method and inspiration for code:
 
 - Richard A. Chechile, _Bayesian Statistics for Experimental
   Scientists,_ MIT Press, Cambdrige, MA (2020), Chapter 7.

"""
function Bayesian_U_test_beta(UE,nC,nE;prior_α=1.,prior_β=1.)
    # Posterior distribution is also Beta, with the parameters given below.  These follow from
    # modelling the likelihood with a modified binomial distribution [Chechile]
    xs = UE / (nE*nC)
    if xs >= .5  x=xs  else x=1−xs end
    nH = (2 * nE * nC) / (nE + nC)
    y5 = (nH^1.1489) / (0.4972 + (nH^1.1489))
    w4 = 0.8 − 1 / (1 + 1.833 * nH)
    w3 = 0.6 − 1 / (1 + 2.111 * nH)
    w2 = 0.4 − 1 / (1 + 2.520 * nH)
    w1 = 0.2 − 1 / (1 + 4.813 * nH)
    y4 = y5 * w4 + (1 − w4) / 2
    y3 = y5 * w3 + (1 − w3) / 2
    y2 = y5 * w2 + (1 − w2) / 2
    y1 = y5 * w1 + (1 − w1) / 2
    Y = [0.5, y1, y2, y3, y4, y5]
    La0 = 252 − 1627 * x + (12500 * x^2 - 15875 * x^3 + 10000 * x^4 - 2500 * x^5) / 3
    La1 = -1050 + 42775 * x / 6 - 38075 * x^2 / 2 + (75125 * x^3 − 48750 * x^4 + 12500 * x^5) / 3
    La2 = 1800 - 12650 * x + (104800 * x^2 - 142250 * x^3 + 95000 * x^4 - 25000 * x^5) / 3
    La3 = -1575 + 11350 * x + (-96575 * x^2 + 134750 * x^3 − 92500 * x^4 + 25000 * x^5) / 3
    La4 = 700 + 14900 * x^2 + 15000 * x^4 - (15425 * x + 63875 * x^3 + 12500 * x^5) / 3
    La5 = -126 + 1879 * x / 2 + (-16625 * x^2 / 2 + 12125 * x^3 - 8750 * x^4 + 2500 * x^5) / 3
    LA = [La0, La1, La2, La3, La4, La5]
    Ωhat = sum(Y .* LA)
    ab = nH * (1.028 + 0.75 * x) + 2
    α = Ωhat * ab
    β = (1 - Ωhat) * ab
    if xs < 0.5
        α, β = β, α
    end
    return (;PΘ_U = Beta(prior_α+α-1, prior_β+β-1))

end

"""
    Bayesian_U_test(XC::Vector{<:Real},XE::Vector{<:Real}; method=nothing,
                    rng=Random.default_rng())

Perform a Bayesian-style Mann-Whitney test.  This test assumes a real
variable is measured under "control" and "experimental" conditions
several times in two groups of identical experiments.  The aim is
to establish to what degree the variable measured under experimental
conditions is likely to exceed the value measured under control
conditions.

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
given value of ``\\Theta``.  This is computed either through a Monte
Carlo method  or using an analytical approximation (modified
binomial distribution); see the `method` parameter below.

Arguments:

 - `XC`: Vector of points measured under _control_ conditions

 - `XE`: Vector of points measured under _experimental_ conditions

 - `method`: specify `:MC` for the Monte Carlo method, or
   `:analytical` for the analytical approximation to the likelihood.
   If `nothing` (default), the Monte Carlo method will be chosen for
   small samples (`length(XC)<=15` and `length(XE)<=15`).  Note that
   `:MC` can be slow for large amounts of data.

 - `rng`: Optional random number generator, to be used if the Monte
   Carlo algorithm is chosen.

Returns: a tuple with

 - `UC` and `UE`: the values of the ``U_C`` and ``U_E`` statistics
   (the number of times that control exceeds experiment, or experiment
   exceeds control respectively).

 - `Θm`: The median of ``\\Theta`` according to the posterior ``P(\\Theta | U_E)``.

 - `Θc`: A tuple given the confidence or credibility interval, namely
   the interval that includes 95% of the probability around the median.

 - `PΘ_U`: the posterior ``P(\\Theta | U_E)`` (a vector if the method
   was Monte Carlo, a `Beta` object from `Distributions` if the
   analytical approximation was used).

 - `Θbins`: a range given the ``\\Theta`` values corresponding to the
   `PΘ_U` vector.  It is returned but arbitrary if `PΘ_U` is a `Beta`
   object.

 - `PH1`: probability that ``\\Theta>1/2``.  This can be regarded as
   the probability of ``H_1``, i.e. the hypothesis that the effect of
   the experimental conditions is to increase ``X``.  Especially for
   small samples, it is not enough as single measure of the effect.
   It is only strong evidence if greater than 95%, with 80% being only
   moderate evidence.  Use in conjunction with the credibility
   interval.

The method, as well as inspiration for the code, come from

 - Richard A. Chechile, _Bayesian Statistics for Experimental
   Scientists,_ MIT Press, Cambdrige, MA (2020), Chapter 7.

"""
function Bayesian_U_test(XC::Vector{<:Real},XE::Vector{<:Real}; method=nothing,
    rng=Random.default_rng())

    UC, UE = compute_UC_UE(XC,XE)
    nC, nE = length(XC), length(XE)
    if method === nothing
        method = length(XC)<=15 || length(XE)<=15 ? :MC : :analytical
    end

    if method == :MC

        r = Bayesian_U_test_MC(UE,nC,nE,bins=200,rng=rng)
        FΘ_U = cumsum(r.PΘ_U)  # cumulative distribution
        i = findlast(x->x<=0.5,FΘ_U)
        Θbins = r.Θbins
        Θm = (Θbins[i]+Θbins[i])/2
        i0 = findlast(x->x<=0.025,FΘ_U)
        i1 = findlast(x->x<=0.975,FΘ_U)
        Θc = ( (Θbins[i0]+Θbins[i0])/2, (Θbins[i1]+Θbins[i1])/2 )
        PH1 = 1−FΘ_U[end÷2]  # probability for E>C

    else

        r = Bayesian_U_test_beta(UE,nC,nE)
        Δ = 1/200
        Θbins = range(start=Δ/2,stop=1-Δ/2,length=200)
        Θm = Distributions.median(r.PΘ_U)
        Θc = (quantile(r.PΘ_U, 0.025), quantile(r.PΘ_U, 0.975))
        PH1 = 1.0 - cdf(r.PΘ_U, 0.5)

    end

    return (;UC, UE, Θm, Θc, PH1, PΘ_U=r.PΘ_U, Θbins)
end
