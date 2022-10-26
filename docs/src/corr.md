# Correlation functions

The functions in this section compute estimates of several space and time correlation functions.  The precised definitions of the different quantities computed are summarised below.  For more details, caveats and discussion of the correlation functions we refer to the review article by T. S. Grigera[^1].


## Time correlations and correlation time

### Time correlation function

The _time (auto-)correlation function_ of a signal ``a(t)`` (which is a noisy quantity that can be found to take different values if it is measured several times, i.e. it is a stochastic process) is defined as
```math
  C(t_0,t) = \langle a^*(t_0) a(t_0+t) \rangle,
```
where ``a^*`` is the complex conjugate and brackets stands for a  history-average (an average over many realizations of the stochastic process, i.e. over many repetitions of the experiment, resetting the initial conditions).

The _connected time correlation_, or _auto-covariance_, is
```math
   C_c(t_0,t) = \left\langle \delta a^*(t_0) \delta a(t_0+t) \right\rangle = 
   C(t_0,t) - \left\langle a^*(t_0) \right\rangle \left\langle a(t_0+t) \right\rangle,
```
where ``\delta a(t)=a(t)-\langle a(t) \rangle``.

Assume the process ``a(t)`` is sampled uniformly in time with interval ``\Delta t``.  Calling ``N`` the number of time samples and ``M``the number of experiments (i.e. different measurements of ``a(t)``  after resetting the initial conditions), we represent these data as ``M`` sequences ``A_i^{(n)}`` with ``i=1,\ldots,N``,  ``n=1,\ldots,M``.  Then the following statistical estimators of the correlation functions can be computed.

In the general case,
```math
\begin{align*}
  \langle a(t_i) \rangle & \approx \overline{A_i} =\frac{1}{M} \sum_{n=1}^M A_i^{(n)}, \\
  C_c(t_i,t_k) & \approx \hat C^{(c)}_{i,k} =\frac{1}{M}\sum_n^M \delta A_i^{(n)}
  \delta A_{i+k}^{(n)}  , \qquad \delta A_i^{(n)} = A_i^{(n)} - \overline{A_i}.
\end{align*}
```

If the process ``a(t)`` is stationary, then the correlation functions time-translation invariant (i.e. independent of ``t_0``), and an estimate can be obtained with a single sample, averaging over the time origin:
```math
\begin{align*}
  \overline{A} &= \frac{1}{N}\sum_{i=1}^N A_i, \\
  \hat C^{(c)}_k &= \frac{1}{N-k} \sum_{j=1}^{N-k} \delta A_j \delta A_{j+k},
  \qquad
  \delta A_j = A_j - \overline{A}.
\end{align*}
```
Note that the estimate obtained this way is noisier the larger the value of ``k``.

If the process is not stationary, then the first estimator, which needs several samples, must be used.  Although we implement this here, note that in the connected case the covariance function `cov` from the `Statistics` package can be used instead and is more convenient.  The present routine is more useful for stationary case.

Refer to the review article[^1] for more details.

### Correlation time

The _correlation time_ is a time scale that measures how separated in time two measurements must be for them to be significantly decorrelated.  A good, though abstract, definition is
```math
  \tau = \lim_{t\to\infty} \frac{t}{-\log C_c(t)}.
```
which picks the slowest (i.e. longest-ranged) exponential decay rate (cf [Correlation length](@ref)).

Since this definition is not directly applicable to finite data, several practical alternatives have been proposed[^1].  We mention here only those that are currently implemented in `BioStatPhys`.

#### Spectral correlation time

The spectral correlation time ``\tau_S`` is defined as the inverse of the frequency ``\omega_0`` such that half the spectral content of the Fourier transform of the connected correlation ``C_c(t)`` is contained in the interval ``[-\omega_0,\omega_0]``.  Translated to the time domain, the definition is
```math
 \int_0^\infty \!\!dt \, \frac{C_c(t)}{C_c(0)} \frac{\sin t/\tau}{t} = \frac{\pi}{4}
```
This is computed by [`correlation_time_spectral`](@ref).

### API

The stationary (or TTI) estimate for a single real or complex sequence is implemented with an algorithm that uses fast-Fourier-transforms, giving ``O(N\log N)`` performance.

```@docs
time_correlation
correlation_time_spectral
```


## Space correlations and correlation length

The space correlation functions of a space-dependent quantity ``a(\mathbf{r})`` are defined as
```math
\begin{align*}
   C(\mathbf{r}) & = \left\langle a(\mathbf{r_0}) a(\mathbf{r_0}+\mathbf{r}) \right\rangle,\\
   C_c(\mathbf{r}) & = \left\langle \delta a(\mathbf{r_0}) \delta a(\mathbf{r_0}+\mathbf{r}) \right\rangle,
\end{align*}
```
where ``\delta a(\mathbf{r}) = a(\mathbf{r}) - \langle a\rangle`` and we are assuming the system is homogeneous (i.e. translation-invariant).  ``C_c(\mathbf{r})`` is called _connected_ correlation in the physics literature, or _auto covariance_ in the mathematical statistics literature.

The [`space_correlation`](@ref) function computes the estimate of ``C(r)`` for an isotropic discrete system,
```math
\begin{align*}
  \hat C_c(r) & = \frac{\sum_{ij} \delta a_i \delta a_j
    \delta(r-r_{ij}) } {\sum_{kl} \delta(r-r_{kl})}.
\end{align*}
```
The average ``\langle a\rangle`` can be estimated with space average or phase average (see below and the review[^1] for details).

### Correlation length

The _correlation length_ is a length scale that measures how far apart two points in space must be taken for them to be significantly decorrelated.  A good, though abstract, definition is
```math
  \xi = \lim_{r\to\infty} \frac{r}{-\log C_c(r)},
```
which picks the slowest (i.e. longest-ranged) exponential decay rate (cf. [Correlation time](@ref)).

Since this definition is not directly applicable to finite data, several alternatives have been proposed, more suited to experimental or numerical determination but respecting the functional dependence with control parameters[^1].  We mention here only those that are currently implemented in `BioStatPhys`.

#### ``r_0``

If the connected correlation has been computed with space averaging, ``C_c(r)`` will have at least one zero, and the location ``r_0`` of the first of these can be used as a proxy of ``\xi``.  We stress that ``r_0`` is __not__ a correlation length, but a useful scale in the case the correlation can be measured for different system sizes ``L`` (or in observation windows of different size[^2]).  ``r_0`` scales with size as[^3]
```math
  \begin{align*}
    r_0 &\sim \xi \log(L/\xi), & L\gg\xi, \\
    r_0 & \sim L, & L\ll \xi.
  \end{align*}
```
The second situation will always be the case in scale-free systems where ``\xi=\infty``.  See [`correlation_length_r0`](@ref).


### API

```@docs
space_correlation
correlation_length_r0
```

## References

[^1]:T. S. Grigera, Correlation functions as a tool to study
     collective behaviour phenomena in biological
     systems. _J. Phys. Complex._ __2,__ 045016 (2021).
     [[DOI](https://doi.org/10.1088/2632-072X/ac2b06)].  This review
     discusses the correlation functions and estimators computed by
     the routines documented here and is the general reference for
     this section.

[^2]:D. A. Martin, T. L. Ribeiro, S. A. Cannas, T. S. Grigera,
     D. Plenz, and D. R. Chialvo. Box scaling as a proxy of finite
     size correlations. _Sci Rep_ __11,__ 15937. (2021)
     [[DOI](http://dx.doi.org/10.1038/s41598-021-95595-2)]

[^3]:A. Cavagna, I. Giardina, and T. S.  Grigera. The physics of
     flocking: Correlation as a compass from experiments to
     theory. _Physics Reports_ __728,__ 1–62
     (2018). [[DOI](http://dx.doi.org/10.1016/j.physrep.2017.11.003)]


