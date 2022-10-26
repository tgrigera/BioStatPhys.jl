# Correlation functions

The functions in this section compute estimates of several space and time correlation functions.  The precised definitions of the different quantities computed are summarised below.  For more details, caveats and discussion of the correlation functions we refer to the review article by T. S. Grigera[^review].



[^review]:

    The correlation functions estimators computed by the routines in this section are discussed in:

    - T. S. Grigera, Correlation functions as a tool to study collective behaviour phenomena in biological systems. _J. Phys. Complex._ __2,__ 045016 (2021).  [[DOI](https://doi.org/10.1088/2632-072X/ac2b06)]


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

Refer to the review article[^review] for more details.

### Correlation time

### API

The stationary (or TTI) estimate for a single real or complex sequence is implemented with an algorithm that uses fast-Fourier-transforms, giving ``O(N\log N)`` performance.

```@docs
time_correlation
```


## Space correlations and correlation length

The space correlation of a space-dependent quantity ``a(\mathbf{r])`` are defined as
```math
\begin{align*}
   C(\mathbf{r}) & = \left\langle a(\mathbf{r_0}) a(\mathbf{r_0}+\mathbf{r}) \right\rangle,\\
   C_c(\mathbf{r}) & = \left\langle \delta a(\mathbf{r_0}) \delta a(\mathbf{r_0}+\mathbf{r}) \right\rangle,
\end{align*}
```
where ``\delta a(\mathbf{r}) = a(\mathbf{r}) - \langle a\rangle`` and we are assuming the system is homogeneous (i.e. translation-invariant).  ``C_c(\mathbf{r})`` is called _connected_ correlation in the physics literature, or _auto covariance_ in the mathematical statistics literature.

The `space_correlation`(@ref) function computes the estimate of ``C(r)`` for an isotropic discrete system,
```math
\begin{align*}
  \hat C_c(r) & = \frac{\sum_{ij} \delta a_i \delta a_j
    \delta(r-r_{ij}) } {\sum_{kl} \delta(r-r_{kl})}.
\end{align*}
```
The average ``\langle a\rangle`` can be estimated with space average or phase average (see below and the review[^review] for details).

### API

```@docs
space_correlation
```
