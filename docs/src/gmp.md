# Generalized Moment Problems

The Generalized Moment Problem (GMP) as described in Lasserre (2009) is a linear optimization problem with countably many constraints in the space of measures. Formally, it takes the form
```math
\begin{array}{c l}
\min_{\color{red}\mu_i} & \sum_i \int {\color{blue}f_i}(x) d{\color{red}\mu_i}(x)\\
& \sum_i \int {\color{blue}h_i^\gamma}(x)d{\color{red}\mu_i}(x) \leq {\color{green}b^\gamma}, \quad \forall \gamma \in\Gamma\\
& \text{supp}({\color{red}\mu_i})\subseteq \{x\in\mathbb{R}^n | {\color{blue}g_i^1}(x)\geq 0,\ldots,{\color{blue}g_i^{m_i}}(x)\geq 0 \}
\end{array}
```
where ``{\color{red}\mu_i}`` are unknown (finite, positive, Borel) measures, ``\Gamma`` is a countable set, and ``{\color{green}b^\gamma}`` are given numbers.
In order to be able to compute approximations to this problem we require ``{\color{blue}f_i}``, ``{\color{blue}h_i^\gamma}``, and ``{\color{blue}g_i^j}`` to be polynomials.

Under some conditions a measure ``{\color{red}\mu}`` supported on some set ``{\color{blue}K}\subseteq \mathbb{R}^n`` is uniquely determined by its sequence of moments, defined by
```math
{\color{red}z_\alpha} = \int_{\color{blue}K} x^\alpha d{\color{red} \mu}(x)
```
where ``x^\alpha := \prod_i x_i^{\alpha_i}``. This is for example the case when ``{\color{blue}K}`` is compact, or when the sequence ``{\color{red}(z_\alpha)_\alpha}`` satisfies the so-called Carleman's condition. Therefore one is often interested in the sequences of moments, rather than the measure itself, motivating the name generalized MOMENT problem.

## Example: Polynomial Optimization

Given a polynomial ``{\color{blue}f}`` and a set ``{\color{blue}K} = \{x\in\mathbb{R}^n | {\color{blue}g^1}(x)\geq 0,\ldots,{\color{blue}g^{m}}(x)\geq 0 \}`` one is interested in the global minimum
```math
\begin{array}{c l}
\min_{\color{red}x} &  {\color{blue}f}({\color{red}x})\\
& {\color{red}x}\in{\color{blue}K}.
\end{array}
```
This problem has a GMP formulation in terms, where we optimize over one single measure  ``{\color{red}\mu}`` and ``\#\Gamma=1``:
```math
\begin{array}{c l}
\min_{\color{red}\mu} &  \int {\color{blue}f}(x) d{\color{red}\mu}(x)\\
& \int {\color{blue}1}d{\color{red}\mu}(x) = {\color{green}1},\\
& \text{supp}({\color{red}\mu})\subseteq{\color{blue}K}
\end{array}
```
Note that instead of optimizing over ``{\color{red}x}\in{\color{blue}K}``, the variable of the GMP is the measure ``{\color{red}\mu}``, hence turning a non-linear problem into a linear one.

The equivalence of both problems can be seen as follows: First note, that for every ``x\in{\color{blue}K}``, the Dirac measure ``\delta_{x}`` is a probability measure supported in ``{\color{blue}K}`` and hence feasible for the GMP formulation. Together with the equality ``{\color{blue}f}(x) = \int {\color{blue}f}d\delta_x`` this shows that the optimal value of the GMP formulation is always smaller than the one of the polynomial optimization problem.

For the other inequality let ``{\color{red}f^\ast}`` denote the global minimum of ``{\color{blue}f}`` on ``{\color{blue}K}`` and let ``{\color{red}\mu^\ast}`` be the optimal measure for the GMP formulation. Then
```math
\int {\color{blue}f} d{\color{red}\mu^\ast}\geq \int{\color{red}f^\ast} d{\color{red}\mu^\ast}= {\color{red}f^\ast}\int {\color{blue}1} d{\color{red}\mu^\ast} = {\color{red}f^\ast} {\color{green}1} = {\color{red}f^\ast},
```
showing that the optimal value of the GMP formulation is always bigger than the one of the polynomial optimization problem - and in combination with the first part, equivalence of both problems.

In MomentOpt.jl this GMP can be modeled as demonstrated in [polopt](https://github.com/lanl-ansi/MomentOpt.jl/blob/master/examples/polopt.jl) in the example folder.

