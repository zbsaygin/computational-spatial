# QSE Inversion

A Julia implementation of the model-inversion procedure that Redding & Rossi-Hansberg (2017) describe in §3.5–3.6 of their *Annual Review* survey, building on Helpman (1998) and the existence/uniqueness machinery from Allen & Arkolakis (2014). I wrote it for a TA section on spatial economics at Chicago.

The companion slide deck covers the underlying theory: [`../slides/01_model_inversion.pdf`](../slides/01_model_inversion.pdf).

## The setup

There's a grid of locations. At each one we observe population $L_n$ and wages $w_n$. We want to recover the *unobservables* — productivity $A_n$ and housing supply $H_n$ — consistent with the general equilibrium of the model (monopolistic competition, iceberg trade costs $d_{ni}$, Cobb–Douglas preferences over a goods aggregate and residential land).

That's "inverting the model," and RRH §3.6 walks through the recipe.

## The recipe

Equation (17) of the paper ties the four quantities together:

$$ w_n^{1-2\sigma}\, A_n^{\sigma-1}\, L_n^{(\sigma-1)(1-\alpha)/\alpha}\, H_n^{-(\sigma-1)(1-\alpha)/\alpha} = \xi. $$

Pick the numeraire $\xi = 1$ and solve for $H$ — that gives a closed form $H = H(A; w, L)$. So once we know $A$, we know $H$. Substitute into equation (16), the labor-market / market-clearing reduced form, and you're left with a single nonlinear system in $A$ alone.

That's the system the script solves with `NonlinearSolve.jl`. I work in log-coordinates $x = \log A$, which keeps $A > 0$ for free without needing constrained optimization.

## Attempt 1: fixed-point iteration

The first thing you'd try is to rearrange eq (16) for $A$, treat the rearrangement as a map $T(A)$, and iterate. **It diverges.** $T$'s Jacobian has spectral radius greater than 1 for any $\sigma > 1$, so equilibria are locally repelling. The contraction theory and the explicit spectral-radius derivation are in [`../slides/03_numerical_methods.pdf`](../slides/03_numerical_methods.pdf).

RRH do cite a contraction result, but it's for the *forward* problem (solve for $L$ given $A, H$, when $\gamma_2/\gamma_1 \in (0,1)$). The inverse problem is a different map; that result doesn't transfer.

The script ships this as a demo behind a flag (`RUN_ITERATION_DEMO`, default `false`). Flip it to `true` and you'll see 50 lines of $\|\Delta A\|^2$ blowing up.

## Attempt 2: pure Newton

So switch to Newton — it has quadratic convergence near a root and doesn't care about contraction properties. Done, right?

Not quite. The bracketed terms in eq (16) carry powers up to $A^{\sigma-1} = A^4$ and $H^{(\sigma-1)^2(1-\alpha)/(\alpha(2\sigma-1))} \approx H^{1.78}$. From the initial guess $x_0 = \log L$, pure Newton's local quadratic model is too aggressive — the solver fails with `retcode = :Unstable` and a residual on the order of $10^3$.

I tried four solvers from `NonlinearSolve.jl` on the same problem:

| Solver | Outcome | $\max_n \lvert (\text{residual})_n \rvert$ |
|---|---|---|
| `TrustRegion(autodiff=AutoForwardDiff())` | converges | $\sim 10^{-14}$ |
| `PseudoTransient(autodiff=AutoForwardDiff())` | converges | $\sim 10^{-14}$ |
| `NewtonRaphson(autodiff=AutoForwardDiff())` | `:Unstable` | $\sim 10^{3}$ |
| `Broyden(autodiff=AutoForwardDiff())` | `:Unstable` | $\sim 10^{1}$ |

The script defaults to `TrustRegion`. Trust-region globalization adapts each step's length based on actual-vs-predicted reduction, which is what saves you when the pure quadratic model is misleading. Deck 03 ([`../slides/03_numerical_methods.pdf`](../slides/03_numerical_methods.pdf)) covers the broader theory of solvers and globalization strategies.

## Outputs

- `figures/inversion/A_grid.png`, `figures/inversion/H_grid.png` — recovered productivity and housing supply
- `figures/data/L_grid.png`, `figures/data/w_grid.png` — input population and wages

## Default parameters

| Symbol | Value | What it is |
|---|---|---|
| $\alpha$ | 0.75 | Cobb–Douglas share on tradables |
| $\sigma$ | 5 | CES elasticity of substitution |
| $F$ | 0.1 | Fixed labor cost per variety |
| $\bar W$ | 1.0 | Reservation utility (numeraire) |
| $T$ | 15 | Grid side; $T^2 = 225$ locations |

To override, edit the keyword args to `build_params(...)`. Note that RRH show the inversion is unique (up to normalization) only when $\sigma(1 - \alpha) > 1$. The defaults give $\sigma(1 - \alpha) = 1.25$.

## A note on the data

The "data" here is `rand(N) .+ 0.5` — synthetic; the script is a self-contained demo.

## References

- Allen, T. & Arkolakis, C. (2014). "Trade and the Topography of the Spatial Economy." *QJE* 129(3): 1085–1140.
- Fujimoto, T. & Krause, U. (1985). "Strong Ergodicity for Strictly Increasing Nonlinear Operators." *J. Math. Econ.* 14(2): 119–125.
- Helpman, E. (1998). "The Size of Regions." In *Topics in Public Economics*.
- Redding, S. J. & Rossi-Hansberg, E. (2017). "Quantitative Spatial Economics." *Annu. Rev. Econ.* 9: 21–58. §3.5–3.6, eqs (16)–(17), p. 37.
