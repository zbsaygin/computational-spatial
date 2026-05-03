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

## Gotcha 1: you can't solve eq (16) by iteration

The first thing you'd try is to rearrange eq (16) for $A$, treat the rearrangement as a map $T(A)$, and just iterate. **It diverges.** And it's not a fluke — there's a clean reason.

Linearize $T$ in $\log A$. The Jacobian $J$ has the all-ones vector $\mathbf{1}$ as an eigenvector with eigenvalue

$$ \lambda_1 \;=\; -\frac{3\sigma - 1}{\sigma - 1}. $$

For $\sigma > 1$, $|\lambda_1| > 1$ — always. (At the default $\sigma = 5$, $|\lambda_1| = 3.5$.) So equilibria are *locally repelling* under $T$, and the iteration explodes — with sign-flipping, since $\lambda_1 < 0$.

Read RRH carefully: they cite a contraction result from Allen & Arkolakis (2014) and Fujimoto & Krause (1985), but it's a contraction for the **forward** problem of solving for $L$ given $A, H$. The inverse problem we're doing here is a different map; the contraction result doesn't transfer.

The script ships this as a demo behind a flag (`RUN_ITERATION_DEMO`, default `false`). Flip it to `true` and you'll see 50 lines of $\|\Delta A\|^2$ blowing up.

## Gotcha 2: pure Newton isn't robust either

So switch to Newton — it has quadratic convergence near a root and doesn't care about contraction properties. Done, right?

Not quite. The bracketed terms in eq (16) carry powers up to $A^{\sigma-1} = A^4$ and $H^{(\sigma-1)^2(1-\alpha)/(\alpha(2\sigma-1))} \approx H^{1.78}$. From the initial guess $x_0 = \log L$, pure Newton's local quadratic model is way too aggressive — the solver bails with `retcode = :Unstable` and a residual on the order of $10^3$.

I tried four solvers from `NonlinearSolve.jl` on the same problem:

| Solver | Outcome | $\max_n \lvert (\text{residual})_n \rvert$ |
|---|---|---|
| `TrustRegion(autodiff=AutoForwardDiff())` | converges | $\sim 10^{-14}$ |
| `PseudoTransient(autodiff=AutoForwardDiff())` | converges | $\sim 10^{-14}$ |
| `NewtonRaphson(autodiff=AutoForwardDiff())` | `:Unstable` | $\sim 10^{3}$ |
| `Broyden(autodiff=AutoForwardDiff())` | `:Unstable` | $\sim 10^{1}$ |

The script defaults to `TrustRegion`. The trust-region globalization adapts each step's length based on actual-vs-predicted reduction, which is what saves you when the pure quadratic model is misleading.

(Embarrassing footnote: my first version of this script just ran `NewtonRaphson()` and used the result without checking `sol.retcode`. The heatmaps looked plausible because the input data was random — random in, random out. Adding `@assert sol.retcode == ReturnCode.Success` and an `@show maximum(abs, eq16(A, p))` is what finally surfaced this.)

## Running

```bash
julia --project=. redding_rossihansberg_2017_inversion.jl
```

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

The "data" here is `rand(N) .+ 0.5` — synthetic, just so the script runs end-to-end. The heatmaps look like noise because the inputs are noise.

## References

- Allen, T. & Arkolakis, C. (2014). "Trade and the Topography of the Spatial Economy." *QJE* 129(3): 1085–1140.
- Fujimoto, T. & Krause, U. (1985). "Strong Ergodicity for Strictly Increasing Nonlinear Operators." *J. Math. Econ.* 14(2): 119–125.
- Helpman, E. (1998). "The Size of Regions." In *Topics in Public Economics*.
- Redding, S. J. & Rossi-Hansberg, E. (2017). "Quantitative Spatial Economics." *Annu. Rev. Econ.* 9: 21–58. §3.5–3.6, eqs (16)–(17), p. 37.
