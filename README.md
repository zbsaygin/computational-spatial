# Computational Spatial Economics

Slides and code from a graduate spatial-economics course at the University of Chicago (2025–2026), where I TA'd. The slides cover the modern quantitative-spatial-economics workflow end to end. The code implements the algorithms from the early decks in Julia — more will follow as I find time.

## Slides

A series of five decks, intended to be read in order. They walk from the basic model-inversion problem up through the frontier methods for solving high-dimensional heterogeneous-agent dynamics.

| # | Deck | Topic |
|---|---|---|
| 01 | [Model Inversion](./slides/01_model_inversion.pdf) | The canonical Helpman (1998) / Allen–Arkolakis (2014) / Redding–Rossi-Hansberg (2017) quantitative spatial model, set up as an inversion problem. |
| 02 | [Counterfactuals](./slides/02_counterfactuals.pdf) | Caliendo, Parro, Rossi-Hansberg & Sarte (2018) — multi-region, multi-sector trade and migration via Eaton–Kortum (2002) Fréchet productivity. Hat algebra in the changes. |
| 03 | [Numerical Methods](./slides/03_numerical_methods.pdf) | Solvers (Newton–Raphson, Broyden, Levenberg–Marquardt, Nelder–Mead), fixed-point theory (Banach, Perov, Blackwell), Anderson acceleration, dampening, Walrasian tâtonnement, shooting and BVP methods. |
| 04 | [Dynamic Spatial Models](./slides/04_dynamic_models.pdf) | Caliendo, Dvorkin & Parro (2019) — dynamic discrete choice for migration with forward-looking workers; dynamic hat algebra; the China-shock counterfactual. |
| 05 | [Master Equation and FAME](./slides/05_master_equation.pdf) | Bilal (2023) and Bilal & Rossi-Hansberg (2023). Fréchet derivatives and adjoints, the master equation, the First-order Approximation to the Master Equation (FAME), and an application to a workers + capitalists + capital + climate model. |

## Code

| Folder | Companion deck | What it does |
|---|---|---|
| [`rrh2017_inversion/`](./rrh2017_inversion) | 01 | Solves equations (16)–(17) of RRH (2017) for productivity $A$ and housing supply $H$ using `NonlinearSolve.jl`. Includes a derivation of why naive iteration on equation (16) diverges — spectral-radius argument, $\lambda_1 = -(3\sigma-1)/(\sigma-1) > 1$ — and which solvers actually converge. |

Each code subfolder ships with its own pinned environment:

```bash
cd rrh2017_inversion
julia --project=. redding_rossihansberg_2017_inversion.jl
```

## References

- Allen, T. & Arkolakis, C. (2014). "Trade and the Topography of the Spatial Economy." *Quarterly Journal of Economics* 129(3): 1085–1140.
- Bilal, A. (2023). "Solving Heterogeneous-Agent Models with the Master Equation." Working paper.
- Bilal, A. & Rossi-Hansberg, E. (2023). "Anticipating Climate Change Across the United States." Working paper.
- Caliendo, L. & Parro, F. (2015). "Estimates of the Trade and Welfare Effects of NAFTA." *Review of Economic Studies* 82(1): 1–44.
- Caliendo, L., Dvorkin, M. & Parro, F. (2019). "Trade and Labor Market Dynamics: General Equilibrium Analysis of the China Trade Shock." *Econometrica* 87(3): 741–835.
- Caliendo, L., Parro, F., Rossi-Hansberg, E. & Sarte, P.-D. (2018). "The Impact of Regional and Sectoral Productivity Changes on the U.S. Economy." *Review of Economic Studies* 85(4): 2042–2096.
- Dornbusch, R., Fischer, S. & Samuelson, P. (1977). "Comparative Advantage, Trade, and Payments in a Ricardian Model with a Continuum of Goods." *American Economic Review* 67(5): 823–839.
- Eaton, J. & Kortum, S. (2002). "Technology, Geography, and Trade." *Econometrica* 70(5): 1741–1779.
- Helpman, E. (1998). "The Size of Regions." In *Topics in Public Economics*, Cambridge University Press.
- Redding, S. J. & Rossi-Hansberg, E. (2017). "Quantitative Spatial Economics." *Annual Review of Economics* 9: 21–58.
