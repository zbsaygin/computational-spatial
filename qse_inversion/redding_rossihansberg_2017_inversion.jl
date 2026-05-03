#=
Model inversion for the canonical Quantitative Spatial Economics (QSE) model.

Given observed population L_n and wages w_n on a grid of locations, recovers
productivity A_n and housing supply H_n consistent with the general
equilibrium of the multi-region monopolistic-competition model.

Reference equations and notation follow:

    Redding, S. J. & Rossi-Hansberg, E. (2017),
    "Quantitative Spatial Economics," Annual Review of Economics 9: 21-58.
    Section 3.5-3.6, equations (16)-(17), p. 37.

Underlying model:        Helpman (1998), "The Size of Regions."
System reduction:        Allen & Arkolakis (2014, QJE).
Existence/uniqueness:    Fujimoto & Krause (1985).

Algorithm:
  1. Use eq (17) to write H = H(A; w, L) in closed form.
  2. Substitute into eq (16) to get a single nonlinear system in A.
  3. Solve by Newton's method with autodiff (NonlinearSolve.jl), using a
     log-transform A = exp(x) to enforce positivity.
  4. Verify the residual ~0 and plot heatmaps of {A, H, L, w}.
=#

using NonlinearSolve   # Newton, Trust-Region, etc., with a uniform AD interface
using ADTypes          # AutoForwardDiff
using LinearAlgebra
using Random
using Plots

Random.seed!(1)

## --- Output paths (figures land next to this script) ------------------------
const inversion_fig_dir = joinpath(@__DIR__, "figures", "inversion")
const data_fig_dir      = joinpath(@__DIR__, "figures", "data")
mkpath(inversion_fig_dir)
mkpath(data_fig_dir)

## --- Demo flag --------------------------------------------------------------
## Set to true to run the divergent inversion-by-iteration block at the bottom
## (50 iterations of a non-contractive fixed-point map; see the aside below).
const RUN_ITERATION_DEMO = false

## --- Build parameter / data bundle ------------------------------------------
## Returned NamedTuple `p` carries everything the equation routines need; it is
## passed through `NonlinearProblem(..., p)` rather than referenced as globals.
function build_params(; T::Int = 15, α = 0.75, σ = 5, F = 0.1, Wbar = 1.0)
    N = T^2

    coords_vec = vec(tuple.(collect(1:T)', reverse(collect(1:T))))
    d = zeros(N, N)
    for i in 1:N, j in 1:N
        x1, y1 = coords_vec[i]
        x2, y2 = coords_vec[j]
        d[i, j] = sqrt((x1 - x2)^2 + (y1 - y2)^2)
    end
    d = 1 .+ d ./ maximum(d)              # iceberg costs in [1, 2]

    w = rand(N) .+ 0.5                    # synthetic data — replace with real (w, L)
    L = rand(N) .+ 0.5

    γ1      = σ * (1 - α) / α
    γ2      = 1 + σ / (σ - 1) - (σ - 1) * (1 - α) / α
    σ_tilde = (σ - 1) / (2σ - 1)

    return (; T, N, α, σ, F, Wbar, γ1, γ2, σ_tilde, w, L, d)
end

const p = build_params()

## --- Equation (17) solved for H ---------------------------------------------
## Paper writes (17) as: w_n^(1-2σ) A_n^(σ-1) L_n^((σ-1)(1-α)/α) H_n^(-(σ-1)(1-α)/α) = ξ.
## Setting ξ = 1 and solving for H gives the closed form below.
function eq17_H(A, p)
    (; σ, α, w, L) = p
    return (w.^(1 - 2σ) .* A.^(σ - 1) .* L.^((σ - 1) * (1 - α) / α)).^(α / ((1 - α) * (σ - 1)))
end

## --- Equation (16) residual (target = 0) ------------------------------------
## After substituting H from eq (17), the equilibrium condition reduces to a
## single equation in A. We return RHS - LHS so a root-finder targets zero.
function eq16(A, p)
    (; σ, α, F, Wbar, γ1, γ2, σ_tilde, L, d) = p
    H   = eq17_H(A, p)
    lhs = L.^(σ_tilde * γ1) .*
          A.^(-(σ - 1)^2 / (2σ - 1)) .*
          H.^(-σ * (σ - 1) * (1 - α) / (α * (2σ - 1)))
    inside = (σ * F)^(-1) .* (σ / (σ - 1))^(1 - σ) .* d .*
             (L.^(σ_tilde * γ2) .*
              A.^(σ * (σ - 1) / (2σ - 1)) .*
              H.^((σ - 1)^2 * (1 - α) / (α * (2σ - 1))))'
    rhs = Wbar .* sum(inside, dims = 2)
    return rhs - lhs
end

## Wrap with log-transform A = exp(x) to enforce positivity in the unconstrained solve.
diff_eq16(x, p) = eq16(exp.(x), p)

## --- Solve for A ------------------------------------------------------------
## Solver choice matters here. The system is highly nonlinear (powers up to A^4,
## H^3 enter eq 16), so pure Newton's local quadratic model from a random
## initial guess is poor — `NewtonRaphson` returns retcode `:Unstable` with
## residual O(10^3). Methods with globalization succeed:
##   TrustRegion       residual ~ 1e-14   ✓ (used here)
##   PseudoTransient   residual ~ 1e-14   ✓
##   NewtonRaphson     residual ~ 1e+3    ✗ (no step control)
##   Broyden           residual ~ 1e+1    ✗ (Jacobian approximation insufficient)
x0   = log.(p.L)                          # initial guess: A_n = L_n
prob = NonlinearProblem(diff_eq16, x0, p)
@time sol = solve(prob, TrustRegion(; autodiff = AutoForwardDiff()))

A_solved = exp.(sol.u)
H_solved = eq17_H(A_solved, p)

## Self-check: equilibrium residual should be at the level of solver tolerance.
@assert sol.retcode == ReturnCode.Success "solver did not converge: retcode=$(sol.retcode)"
@show maximum(abs, eq16(A_solved, p))

## --- Plot heatmaps ----------------------------------------------------------
const T = p.T
A_grid = reshape(A_solved, T, T) ./ maximum(A_solved)
H_grid = reshape(H_solved, T, T) ./ maximum(H_solved)
L_grid = reshape(p.L,      T, T) ./ maximum(p.L)
w_grid = reshape(p.w,      T, T) ./ maximum(p.w)

function plot_grid(M, title_str; T::Int = size(M, 1))
    return heatmap(1:T, 1:T, M;
        xlabel = "x", ylabel = "y", title = title_str,
        color = :viridis, aspect_ratio = :equal,
        xticks = 1:T, yticks = 1:T,
        xlims = (0.5, T + 0.5), ylims = (0.5, T + 0.5))
end

pA = plot_grid(A_grid, "Model Inversion: A")
pH = plot_grid(H_grid, "Model Inversion: H")
pL = plot_grid(L_grid, "Data: L")
pw = plot_grid(w_grid, "Data: w")

display(pA); savefig(pA, joinpath(inversion_fig_dir, "A_grid.png"))
display(pH); savefig(pH, joinpath(inversion_fig_dir, "H_grid.png"))
display(pL); savefig(pL, joinpath(data_fig_dir,      "L_grid.png"))
display(pw); savefig(pw, joinpath(data_fig_dir,      "w_grid.png"))

## --- Aside: can the inversion be done by simple iteration? ------------------
## Rearranging eq (16) for A yields a fixed-point map T(A). It is NOT a
## contraction. Linearizing in log A, the Jacobian J of T satisfies
##
##     J · 𝟙 = λ₁ · 𝟙       with    λ₁ = -(3σ - 1) / (σ - 1)
##
## (the all-ones vector is a uniform-direction eigenvector). For σ > 1,
## |λ₁| = (3σ-1)/(σ-1) > 1 always (= 3.5 for σ = 5), so the equilibrium is
## locally repelling under T and iteration diverges with sign-flipping. RRH
## (2017) §3.5 establishes a contraction result for the *equilibrium* problem
## (solve for L given A, H, when γ₂/γ₁ ∈ (0,1)) — that does not apply to the
## inverse problem. The trust-region solver (above) sidesteps this entirely
## by globalizing each step.
function iterate_A_eq16(A, p)
    (; σ, α, F, Wbar, γ1, γ2, σ_tilde, L, d) = p
    H = eq17_H(A, p)
    inside = (σ * F)^(-1) .* (σ / (σ - 1))^(1 - σ) .* d .*
             (L.^(σ_tilde * γ2) .*
              A.^(σ * (σ - 1) / (2σ - 1)) .*
              H.^((σ - 1)^2 * (1 - α) / (α * (2σ - 1))))'
    rhs = Wbar .* sum(inside, dims = 2)
    lhs_part = L.^(σ_tilde * γ1) .*
               H.^(-σ * (σ - 1) * (1 - α) / (α * (2σ - 1)))
    return (rhs ./ lhs_part).^((2σ - 1) / (-(σ - 1)^2))
end

if RUN_ITERATION_DEMO
    let
        A = p.L ./ sum(p.L)
        tol = 1e-6
        max_iter = 50
        iter = 1
        err = 1.0
        while iter < max_iter && err > tol
            A_new = vec(iterate_A_eq16(A, p))
            err = (A_new - A)' * (A_new - A)
            println("iteration $iter, ‖ΔA‖² = $err")
            A = A_new
            iter += 1
        end
    end
end
